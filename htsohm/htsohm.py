import os
import sys
from math import sqrt
from datetime import datetime

import numpy as np
from sqlalchemy.sql import func, or_
from sqlalchemy.orm.exc import FlushError
from sqlalchemy.exc import IntegrityError
from sqlalchemy.sql import text
import yaml

import htsohm
from htsohm import config
from htsohm.db import engine, session, Material, MutationStrength, Structure
from htsohm.material_files import generate_material, mutate_material
from htsohm import simulation

def materials_in_generation(run_id, generation):
    """Count number of materials in a generation.

    Args:
        run_id (str): identification string for run.
        generation (int): iteration in overall bin-mutate-simulate rountine.

    Returns:
        Number(int) of materials in a particular generation that are present in
        the database (the final step in bin-mutate-simulate routine).

    """
    return session.query(Material).filter(
        Material.run_id == run_id,
        Material.generation == generation
    ).count()

def last_generation(run_id):
    """Finds latest generation present in database.

    Args:
        run_id (str): identification string for run.

    Returns:
        Last generation(int) to be included in database.

    """
    return session.query(func.max(Material.generation)).filter(
        Material.run_id == run_id,
    )[0][0]

def calc_bin(value, bound_min, bound_max, bins):
    """Find bin in parameter range.

    Args:
        value (float): some value, the result of a simulation.
        bound_min (float): lower limit, defining the parameter-space.
        bound_max (float): upper limit, defining the parameter-space.
        bins (int): number of bins used to subdivide parameter-space.

    Returns:
        Bin(int) corresponding to the input-value.

    """
    step = (bound_max - bound_min) / bins
    assigned_bin = (value - bound_min) // step
    assigned_bin = min(assigned_bin, bins-1)
    assigned_bin = max(assigned_bin, 0)
    return int(assigned_bin)

def select_parent(run_id, max_generation, generation_limit):
    """Use bin-counts to preferentially select a list of 'rare' parents.

    Args:
        run_id (str): identification string for run.
        max_generation (int): latest generation to include when counting number
            of materials ub each bin.
        generation_limit (int): number of materials to query in each generation
            (as materials are added to database they are assigned an index
            within the generation to bound the number of materials in each
            generation).

    Returns:
        The material id(int) corresponding to some parent-material selected
        from database with a bias favoring materials in bins with the lowest
        counts.

    """
    # Determine which dimensions are being binned
    simulations = config['material_properties']
    queries = []
    if 'gas_adsorption_0' in simulations or 'gas_adsorption_1' in simulations:
        queries.append( getattr(Material, 'gas_adsorption_bin') )
    if 'surface_area' in simulations:
        queries.append( getattr(Material, 'surface_area_bin') )
    if 'helium_void_fraction' in simulations:
        queries.append( getattr(Material, 'void_fraction_bin') )

    # Count the number of materials in each bin
    bins_and_counts = session \
        .query(
            func.count(Material.id),
            *queries
        ) \
        .filter(
            Material.run_id == run_id,
            or_(Material.retest_passed == True, Material.retest_passed == None),
            Material.generation <= max_generation,
            Material.generation_index < generation_limit,
        ) \
        .group_by(*queries).all()[1:]

    # Print bin-counts in interactive mode
    if config['selection_mode'] == 'manual':
        print('\n  BIN\t\t|  COUNT')
        print('----------------+-----------')
        for bin_and_count in bins_and_counts:
            print('  {}\t|  {}'.format(bin_and_count[1:len(queries) + 2], bin_and_count[0]))
        print('\n')
    
        input('Press Enter to calculate normalized weights...')

    # Extract only bin coordinate from bins_and_counts query
    bins = []
    for i in bins_and_counts:
        some_bin = {}
        for j in range(len(queries)):
            some_bin[queries[j]] = i[j + 1]
        bins.append(some_bin)

    # Assign weight to each bin
    total = sum([i[0] for i in bins_and_counts])
    weights = [ total / float(i[0]) for i in bins_and_counts ]
    normalized_weights = [ weight / sum(weights) for weight in weights ]

    # Print bins and weights in interactive mode
    if config['selection_mode'] == 'manual':
        print('\n  BIN\t\t|  COUNT\t|  WEIGHT\t|  NORMALIZED-WEIGHT')
        print('----------------+---------------+---------------+---------------------')
        for i in range(len(bins_and_counts)):
            print('  {}\t|  {}\t\t|  {}\t\t|  {}'.format(
                bins_and_counts[i][1:len(queries) + 1], bins_and_counts[i][0],
                weights[i], normalized_weights[i]))
        print('\n')
    
    # Have user select a parent-bin
        ga_bin = input('Press select a gas adsorption bin :\t')
        sa_bin = input('Press select a surface area bin :\t')
        vf_bin = input('Press select a void fraction bin :\t')
    
    # Query parent-ids corresponding to selected bin
        parent_query = session \
            .query(Material.uuid) \
            .filter(
                Material.run_id == run_id,
                or_(Material.retest_passed == True, Material.retest_passed == None),
                Material.gas_adsorption_bin == ga_bin,
                Material.surface_area_bin == sa_bin,
                Material.void_fraction_bin == vf_bin,
                Material.generation <= max_generation,
                Material.generation_index < generation_limit,
            ).all()
        potential_parents = [i[0] for i in parent_query]
    
        print('Potential parents :')
        for i in potential_parents:
            print('\t{}'.format(i))
    
    # Have user select parent-pseudomaterial
        parent_uuid = input('Please select a parent UUID :\t')
        parent_id = session.query(Material.id).filter(Material.uuid==parent_uuid).one()

    # Weighted-selection for parent bin (non-interactive mode)
    else:
        parent_bin = np.random.choice(bins, p = normalized_weights)
        parent_queries = [i == parent_bin[i] for i in queries]
        parent_query = session \
            .query(Material.id) \
            .filter(
                Material.run_id == run_id,
                or_(Material.retest_passed == True, Material.retest_passed == None),
                *parent_queries,
                Material.generation <= max_generation,
                Material.generation_index < generation_limit,
            ).all()
        potential_parents = [i[0] for i in parent_query]
    # Randomly-select parent-pseudomaterial (non-interactive mode)
        parent_id = int(np.random.choice(potential_parents))

    return parent_id

def run_all_simulations(material):
    """Simulate helium void fraction, gas loading, and surface area.

    Args:
        material (sqlalchemy.orm.query.Query): material to be analyzed.

    Depending on properties specified in config, adds simulated data for helium
    void fraction, gas loading, heat of adsorption, surface area, and
    corresponding bins to row in database corresponding to the input-material.
        
    """
    simulations = config['material_properties']

    print('\n===============================|')
    print('UUID :\t{}  |'.format(material.uuid))
    print('----------------+---------------+-------------|')
    print('  PROPERTY\t|  VALUE\t|  BIN\t      |')
    print('----------------+---------------+-------------|')

    ############################################################################
    # run helium void fraction simulation
    if 'helium_void_fraction' in simulations:
        if config['artificial_data'] == 'off':
            results = simulation.helium_void_fraction.run(
                material.run_id, material)
            material.update_from_dict(results)
            material.void_fraction_bin = calc_bin(
                material.vf_helium_void_fraction,
                *config['helium_void_fraction']['limits'],
                config['number_of_convergence_bins']
            )
        else:
            results = {}
            results['vf_helium_void_fraction'] = material.artificial_void_fraction()
            material.update_from_dict(results)
            material.void_fraction_bin = calc_bin(
                material.vf_helium_void_fraction,
                *config['helium_void_fraction']['limits'],
                config['number_of_convergence_bins']
            )
            print('void fraction\t|  {}\t|  {}\t      |'.format(
                round(results['vf_helium_void_fraction'], 4),
                material.void_fraction_bin))

    else:
        material.void_fraction_bin = 0
    ############################################################################
    # run gas loading simulation
    if 'gas_adsorption_0' in simulations and 'gas_adsorption_1' not in simulations:
        if config['artificial_data'] == 'off':
            arguments = [material.run_id, material]
            if 'helium_void_fraction' in simulations:
                arguments.append(material.vf_helium_void_fraction)
            results = simulation.gas_adsorption_0.run(*arguments)
            material.update_from_dict(results)
            material.gas_adsorption_bin = calc_bin(
                material.ga0_absolute_volumetric_loading,
                *config['gas_adsorption_0']['limits'],
                config['number_of_convergence_bins']
            )
        else:
            results = {}
            results['ga0_absolute_volumetric_loading'] = material.artificial_gas_adsorption()
            material.update_from_dict(results)
            material.gas_adsorption_bin = calc_bin(
                material.ga0_absolute_volumetric_loading,
                *config['gas_adsorption_0']['limits'],
                config['number_of_convergence_bins']
            )
            print('gas adsorption\t|  {}\t|  {}\t      |'.format(
                round(results['ga0_absolute_volumetric_loading'], 4),
                material.gas_adsorption_bin))

    elif 'gas_adsorption_0' in simulations and 'gas_adsorption_1' in simulations:
        if config['artificial_data'] == 'off':
            arguments = [material.run_id, material]
            if 'helium_void_fraction' in simulations:
                arguments.append(material.vf_helium_void_fraction)
            results = simulation.gas_adsorption_0.run(*arguments)
            material.update_from_dict(results)
            results = simulation.gas_adsorption_1.run(*arguments)
            material.update_from_dict(results)
            material.gas_adsorption_bin = calc_bin(
                abs(material.ga0_absolute_volumetric_loading - material.ga1_absolute_volumetric_loading),
                *config['gas_adsorption_0']['limits'],
                config['number_of_convergence_bins']
            )
        else:
            results = {}
            results['ga0_absolute_volumetric_loading'] = material.artificial_gas_adsorption()
            results['ga1_absolute_volumetric_loading'] = material.artificial_gas_adsorption()
            material.update_from_dict(results)
            material.gas_adsorption_bin = calc_bin(
                material.ga0_absolute_volumetric_loading,
                *config['gas_adsorption_0']['limits'],
                config['number_of_convergence_bins']
            )
            print('gas adsorption\t|  {}\t|  {}\t      |'.format(
                round(results['ga0_absolute_volumetric_loading'], 4),
                material.gas_adsorption_bin))

    elif 'gas_adsorption_0' not in simulations:
        material.gas_adsorption_bin = 0
    ############################################################################
    # run surface area simulation
    if 'surface_area' in simulations:
        if config['artificial_data'] == 'off':
            results = simulation.surface_area.run(
                    material.run_id, material)
            material.update_from_dict(results)
            material.surface_area_bin = calc_bin(
                material.sa_volumetric_surface_area,
                *config['surface_area']['limits'],
                config['number_of_convergence_bins']
            )
        else:
            results = {}
            results['sa_volumetric_surface_area'] = material.artificial_surface_area()
            material.update_from_dict(results)
            material.surface_area_bin = calc_bin(
                material.sa_volumetric_surface_area,
                *config['surface_area']['limits'],
                config['number_of_convergence_bins']
            )
            print('surface area\t|  {}\t|  {}\t      |'.format(
                round(results['sa_volumetric_surface_area'], 4),
                material.surface_area_bin))

    else:
        material.surface_area_bin = 0

def retest(m_orig, retests, tolerance):
    """Reproduce simulations  to prevent statistical errors.

    Args:
        m_orig (sqlalchemy.orm.query.Query): material to retest.
        retests (int): number of times to reproduce each simulation.
        tolerance (float): acceptance criteria as percent deviation from
            originally calculated value.

    Queries database to determine if there are remained retests to  be run.
    Updates row in database with total number of retests and results.

    """
    m = m_orig.clone()

    run_all_simulations(m)
    print('\n\nRETEST_NUM :\t%s' % m_orig.retest_num)
    print('retests :\t%s' % retests)

    simulations = config['material_properties']

    # requery row from database, in case someone else has changed it, and lock it
    # if the row is presently locked, this method blocks until the row lock is released
    session.refresh(m_orig, lockmode='update')
    if m_orig.retest_num < retests:
        if 'gas_adsorption_0' in simulations:
            m_orig.retest_gas_adsorption_0_sum += m.ga0_absolute_volumetric_loading
        if 'gas_adsorption_1' in simulations:
            m_orig.retest_gas_adsorption_1_sum += m.ga1_absolute_volumetric_loading
        if 'surface_area' in simulations:
            m_orig.retest_surface_area_sum += m.sa_volumetric_surface_area
        if 'helium_void_fraction' in simulations:
            m_orig.retest_void_fraction_sum += m.vf_helium_void_fraction
        m_orig.retest_num += 1
        session.commit()        

    if m_orig.retest_num >= retests:
        try:
            m_orig.retest_passed = m.calculate_retest_result(tolerance)
            print('\nRETEST_PASSED :\t%s' % m_orig.retest_passed)
            session.commit()
        except ZeroDivisionError as e:
            print('WARNING: ZeroDivisionError - material.calculate_retest_result(tolerance)')

def get_all_parent_ids(run_id, generation):
    return [e[0] for e in session.query(Material.parent_id) \
            .filter(Material.run_id == run_id, Material.generation == generation) \
            .distinct() if e[0] != None]

def calculate_percent_children_in_bin(run_id, generation, bin_coordinate):
    """Find the fraction of children generated by all parents from a particular
    bin who also belong in the same bin as their parent.

    Args:
        run_id (str): run identification string.
        generation (int): interation in overall process.
        bin (list int): [gas_adsorption_bin, surface_area_bin, void_fraction_bin].

    Returns:
        Fraction of children in the same bin as their parent for all parents in
        a particular bin.
    """
    sql = text("""
        select
            m.gas_adsorption_bin,
            m.surface_area_bin,
            m.void_fraction_bin,
            (
                m.gas_adsorption_bin = p.gas_adsorption_bin and
                m.surface_area_bin = p.surface_area_bin and
                m.void_fraction_bin = p.void_fraction_bin
            ) as in_bin
        from materials m
        join materials p on (m.parent_id = p.id)
        where m.generation = :gen
            and m.run_id = :run_id
            and p.gas_adsorption_bin = :ga_bin
            and p.surface_area_bin = :sa_bin
            and p.void_fraction_bin = :vf_bin
        """)

    rows = engine.connect().execute(
        sql,
        gen=generation,
        run_id=run_id,
        ga_bin=bin_coordinate[0],
        sa_bin=bin_coordinate[1],
        vf_bin=bin_coordinate[2]
    ).fetchall()

    if (config['interactive_mode'] == 'on' and 
            len([r for r in rows if r.in_bin]) != 0):
        print('Number of children in the parent-bin :\t{}' \
                .format(len([r for r in rows if r.in_bin])))
        print('Total number of children :\t\t{}'.format(len(rows)))

    return len([ r for r in rows if r.in_bin ]) / len(rows)

def double_check():
    while True:
        proceed = input('Are you sure? [y/n] :\t')
        if proceed in ['n', 'N']:
            return False
        elif proceed in ['y', 'Y']:
            return True
        else:
            print('Please enter y/n.')

def set_variable(message, variable_name):
    while True:
        value = input(message)
        print('Setting {} to {} ...'.format(variable_name, value))
        if double_check():
            return value

def calculate_mutation_strength(run_id, generation, mutation_strength_bin):
    """Query mutation_strength for bin and adjust as necessary.

    Args:
        run_id (str): identification string for run.
        generation (int): iteration in bin-mutate-simulate routine.
        parent (sqlalchemy.orm.query.Query): parent-material corresponding to
            the bin being queried.

    Returns:
        mutation_strength.strength (float): mutation strength to be used for
        parents in the bin being queried. If the fraction of children from
        previous generation which populate the SAME bin as their parent is LESS
        THAN 10% then the mutation strength is REDUCED BY 5%. If the fraction
        of these children populating the SAME bin as their parent is GREATER
        THAN 50% then the mutation strength is INCREASED BY 5%.

    """
    mutation_strength_key = [run_id, generation] + mutation_strength_bin
    mutation_strength = session.query(MutationStrength).get(mutation_strength_key)

    if mutation_strength:
        print("Mutation strength already calculated for this bin and generation.")
    else:
        print("Calculating mutation strength...")
        mutation_strength = MutationStrength.get_prior(*mutation_strength_key).clone()
        if config['interactive_mode'] == 'on':
            print('Prior mutation strength :\t{}'.format(mutation_strength.strength))
        mutation_strength.generation = generation

        try:
            fraction_in_parent_bin = calculate_percent_children_in_bin(run_id, generation, mutation_strength_bin)

            if config['interactive_mode'] != 'on':
                if fraction_in_parent_bin < 0.1 and mutation_strength.strength - 0.05 > 0:
                    mutation_strength.strength -= 0.05
                elif fraction_in_parent_bin > 0.5 and mutation_strength.strength + 0.05 < 1:
                    mutation_strength.strength += 0.05

            elif config['interactive_mode'] == 'on':
                print('\tFor adaptive mutation strength(s), DECREASE strength \n' +
                        '\tby 5% if the fraction of children in parent-bin is \n' +
                        '\tLESS THAN 0.1; INCREASE strength by 5% if the \n' +
                        '\tfraction of children in parent-bin is GREATER THAN \n' +
                        '\t0.5.')
                mutation_strength.strength = set_variable(
                    'Please Enter new mutation strength :\t', 'mutation_strength.strength')
        except ZeroDivisionError:
            print("No prior generation materials in this bin with children.")

        try:
            session.add(mutation_strength)
            session.commit()
        except (FlushError, IntegrityError) as e:
            print("Somebody beat us to saving a row with this generation. That's ok!")
            session.rollback()
            # it's ok b/c this calculation should always yield the exact same result!
    sys.stdout.flush()
    return mutation_strength.strength

def evaluate_convergence(run_id, generation):
    '''Determines convergence by calculating variance of bin-counts.
    
    Args:
        run_id (str): identification string for run.
        generation (int): iteration in bin-mutate-simulate routine.

    Returns:
        bool: True if variance is less than or equal to cutt-off criteria (so
            method will continue running).
    '''
    simulations = config['material_properties']
    query_group = []
    if 'gas_adsorption' in simulations:
        query_group.append( getattr(Material, 'gas_adsorption_bin') )
    if 'surface_area' in simulations:
        query_group.append( getattr(Material, 'surface_area_bin') )
    if 'helium_void_fraction' in simulations:
        query_group.append( getattr(Material, 'void_fraction_bin') )

    bin_counts = session \
        .query(func.count(Material.id)) \
        .filter(
            Material.run_id == run_id, Material.generation < generation,
            Material.generation_index < config['children_per_generation']
        ) \
        .group_by(*query_group).all()
    bin_counts = [i[0] for i in bin_counts]    # convert SQLAlchemy result to list
    variance = sqrt( sum([(i - (sum(bin_counts) / len(bin_counts)))**2 for i in bin_counts]) / len(bin_counts))
    print('\nCONVERGENCE:\t%s\n' % variance)
    sys.stdout.flush()
    return variance <= config['convergence_cutoff_criteria']

def print_block(string):
    print('{0}\n{1}\n{0}'.format('=' * 80, string))

def worker_run_loop(run_id):
    """
    Args:
        run_id (str): identification string for run.

    Writes seed generation and simulates properties, then manages overall
    bin-mutate-simualte routine until convergence cutt-off or maximum
    number of generations is reached.

    """
    print('CONFIG\n{0}'.format(config))

    gen = last_generation(run_id) or 0

    converged = False
    while not converged:
        print_block('GENERATION {}'.format(gen))
        size_of_generation = config['children_per_generation']

        while materials_in_generation(run_id, gen) < size_of_generation:
            if gen == 0:
                print("writing new seed...")
                material = generate_material(run_id, config['number_of_atom_types'])
            else:
                print_block('Selecting parent...')
                
                parent_id = select_parent(run_id, max_generation=(gen - 1),
                        generation_limit=config['children_per_generation'])

                parent_material = session.query(Material).get(parent_id)

                # run retests until we've run enough
                print_block('Running retest...')
                while parent_material.retest_passed is None:
                    print("running retest...")
                    print("Date :\t%s" % datetime.now().date().isoformat())
                    print("Time :\t%s" % datetime.now().time().isoformat())
                    retest(parent_material, config['retests']['number'], 
                            config['retests']['tolerance'])
                    session.refresh(parent_material)

                if not parent_material.retest_passed:
                    print("parent failed retest. restarting with parent selection.")
                    continue

                # obtain mutation strength
                if config['mutation_strength_method'] == 'flat':
                    mutation_strength = config['initial_mutation_strength']
                else:
                    mutation_strength_key = [run_id, gen] + parent_material.bin
                    mutation_strength = MutationStrength \
                            .get_prior(*mutation_strength_key).clone().strength
                
                # mutate material
                print_block('Mutating pseudomaterial..')
                material = mutate_material(parent_material, mutation_strength, gen)
            run_all_simulations(material)
            session.add(material)
            session.commit()

            material.generation_index = material.calculate_generation_index()
            if material.generation_index < config['children_per_generation']:
                print_block('ADDING MATERIAL {}'.format(material.uuid))
                session.add(material)

            if config['interactive_mode'] == 'on':
                print("\n\nGeneration :\t{}\nMat. no. :\t{} / {}".format(gen,
                    material.generation_index + 1, config['children_per_generation']))
                input("\nPress Enter to continue...\n")

            if config['mutation_strength_method'] != 'flat':
                if config['interactive_mode'] != 'on':
                    if material.generation_index == config['children_per_generation'] - 1 and gen > 0:
                    # standard calculation of mutation strengths for all accessed bins
                        if config['annealing_on'] != 'on' or gen % config['annealing_frequency'] != 0:
                            parent_ids = get_all_parent_ids(run_id, gen)
                            print_block('CALCULATING MUTATION STRENGTHS')
                            ms_bins = []
                            for parent_id in parent_ids:
                                parent_bin = session.query(Material).get(parent_id).bin
                                if parent_bin not in ms_bins:
                                    print('Calculating bin-mutation-strength for bin : {}' \
                                            .format(parent_bin))
                                    calculate_mutation_strength(run_id, gen + 1, parent_bin)
                                ms_bins.append(parent_material.bin)
                        # annealing to reset all mutation strengths to initial value
                        elif config['annealing_on'] == 'on' and gen % config['annealing_frequency'] == 0:
                            all_accessed_bin_tuples = session \
                                    .query(
                                        Material.gas_adsorption_bin,
                                        Material.surface_area_bin,
                                        Material.void_fraction_bin) \
                                    .filter(
                                        Material.run_id == run_id,
                                        Material.retest_passed != False,
                                        Material.generation_index < config['children_per_generation']) \
                                    .distinct().all()
                            all_accessed_bins = [ [e[0], e[1], e[2]]
                                    for e in all_accessed_bin_tuples]
                            print_block('ANNEALING WITH MUTATION STRENGTH :\t{}' \
                                    .format(config['initial_mutation_strength']))
                            for some_bin in all_accessed_bins:
                                print('Annealing bin :\t{}'.format(some_bin))
                                mutation_strength = MutationStrength(run_id, gen + 1, *some_bin)
                                mutation_strength.strength = config['annealing_strength']
                                session.add(mutation_strength)

                if config['interactive_mode'] == 'on':
                    if material.generation_index == config['children_per_generation'] - 1 and gen > 0:
                        parent_ids = get_all_parent_ids(run_id, gen)
                        ms_bins = []
                        for parent_id in parent_ids:
                            parent_bin = session.query(Material).get(parent_id).bin
                            if parent_bin not in ms_bins:
                                ms_bins.append(parent_bin)
                        print('The following bins contain parents for the preceding generation:')
                        for some_bin in ms_bins:
                            print('\t{}'.format(some_bin))
                        input('Press Enter to begin calculating new mutation strength(s) :')
                        counter = 1
                        for some_bin in ms_bins:
                            print('\nCalculating strength for bin :\t{}'.format(some_bin))
                            print('Bin {} / {}'.format(counter, len(ms_bins)))
                            calculate_mutation_strength(run_id, gen + 1, parent_bin)
                            input('Press Enter to continue :')
            else:
                # delete excess rows
                # session.delete(material)
                pass
            session.commit()
            sys.stdout.flush()
        gen += 1
        converged = evaluate_convergence(run_id, gen)
