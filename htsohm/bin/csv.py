#!/usr/bin/env python3
import csv
import sys

import click
from sqlalchemy.orm import joinedload

from htsohm import db
from htsohm.db import Material

@click.command()
@click.argument('database-path', type=click.Path())
def output_csv(database_path):
    db.init_database(db.get_sqlite_dbcs(database_path))
    session = db.get_session()

    mats = session.query(Material) \
        .options(joinedload("structure").joinedload("lennard_jones")) \
        .options(joinedload("gas_loading")) \
        .options(joinedload("void_fraction")).all()

    f = csv.writer(sys.stdout, lineterminator="\n")
    f.writerow(["id", "a", "b", "c", "sigma", "epsilon", "void_fraction", "absolute_volumetric_loading", "absolute_volumetric_loading_error"])
    for m in mats:
        f.writerow([m.id, m.structure.a, m.structure.b, m.structure.c,
            m.structure.lennard_jones[0].sigma, m.structure.lennard_jones[0].epsilon,
            m.void_fraction[0].void_fraction,
            m.gas_loading[0].absolute_volumetric_loading, m.gas_loading[0].absolute_volumetric_loading_error
        ])

if __name__ == '__main__':
    output_csv()
