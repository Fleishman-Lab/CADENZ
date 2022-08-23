#!/usr/bin/env python3

"""
Parsing and manipulation of a torsion database file, outut from AbDEsign
"""

__author__ = 'Rosalie Lipsh'

import os
import gzip
import tempfile
import pandas as pd
from typing import List
from Bio.SeqUtils import seq1, seq3
from Bio.PDB.PDBParser import PDBParser

def parse_pdb_file(pdb_path):
    """Reads a pdb file and return a structure object"""
    if isinstance(pdb_path, Bio.PDB.Structure.Structure):
        return pdb_path
    if pdb_path.endswith('.gz'):
        with gzip.open(pdb_path, 'rb') as f:
            file_content = f.read()
            with tempfile.NamedTemporaryFile() as tmphandler:
                tmphandler.write(file_content)
                tmphandler.seek(0)
                pdb = PDBParser().get_structure('', tmphandler.name)
    else:
        pdb = PDBParser().get_structure('',pdb_path)
    return pdb

class SpliceDB:
    def __init__(self, db):
        self.original_path = db if os.path.isfile(db) else ''
        self.__set_db(self.original_path)

    def __set_db(self, path):
        """Parses a torsions db file into a pandas table"""
        db = open(path).readlines()
        if len(db) != 1:
            raise ValueError('splice DB {} length isn\'t 1'.format(path))
        db = db[0].split()
        db = [db[i:i+4] for i in range(0, len(db), 4)]
        self.start, self.end, self.chainbreak, self.name = db.pop()
        self.start = int(self.start)
        self.end = int(self.end)
        self.chainbreak = int(self.chainbreak)
        self.db = pd.DataFrame(db, columns=['phi', 'psi', 'omega', 'res'])

    def update_db(self, design):
        """
        Changes the sequence of the db
        :param design: a splice out design which was generated together with
        this database.
        """
        design = parse_pdb_file(design)
        chain = list(design.get_chains())[0].id
        design_residues = design[0][chain]
        seq = [design_residues[i].resname
               for i in range(self.start, self.start + len(self.db))]
        assert len(seq) == len(self.db)
        self.db['res'] = seq

    def write_db(self, path=''):
        """
        :param path: path to write the new database to
        :return:
        """
        if not path:
            path = self.original_path + '_new.db'
        db = ' '.join([' '.join(r) for i, r in self.db.iterrows()])
        db += f' {self.start} {self.end} {self.chainbreak} {self.name}\n'
        open(path, 'w').write(db)

    def get_sequence(self):
        """Returns the sequence of the database (the forth element of each
        residue"""
        seq = ''.join([seq1(r) for r in self.db['res']])
        return seq

    def __len__(self):
        return self.db.shape[0]

def get_chimera_limits(blades: List[SpliceDB]) -> dict:
    """
    :param blades: a list of SpliceDB objects of all spliced segments in the
    protein
    :return: data frame of ['name', 'start', 'end'] for each inserted blade
    with the actual start and end [start, end] of the segments
    """
    limits = pd.DataFrame(columns=['name', 'start', 'end'])
    limits.loc[len(limits)] = [blades[0].name, blades[0].start,
                               blades[0].start + len(blades[0]) - 1]
    prev_end = limits.iloc[-1]['end']
    for i in range(1, len(blades)):
        fixed_template = blades[i].start - blades[i-1].end - 1
        new_start = prev_end + 1 + fixed_template
        new_end = new_start + len(blades[i]) - 1
        limits.loc[len(limits)] = [blades[i].name, new_start, new_end]
        prev_end = new_end
    return limits

if __name__ == '__main__':
    # Example of usage to calculate limits in a chimeric protein
    s1 = SpliceDB('name_s1.db')  # name_s1.db is tosrion db from AbDesign
    s2 = SpliceDB('name_s2.db')
    s3 = SpliceDB('name_s3.db')
    s4 = SpliceDB('name_s4.db')
    dbs = [s1, s2, s3, s4]
    get_chimera_limits(dbs)
