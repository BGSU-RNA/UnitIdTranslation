from __future__ import with_statement

import os
import csv
import sys
from collections import defaultdict

from pdbx.reader.PdbxParser import PdbxReader

from Bio.PDB.PDBParser import PDBParser

from pprint import pprint

# NEW_FIELDS = ['pdb', 'model', 'chain', 'unit', 'number', 'atom_name', 'alt_id',
#               'insertion', 'operator']

NEW_FIELDS = ['pdb', 'model', 'chain', 'unit', 'number', 'alt_id', 'insertion',
              'operator']

OLD_FIELDS = ['pdb', 'type', 'model', 'chain', 'number', 'unit', 'insertion']


def table(block, name):
    entries = {}
    for row in rows(block, name):
        entries[row['id']] = row
    return entries


def rows(block, name):
    entry = block.getObj(name)
    columns = entry.getItemNameList()
    columns = [column.replace('_' + name + '.', '') for column in columns]
    for index in range(entry.getRowCount()):
        yield dict(zip(columns, entry.getRow(index)))


def build_translation_table(filename):
    """Given a file object of a CIF file this will produce a translation table
    usable by translate
    """
    translation_table = {}
    data = []
    with open(filename, 'r') as raw:
        parser = PdbxReader(raw)
        parser = parser.read(data)
    pdb = data[0]
    pdb_id = pdb.getName()

    translation_table[pdb_id] = {}
    # pdb = parser.GetBlock(pdb_id)
    operator_table = table(pdb, 'pdbx_struct_oper_list')

    for gen_row in rows(pdb, 'pdbx_struct_assembly_gen'):
        assembly_id = gen_row['assembly_id']
        translation_table[pdb_id][assembly_id] = {}
        models = gen_row['oper_expression'].split(',')

        for model in models:
            name = operator_table[model]['name']
            translation_table[pdb_id][assembly_id][model] = name
    return translation_table


def load_translation_table(raw):
    names = ['pdb', 'ba', 'model', 'sym_name', 'sym_frac']
    reader =  csv.DictReader(raw, names, delimiter='\t')
    data = defaultdict(lambda : defaultdict(dict))
    for entry in reader:
        data[entry['pdb']][entry['ba']][entry['model']] = entry['sym_name']
    return data


def translate(old_ids, table):
    new = []
    for oid in old_ids:
        nid = dict(oid)
        ntype = nid['type']
        del nid['type']
        if ntype != 'AU':
            ba = ntype[-1]
            name = table[nid['pdb']][ba][nid['model']]
            nid['operator'] = name
        new.append(nid)
    return new


def old_residue_ids(raw, filename):
    parser = PDBParser()
    path, ext = os.path.splitext(filename)
    pdb_id = os.path.basename(path)
    structure = parser.get_structure(pdb_id, raw)
    data = []

    pdb_type = 'AU'
    if ext != '.pdb':
        pdb_type = 'BA' + filename[-1]

    for model in structure:
        # BioPython seems to start number models at 0, but it should start
        # at 1.
        model_id = str(model.get_id() + 1)
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_id = residue.get_id()
                data.append({
                    'pdb': pdb_id,
                    'type': pdb_type,
                    'model': model_id,
                    'chain': chain_id,
                    'number': str(res_id[1]),
                    'unit': residue.resname.strip(),
                    'insertion': res_id[2].rstrip()
                })

    return data


def as_new_id(new_id):
    fields = NEW_FIELDS
    if not new_id.get('operator', None) or new_id['operator'] == '1_555':
        fields = fields[:-1]
        if not new_id.get('insertion', None):
            fields = fields[:-2]

    return '_'.join([new_id.get(part, '') for part in fields])


def as_old_id(old_id):
    return '_'.join([old_id.get(part, '') for part in OLD_FIELDS])


def get_id_correspondences(pdb_file, cif_file):
    translation_table = build_translation_table(cif_file)
    with open(pdb_file, 'r') as raw:
        old_ids = old_residue_ids(raw, pdb_file)

    new_ids = translate(old_ids, translation_table)
    correspondences = []
    for index, old_id in enumerate(old_ids):
        correspondences.append((as_old_id(old_id), as_new_id(new_ids[index])))
    return correspondences


def main(pdb_file, cif_file):
    translated = get_id_correspondences(pdb_file, cif_file)
    pprint(translated)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
