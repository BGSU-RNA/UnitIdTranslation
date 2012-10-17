#!/usr/bin/env python

from __future__ import with_statement

import os
import csv
import sys
from itertools import ifilter

from CorePyWrap import ParseCif

from Bio.PDB.PDBParser import PDBParser

from pprint import pprint


NEW_FIELDS = ['pdb', 'model', 'chain', 'unit', 'number', 'atom_name', 'alt_id',
              'insertion', 'operator']

OLD_FIELDS = ['pdb', 'type', 'model', 'chain', 'number', 'unit', 'insertion']


def table(block, name):
    table = block.GetTable(name)
    entries = {}
    for row in rows(table):
        entries[row['id']] = row
    return entries


def row(table, index):
    columns = table.GetColumnNames()
    return dict(zip(columns, table.GetRow(index)))


def rows(table):
    columns = table.GetColumnNames()
    for index in range(table.GetNumRows()):
        data = dict(zip(columns, table.GetRow(index)))
        yield data


def build_translation_table(filename):
    """Given a file object of a CIF file this will produce a translation table
    usable by translate
    """
    translation_table = {}
    parser = ParseCif(filename)

    for pdb_id in parser.GetBlockNames():
        translation_table[pdb_id] = {}
        pdb = parser.GetBlock(pdb_id)
        gen_assemblies = pdb.GetTable('pdbx_struct_assembly_gen')
        operator_table = table(pdb, 'pdbx_struct_oper_list')

        for gen_row in rows(gen_assemblies):
            assembly_id = gen_row['assembly_id']
            translation_table[pdb_id][assembly_id] = {}
            models = gen_row['oper_expression'].split(',')

            for model in models:
                name = operator_table[model]['name']
                translation_table[pdb_id][assembly_id][model] = name
    return translation_table


def load_translation_table(raw):
    names = ['pdb', 'ba', 'model', 'sym_name', 'sym_frac']
    reader =  DictReader(raw, names, delimiter='\t')
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
                name = residue.resname.strip()
                res_id = residue.get_id()
                if residue.has_id("C1'"):
                    data.append({
                        'pdb': pdb_id,
                        'type': pdb_type,
                        'model': model_id,
                        'chain': chain_id,
                        'number': str(res_id[1]),
                        'unit': name,
                        'insertion': res_id[2].rstrip()
                    })

    return data


def as_new_id(new_id):
    fields = NEW_FIELDS
    if not new_id.get('operator', None) or new_id['operator'] == '1_555':
        fields = fields[:-1]
        if not new_id.get('insertion', None):
            fields = fields[:-3]

    return '_'.join([new_id.get(part, '') for part in fields])


def as_old_id(old_id):
    return '_'.join([old_id.get(part, '') for part in OLD_FIELDS])


def get_id_correspondences(pdb_file, cif_file):
    translation_table = build_translation_table(cif_file)
    with open(pdb_file, 'r') as raw:
        old_ids = old_residue_ids(cif_file)

    new_ids = tranlsation(old_ids, translation_table)
    correspondences = []
    for index, old_id in enumerate(old_ids):
        correspondences.append((as_old_id(old_id), as_new_id(new_ids[index])))
    return correpondecies


def main(cif_file):
    translated = get_id_correspondences(cif_file)
    pprint(translated)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

# def sort(a):
#     parts = a.split('_')
#     if '__' in a:
#         return '_'.join([parts[0], parts[1], parts[2], parts[3], parts[4], ''])
#     return '_'.join([parts[0], parts[2], parts[3], parts[4], parts[5], parts[6]])

# def make_matches(target, possible):
#     found = list(possible)
#     if len(found) == 0:
#         return [(target, None)]
#     if len(found) > 1:
#         raise MappingError("%s maps to more than id" % target)
#     return [(target, match) for match in found]

# def get_assemblies(block):
#     gen_assemblies = block.GetTable('pdbx_struct_assembly_gen')
#     operator_table = table(block, 'pdbx_struct_oper_list')
#     assemblies = {}
#     for gen_row in rows(gen_assemblies):
#         ids = gen_row['asym_id_list'].split(',')
#         oper_ids = gen_row['oper_expression'].split(',')
#         oper_names = [operator_table[oper_id]['name'] for oper_id in oper_ids]
#         for asym_id in ids:
#             if asym_id in assemblies:
#                 raise ValueError("Duplicate assembly_gen")
#             assemblies[asym_id] = {
#                 'ba_id': gen_row['assembly_id'],
#                 'assembly_id': gen_row['assembly_id'],
#                 'operator_names': oper_names
#             }
#     return assemblies

# def new_residue_ids(cif_file):
#     cif = ParseCif(cif_file)
#     blocks = cif.GetBlockNames()
#     parts = [
#       ('pdbx_PDB_model_num', 'model'),
#       ('auth_asym_id', 'chain'),
#       ('label_comp_id', 'unit'),
#       ('auth_seq_id', 'number'),
#       # TODO: Note we are assuming residues are are not in alternate locations
#       # 'label_alt_id', Should be ok I think
#       ('pdbx_PDB_ins_code', 'insertion'),
#     ]
#
#     ids = []
#     seen = set()
#     for pdb_name in blocks:
#         block = cif.GetBlock(pdb_name)
#         assemblies = get_assemblies(block)
#         atoms = block.GetTable('atom_site')
#         for atom in rows(atoms):
#             data = {}
#             for part in parts:
#                 value = atom[part[0]]
#                 if value == '?' or value == '.':
#                     value = ''
#                 data[part[1]] = value
#
#             assembly = assemblies[atom['label_asym_id']]
#             for name in assembly['operator_names']:
#                 add = dict(data)
#                 add['pdb'] = pdb_name
#                 add['operator'] = name
#                 add['ba_id'] = str(assembly['ba_id'])
#                 ids.append(add)
#
#     # Uniqueify atom ids to get residue ids
#     seen = {}
#     order = sorted(ids[0].keys())
#     for atom in ids:
#         data = '_'.join([atom[key] for key in order])
#         if data not in seen:
#             seen[data] = atom
#
#     return seen.values()

# The translation table only builds stuff for biological assemblies, we can
# assume that all AU's are 1_555.
# def build_translation_table(new_ids):
#     translation = defaultdict(lambda : defaultdict(dict))
#     for nid in new_ids:
#         pdb = nid['pdb']
#         ba = nid['ba_id']
#         model = nid[]
#         name = nid['operator']
#         translation[pdb][ba][model] = name
#     return dict(translation)

# In the case of asymmetric units comparing is easy:
#   An old and new ID are the same if the pdb, model, chain, number, insertion
#   code are
# In other cases two ids are the same if the above is true but:
#  The biological assembly number must be the same as the
#  pdbx_struct_assembly_gen.assembly_id - 1, which is ba_id in the dict
# def translate(old, new):
#     mappings = []
#     def matches(old, key):
#         return lambda new: old[key] == new[key]
#
#     def ba_match(old):
#         ba = int(old['type'][-1])
#         return lambda new: ba == new['ba_id']
#
#     for old_id in old:
#         possible = ifilter(matches(old_id, 'pdb'), new)
#         possible = ifilter(matches(old_id, 'model'), possible)
#         possible = ifilter(matches(old_id, 'chain'), possible)
#         possible = ifilter(matches(old_id, 'number'), possible)
#         possible = ifilter(matches(old_id, 'insertion'), possible)
#         possible = ifilter(matches(old_id, 'unit'), possible)
#
#         if old_id['type'] == 'AU':
#             mappings.extend(make_matches(old_id, possible))
#
#         else:
#           possible = ifilter(ba_match(old_id), possible)
#           mappings.extend(make_matches(old_id, possible))
#
#     return mappings
