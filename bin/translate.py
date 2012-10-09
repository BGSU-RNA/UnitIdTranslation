#!/usr/bin/env python

import sys

from CorePyWrap import ParseCif

def parse_cif(cif_file):
    cif = ParseCif(cif_file)
    blocks = cif.GetBlockNames()
    parts = [
      ('atom_site', 'pdbx_PDB_model_num'),
      ('atom_site', 'auth_asym_id'),
      ('atom_site', 'label_comp_id'),
      ('atom_site', 'auth_seq_id')
      # ('atom_site', 'label_atom_id'),
      # ('atom_site', 'label_alt_id'),
      # ('atom_site', 'pdbx_PDB_ins_code', '')
    ]
    ids = []
    for name in blocks:
        block = cif.GetBlock(name)
        operators = block.GetTable('pdbx_struct_oper_list').GetColumn('name')
        operator = operators[0]
        for part in parts:
            table = block.GetTable(part[0])
            data = table.GetColumn(part[1])
            for index, value in enumerate(data):
                if len(ids) <= index:
                    ids.append([name])
                if value == '?':
                    if len(part) > 2:
                        value = part[2]
                    else:
                        raise ValueError("Can't place default for: %s, %s" % part)
                if value:
                    ids[index].append(value)
        for unit_id in ids:
            unit_id.append(operator)
    ids = ['_'.join(unit_id) for unit_id in ids]
    return ids



def main(cif_file):
    ids = parse_cif(cif_file)
    print(ids)
        


if __name__ == '__main__':
    main(sys.argv[1])
