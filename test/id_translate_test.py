import sys
import unittest

from id_translate import *


OLD = [
    {
        'pdb': '2AW7',
        'type': 'AU',
        'model': '1',
        'chain': 'A',
        'number': '1',
        'unit': 'A',
        'insertion': ''
    },
    {
        'pdb': '2AW7',
        'type': 'AU',
        'model': '1',
        'chain': 'A',
        'number': '2',
        'unit': 'C',
        'insertion': ''
    },
    {
        'pdb': '1D4R',
        'type': 'BA1',
        'model': '1',
        'chain': 'A',
        'number': '3',
        'unit': 'A',
        'insertion': ''
    },
    {
        'pdb': '1D4R',
        'type': 'BA2',
        'model': '2',
        'chain': 'A',
        'number': '2',
        'unit': 'C',
        'insertion': 'a'
    },
]

NEW = [
    {
        'pdb': '2AW7',
        'model': '1',
        'chain': 'A',
        'number': '1',
        'unit': 'A',
        'insertion': ''
    },
    {
        'pdb': '2AW7',
        'model': '1',
        'chain': 'A',
        'number': '2',
        'unit': 'C',
        'insertion': ''
    },
    {
        'pdb': '1D4R',
        'model': '1',
        'chain': 'A',
        'number': '3',
        'unit': 'A',
        'insertion': '',
        'operator': '1_555'
    },
    {
        'pdb': '1D4R',
        'model': '2',
        'chain': 'A',
        'number': '2',
        'unit': 'C',
        'insertion': 'a',
        'operator': '6_555'
    },
]


class TestTranslate(unittest.TestCase):
    def setUp(self):
        pass

    def test_build_translation_table(self):
        ans = {
            '1D4R': {
                '1': { '1': '1_555', },
                '2': { '1': '1_555', '2': '6_555' }
            }
        }
        val = build_translation_table('test/files/1D4R.cif')
        self.assertEqual(val, ans)

    def test_translate_au(self):
        ids = OLD[0:2]
        ans = NEW[0:2]
        val = translate(ids, None)
        self.assertEqual(val, ans)

    def test_translate_ba(self):
        table = {
            '1D4R': {
                '1': { '1': '1_555', },
                '2': { '1': '1_555', '2': '6_555' }
            }
        }
        ids = OLD[3:5]
        ans = NEW[3:5]
        val = translate(ids, table)
        self.assertEqual(val, ans)

    def test_translate(self):
        table = {
            '1D4R': {
                '1': { '1': '1_555', },
                '2': { '1': '1_555', '2': '6_555' }
            }
        }
        val = translate(OLD, table)
        self.assertEqual(val, NEW)

    def test_new_id_string(self):
        val = [as_new_id(data) for data in NEW]
        ans = ['2AW7_1_A_A_1', '2AW7_1_A_C_2', '1D4R_1_A_A_3', 
               '1D4R_2_A_C_2___a_6_555']
        self.assertEqual(val, ans)

    def test_old_id_string(self):
        val = [as_old_id(data) for data in OLD]
        ans = ['2AW7_AU_1_A_1_A_', '2AW7_AU_1_A_2_C_', '1D4R_BA1_1_A_3_A_',
               '1D4R_BA2_2_A_2_C_a']
        self.assertEqual(val, ans)
