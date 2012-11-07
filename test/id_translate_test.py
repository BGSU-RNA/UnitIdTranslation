import sys
import unittest
import warnings

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
               '1D4R_2_A_C_2__a_6_555']
        self.assertEqual(val, ans)

    def test_old_id_string(self):
        val = [as_old_id(data) for data in OLD]
        ans = ['2AW7_AU_1_A_1_A_', '2AW7_AU_1_A_2_C_', '1D4R_BA1_1_A_3_A_',
               '1D4R_BA2_2_A_2_C_a']
        self.assertEqual(val, ans)


class Test3D2S(unittest.TestCase):
    def test_build_translation_table(self):
        val = build_translation_table('test/files/3D2S.cif')
        ans = {
            '3D2S': {
                '1': { '1': '2_645', '2': '1_555', '3': '1_554', },
                '2': { '2': '1_555' },
                '3': { '2': '1_555', '3': '1_554', '4': '2_545', },
                '4': { '2': '1_555' },
            }
        }
        self.assertEqual(val, ans)


class Test3DSU(unittest.TestCase):
    def test_build_translation_table(self):
        val = build_translation_table('test/files/3DSU.cif')
        ans = {
            '3DSU': {
                '1': { '1': '1_555' }
            }
        }
        self.assertEqual(val, ans)


class Test3QSU(unittest.TestCase):
    def test_build_translation_table(self):
        val = build_translation_table('test/files/3QSU.cif')
        ans = {
            '3QSU': {
                '1': { '1': '1_556', '2': '1_555' },
                '2': { '1': '1_556', '2': '1_555' },
                '3': { '2': '1_555', '3': '2_555', '4': '3_555' }
            }
        }
        self.assertEqual(val, ans)


class Test1AQ4(unittest.TestCase):
    def test_fails_translation_table(self):
        func = lambda: build_translation_table('test/files/1AQ4.cif')
        self.assertRaises(LooksLikeAVirusStructureError, func)


class Test2Z4L(unittest.TestCase):
    # In this case the method uses nested default dicts so that all access to
    # the table produce '1_555'. I don't test that it uses a defaultdict
    # directly since this is likely to change in the future.
    def test_builds_table_with_missing_assembly(self):
        val = build_translation_table('test/files/2Z4L.cif')
        self.assertEquals(val['1']['1'], '1_555')
        self.assertEquals(val['3']['1'], '1_555')
        self.assertEquals(val['1']['10'], '1_555')

    # Smoke test to see if this will crash when generating IDs. It was having
    # trouble before so I'm adding a silly test for it.
    def test_does_not_crash_making_pdb_ids(self):
        # BioPython spits out annoying warnings, so we silence them.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            val = get_id_correspondences('test/files/2Z4L.pdb',
                                         'test/files/2Z4L.cif')
            self.assertEquals(val, val)

    def test_does_not_crash_making_pdb_ids(self):
        # BioPython spits out annoying warnings, so we silence them.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            val = get_id_correspondences('test/files/2Z4L.pdb1',
                                         'test/files/2Z4L.cif')
            self.assertEquals(val, val)
