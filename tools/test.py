#!/usr/bin/env python

import unittest
import tempfile, shutil, os, filecmp

from swarm_util import *
import fix_datafile
import slice_data
import merge_data

PATH = os.path.dirname(os.path.abspath(__file__))

CORRECT_FILE = os.path.join(PATH, './testdata/correct.ptcl')
INCORRECT_LONG_FILE = os.path.join(PATH, './testdata/incorrect_long.ptcl')
INCORRECT_SHORT_FILE = os.path.join(PATH, './testdata/incorrect_short.ptcl')

class Test(unittest.TestCase):

    def create_temporay_file(self):
        tmp_file = tempfile.mkstemp()[1]
        return tmp_file

    def create_temporary_copy(self, path):
        tmp_file = self.create_temporay_file()
        shutil.copy(path, tmp_file)
        os.chmod(tmp_file, 0o600)
        return tmp_file

    def test_create_temporary_copy(self):
        target_file = self.create_temporary_copy(CORRECT_FILE)
        cmp = filecmp.cmp(CORRECT_FILE, target_file)
        self.assertTrue(cmp)
        os.remove(target_file)

    def test_check_data_correct(self):
        sdm = SwarmDataManager(CORRECT_FILE)
        sdm.check_metadata()
        sdm.close()

    def test_check_data_incorrect_long(self):
        sdm = SwarmDataManager(INCORRECT_LONG_FILE)
        self.assertRaises(Exception, sdm.check_metadata)
        sdm.close()

    def test_check_data_incorrect_short(self):
        sdm = SwarmDataManager(INCORRECT_SHORT_FILE)
        self.assertRaises(Exception, sdm.check_metadata)
        sdm.close()

    def test_fix_data_correst(self):
        target_file = self.create_temporary_copy(CORRECT_FILE)
        fix_datafile.fix_data(target_file)
        self.assertTrue(filecmp.cmp(CORRECT_FILE, target_file))
        os.remove(target_file)

    def test_fix_data_incorrest_long(self):
        target_file = self.create_temporary_copy(INCORRECT_LONG_FILE)
        fix_datafile.fix_data(target_file)
        self.assertTrue(filecmp.cmp(CORRECT_FILE, target_file))
        os.remove(target_file)

    def test_fix_data_incorrest_short(self):
        target_file = self.create_temporary_copy(INCORRECT_SHORT_FILE)
        fix_datafile.fix_data(target_file)
        sdm = SwarmDataManager(target_file)
        sdm.check_metadata()
        sdm.close()
        os.remove(target_file)

    def test_slice_data(self):
        target_file = self.create_temporay_file()
        slice_data.slice_data(CORRECT_FILE, target_file, 100, 50)
        sdm = SwarmDataManager(target_file)
        sdm.check_metadata()
        sdm.close()
        os.remove(target_file)

    def test_slice_data_invalid_step(self):
        target_file = self.create_temporay_file()
        self.assertRaises(Exception, lambda: slice_data.slice_data(CORRECT_FILE, target_file, 100, 10000))
        os.remove(target_file)

    def test_slice_data_invalid_start(self):
        target_file = self.create_temporay_file()
        self.assertRaises(Exception, lambda: slice_data.slice_data(CORRECT_FILE, target_file, 10000, 10))
        os.remove(target_file)

    def test_split_merge_data(self):
        target_file1 = self.create_temporay_file()
        target_file2 = self.create_temporay_file()
        target_file3 = self.create_temporay_file()
        target_file4 = self.create_temporay_file()
        slice_data.slice_data(CORRECT_FILE, target_file1, 0, 500)
        slice_data.slice_data(CORRECT_FILE, target_file2, 500, 1000)
        merge_data.merge_data([target_file1, target_file2], target_file3)
        slice_data.slice_data(CORRECT_FILE, target_file4, 0, 1500)
        result = filecmp.cmp(target_file3, target_file4)
        self.assertTrue(result)

    def test_split_merge_data_2(self):
        sdm = SwarmDataManager(CORRECT_FILE)
        N = sdm.N
        sdm.close()
        target_file1 = self.create_temporay_file()
        target_file2 = self.create_temporay_file()
        target_file3 = self.create_temporay_file()
        slice_data.slice_data(CORRECT_FILE, target_file1, 0, N//2)
        slice_data.slice_data(CORRECT_FILE, target_file2, N//2, N - N//2)
        merge_data.merge_data([target_file1, target_file2], target_file3)
        result = filecmp.cmp(CORRECT_FILE, target_file3)
        self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()
