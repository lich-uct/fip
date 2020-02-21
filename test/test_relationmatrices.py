#!/usr/bin/env python3
import unittest

import pandas as pd
import numpy as np

from fiprofiling import relationmatrices as rm

class TestRM(unittest.TestCase):
    def test_instantiation(self):
        relmatrix = rm.RM()
        self.assertIsInstance(relmatrix, rm.RM)

    def test_from_fingerprints(self):
        with self.assertRaises(NotImplementedError):
            rm.RM.from_fingerprints(None)

    def test_from_dataframe(self):
        relmatrix = rm.RM(pd.DataFrame([[1]*4]*4))
        self.assertIsInstance(relmatrix, rm.RM)

    def test_is_square(self):
        relmatrix = rm.RM(pd.DataFrame([[1]*4]*4))
        self.assertTrue(relmatrix.is_square())
        relmatrix = rm.RM(pd.DataFrame([[1]*4]*3))
        self.assertFalse(relmatrix.is_square())

    def test_is_symmetric(self):
        relmatrix = rm.RM(pd.DataFrame(
            [[1, 2, 3],
             [2, 4, 5],
             [3, 5, 8]]
        ))
        self.assertTrue(relmatrix.is_symmetric())
        relmatrix = rm.RM(pd.DataFrame(
            [[1, 2, 3],
             [9, 4, 5],
             [3, 5, 8]]
        ))
        self.assertFalse(relmatrix.is_symmetric())

    def test_basic_mean(self):
        rmlist = [[1, 2, 3],
                  [9, 4, 5],
                  [3, 5, 8]]
        relmatrix = rm.RM(pd.DataFrame(rmlist))
        mean = float(sum([sum(x) for x in rmlist])) / sum([len(x) for x in rmlist])
        self.assertEqual(relmatrix.mean(offdiag_only=False), mean)

    def test_offdiag_only_mean(self):
        rmlist = [[1, 2, 3],
                  [9, 4, 5],
                  [3, 5, 8]]
        relmatrix = rm.RM(pd.DataFrame(rmlist))
        offdiag = [2, 3, 9, 5, 3, 5]
        mean = float(sum(offdiag)) / len(offdiag)
        self.assertEqual(relmatrix.mean(offdiag_only=True), mean)

    def test_difference(self):
        relmatrix1 = rm.RM(pd.DataFrame([[1, 2, 3],
                                         [2, 4, 5],
                                         [3, 5, 8]]),
                                   num_datapoints=10)
        relmatrix2 = rm.RM(pd.DataFrame([[1, 2, 0],
                                         [9, 4, 5],
                                         [3, 7, 8]]),
                                   num_datapoints=9)
        expected_result_rm = rm.RM(pd.DataFrame([[0, 0, 3],
                                                 [-7, 0, 0],
                                                 [0, -2, 0]]),
                                           num_datapoints=9)
        diffmatrix = relmatrix1.difference(relmatrix2)
        self.assertTrue((diffmatrix.df == expected_result_rm.df).all().all())
        self.assertEqual(diffmatrix.num_datapoints, expected_result_rm.num_datapoints)

    def test_distance(self):
        relmatrix1 = rm.RM(pd.DataFrame([[1, 2, 3],
                                         [2, 4, 5],
                                         [3, 5, 8]]),
                                   num_datapoints=10)
        relmatrix2 = rm.RM(pd.DataFrame([[1, 2, 0],
                                         [9, 4, 5],
                                         [3, 7, 8]]),
                                   num_datapoints=9)
        expected_result_rm = rm.RM(pd.DataFrame([[0, 0, 3],
                                                 [-7, 0, 0],
                                                 [0, -2, 0]]),
                                           num_datapoints=9)
        distance = relmatrix1.distance(relmatrix2, type='euclidean')
        self.assertEqual(distance, relmatrix2.distance(relmatrix1, type='euclidean'),
                         "|AB| == |BA|")
        self.assertEqual(distance, np.sqrt(sum([x**2 for x in (3, -7, -2)])),
                         "correct euclidean dist")

    def test_triangular(self):
        fullmatrix = rm.RM(pd.DataFrame([[1, 2, 3],
                                         [2, 4, 5],
                                         [3, 5, 8]]),
                                   num_datapoints=10)
        upper_triangular_df = pd.DataFrame([[1, 2, 3],
                                            [np.nan, 4, 5],
                                            [np.nan, np.nan, 8]])
        upper_triangular_rm = fullmatrix.triangular(strict=False)
        self.assertTrue((upper_triangular_rm.df.fillna(0)
                         == upper_triangular_df.fillna(0)).all().all())
        strict_upper_triangular_df = pd.DataFrame([[np.nan, 2, 3],
                                                   [np.nan, np.nan, 5],
                                                   [np.nan, np.nan, np.nan]])
        upper_triangular_rm = fullmatrix.triangular(strict=True)
        self.assertTrue((upper_triangular_rm.df.fillna(0)
                         == strict_upper_triangular_df.fillna(0)).all().all())

class TestCORM(unittest.TestCase):
    def test_instantiation(self):
        relmatrix = rm.CORM()
        self.assertIsInstance(relmatrix, rm.CORM)
        self.assertIsInstance(relmatrix, rm.RM)

    def test_from_dataframe(self):
        relmatrix = rm.CORM(pd.DataFrame([[1]*4]*4))
        self.assertIsInstance(relmatrix, rm.CORM)
        self.assertTrue((relmatrix.df == pd.DataFrame([[1]*4]*4)).all().all())

    def test_from_fingerprints(self):
        fps = ['00000000',
               '00001111',
               '00110011',
               '01010101',
               '11111111']
        relmatrix = rm.CORM.from_fingerprints(fps, fpformat='bintext')
        coorm = rm.CORM(
            pd.DataFrame([[1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 2, 1, 2, 1, 2, 1, 2],
                          [1, 1, 2, 2, 1, 1, 2, 2],
                          [1, 2, 2, 3, 1, 2, 2, 3],
                          [1, 1, 1, 1, 2, 2, 2, 2],
                          [1, 2, 1, 2, 2, 3, 2, 3],
                          [1, 1, 2, 2, 2, 2, 3, 3],
                          [1, 2, 2, 3, 2, 3, 3, 4]],
                          index=[str(x) for x in range(8)],
                          columns=[str(x) for x in range(8)]),
            num_datapoints=len(fps))
        self.assertTrue(relmatrix.is_equal(coorm))

    def test_add_fingerprint(self):
        fps = ['00000000',
               '00001111',
               '00110011',
               '01010101',
               '11111111']
        relmatrix = rm.CORM.from_fingerprints([fps[0]], fpformat='bintext')
        for fp in fps[1:]:
            relmatrix.add_fingerprint(fp, fpformat='bintext')
        coorm = rm.CORM(
            pd.DataFrame([[1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 2, 1, 2, 1, 2, 1, 2],
                          [1, 1, 2, 2, 1, 1, 2, 2],
                          [1, 2, 2, 3, 1, 2, 2, 3],
                          [1, 1, 1, 1, 2, 2, 2, 2],
                          [1, 2, 1, 2, 2, 3, 2, 3],
                          [1, 1, 2, 2, 2, 2, 3, 3],
                          [1, 2, 2, 3, 2, 3, 3, 4]],
                          index=[str(x) for x in range(8)],
                          columns=[str(x) for x in range(8)]),
            num_datapoints=len(fps))
        self.assertTrue(relmatrix.is_equal(coorm))

    def test_addition(self):
        coorm1_len = 10
        coorm1 = rm.CORM(
            pd.DataFrame([[1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 2, 1, 2, 1, 2, 1, 2],
                          [1, 1, 2, 2, 1, 1, 2, 2],
                          [1, 2, 2, 3, 1, 2, 2, 3],
                          [1, 1, 1, 1, 2, 2, 2, 2],
                          [1, 2, 1, 2, 2, 3, 2, 3],
                          [1, 1, 2, 2, 2, 2, 3, 3],
                          [1, 2, 2, 3, 2, 3, 3, 4]]),
            num_datapoints=coorm1_len)
        coorm2_len = 11
        coorm2 = rm.CORM(
            pd.DataFrame([[1, 1, 1, 1, 1, 1, 1, 1],
                          [1, 2, 1, 2, 1, 2, 1, 2],
                          [1, 1, 4, 2, 1, 1, 2, 2],
                          [1, 2, 2, 3, 1, 2, 2, 3],
                          [1, 1, 1, 1, 3, 2, 2, 2],
                          [1, 2, 1, 2, 2, 3, 2, 3],
                          [1, 1, 3, 2, 2, 2, 8, 3],
                          [1, 2, 2, 3, 2, 3, 3, 4]]),
            num_datapoints=coorm2_len)
        sumcoorm = coorm1 + coorm2
        self.assertEqual(sumcoorm.num_datapoints,
                         coorm1.num_datapoints + coorm2.num_datapoints)
        self.assertTrue((sumcoorm.df == (coorm1.df + coorm2.df)).all().all())


class TestCOPRM(unittest.TestCase):
    def test_instantiation(self):
        relmatrix = rm.COPRM()
        self.assertIsInstance(relmatrix, rm.COPRM)
        self.assertIsInstance(relmatrix, rm.RM)

    def test_from_CORM(self):
        fps = ['00000000',
               '00001111',
               '00110011',
               '01010101',
               '11111111']
        coorm = rm.CORM.from_fingerprints(fps, fpformat='bintext')
        probrm = rm.COPRM.from_CORM(coorm)
        probrm2= rm.COPRM(
            pd.DataFrame([[0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
                          [0.2, 0.4, 0.2, 0.4, 0.2, 0.4, 0.2, 0.4],
                          [0.2, 0.2, 0.4, 0.4, 0.2, 0.2, 0.4, 0.4],
                          [0.2, 0.4, 0.4, 0.6, 0.2, 0.4, 0.4, 0.6],
                          [0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4],
                          [0.2, 0.4, 0.2, 0.4, 0.4, 0.6, 0.4, 0.6],
                          [0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.6, 0.6],
                          [0.2, 0.4, 0.4, 0.6, 0.4, 0.6, 0.6, 0.8]],
                          index=[str(x) for x in range(8)],
                          columns=[str(x) for x in range(8)]),
            num_datapoints=len(fps))
        self.assertTrue(probrm.is_equal(probrm2))


class TestPMIRM(unittest.TestCase):
    def test_instantiation(self):
        pmirm = rm.PMIRM()
        self.assertIsInstance(pmirm, rm.PMIRM)
        self.assertIsInstance(pmirm, rm.RM)

    def test_from_CORM(self):
        fps = ['00000000',
               '00001111',
               '00110011',
               '01010101',
               '11111111']
        coorm = rm.CORM.from_fingerprints(fps, fpformat='bintext')
        pmirm = rm.PMIRM.from_CORM(coorm)
        probrm = rm.COPRM.from_CORM(coorm)
        for a in range(pmirm.df.shape[0]):
            for b in range(pmirm.df.shape[1]):
                if a == b:
                    independent_p_a_b = probrm.df.iloc[a, b]
                else:
                    independent_p_a_b = probrm.df.iloc[a, a]*probrm.df.iloc[b, b]
                p_a_b = probrm.df.iloc[a, b]
                pmirm_a_b = np.log2(p_a_b/independent_p_a_b)
                self.assertEqual(pmirm.df.iloc[a, b], pmirm_a_b)


class TestZPMIRM(unittest.TestCase):
    def test_instantiation(self):
        zpmirm = rm.ZPMIRM()
        self.assertIsInstance(zpmirm, rm.ZPMIRM)
        self.assertIsInstance(zpmirm, rm.PMIRM)
        self.assertIsInstance(zpmirm, rm.RM)

    def test_from_pmirm(self):
        fps = ['00000000',
               '00001111',
               '00110011',
               '01010101',
               '11111111']
        coorm = rm.CORM.from_fingerprints(fps, fpformat='bintext')
        pmirm = rm.PMIRM.from_CORM(coorm)
        zpmirm = rm.ZPMIRM.from_PMIRM(pmirm)
        for a in range(zpmirm.df.shape[0]):
            for b in range(zpmirm.df.shape[1]):
                if a == b:
                    self.assertEqual(zpmirm.df.iloc[a, b], 0)


if __name__ == '__main__':
    unittest.main()
