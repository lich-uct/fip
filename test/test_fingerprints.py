#!/usr/bin/env python3
import unittest

from fiprofiling import chem
from fiprofiling import fingerprints as fps

class TestFingerprints(unittest.TestCase):
    def test_boollist2fp(self):
        fp = fps.boollist2fp([False, False, True, True])
        self.assertEqual(len(fp), 4)
        self.assertEqual(set(fps.fp2onbits(fp)), set([2, 3]))

    def test_bintext2fp(self):
        fp = fps.bintext2fp("0101")
        self.assertEqual(len(fp), 4)
        self.assertEqual(set(fps.fp2onbits(fp)), set([1, 3]))

    def test_hextext2fp(self):
        hextext = '000000000800443a0445b4c83be9bacb9b632fef1d'
        fp = fps.hextext2fp(hextext)
        self.assertEqual(len(fp), len(hextext)*4)

    def test_fp2hextext(self):
        hextext = '000000000800443a0445b4c83be9bacb9b632fef1d'
        fp = fps.hextext2fp(hextext)
        self.assertEqual(fps.fp2hextext(fp), hextext)

    def test_bintext2fp2bintext(self):
        bintext = "001101010101010001010"
        fp = fps.bintext2fp(bintext)
        self.assertEqual(fps.fp2bintext(fp), bintext)

    def test_fp2nparray(self):
        bintext = "001101010101010001010"
        fp = fps.bintext2fp(bintext)
        array = fps.fp2nparray(fp)
        for char, arraynum in zip(bintext, array):
            self.assertEqual(arraynum, int(char))

    def test_default_tanimoto(self):
        bintexts_results = [("0011", "0101", 1.0/3),
                            ("0011", "0000", 0),
                            ("0001", "0001", 1),
                            ("1111", "0101", 0.5),
                            ("0101", "1111", 0.5)]
        for bintext1, bintext2, expected_result in bintexts_results:
            fp1 = fps.bintext2fp(bintext1)
            fp2 = fps.bintext2fp(bintext2)
            self.assertAlmostEqual(fps.fp_tanimoto_similarity(fp1, fp2),
                                   expected_result)
            self.assertAlmostEqual(fps.fp_tanimoto_similarity(fp1, fp2),
                                   fps.fp_tanimoto_similarity(fp2, fp1))

    def test_fps2avg_tanimoto(self):
        bintexts = ("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
                    "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
        avg_tani = fps.fps2avg_tanimoto([fps.bintext2fp(bintext) for bintext in bintexts])
        self.assertAlmostEqual(avg_tani, 0.2916666667)

    def test_fps2avg_on_bits(self):
        bintexts = ("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
                    "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
        avg_on_bits = fps.fps2avg_on_bits([fps.bintext2fp(bintext) for bintext in bintexts])
        self.assertAlmostEqual(avg_on_bits, 2.0)

    def test_fps2avg_bit_value(self):
        bintexts = ("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
                    "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
        avg_bit_value = fps.fps2avg_bit_value([fps.bintext2fp(bintext)
                                               for bintext in bintexts])
        self.assertAlmostEqual(avg_bit_value, 0.5)

if __name__ == '__main__':
    unittest.main()
