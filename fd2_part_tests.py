import fd2_functions
import numpy as np
import unittest
import fd2_part_functions as fd2part


class TestParticleFunctions(unittest.testcase):
  def testExit(self):
    v1 = 1.0
    v2 = 1.00000006
    v3 = 2.0
    v4 = -3.0
    v5 = -2.0
    v6 = -2.000000006

# types of assertions
# self.assertEqual(x,2)
# self.assertAlmostEqual(x,1./5.)
# self.assertRaises(Exception)
#   only works if there is a function inside that does the following
#     raise Exception ('ya messed up')

if __name__ = "__main__":
  unittest.main()

