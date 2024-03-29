import fd2_functions as fd2 
import numpy as np
import unittest
import fd2_part_functions as pf


class TestParticleFunctions(unittest.TestCase):
  def testExit(self):
    v1 = 1.0
    v2 = 1.00000006
    v3 = 2.0
    v4 = -3.0
    v5 = -2.0
    v6 = -2.000000006


    c1 = pf.CheckExit(v1, v2, v3)
    self.assertEqual(c1, 1 )

    c2 = pf.CheckExit(v4, v5, v6)
    self.assertEqual(c2, 0 )

    c3 = pf.CheckExit(v1, v4, v2)
    self.assertEqual(c3, "flag")

    c4 = pf.CheckExit(v4, v1, v3)
    self.assertEqual(c4, 1)

    c5 = pf.CheckExit(v4, v1, v5)
    self.assertEqual(c5, 0)

    c6 = pf.CheckExit(0.0000005, v3, v1)
    self.assertEqual(c6,1)

    c7 = pf.CheckExit(0.0000005, v4, v1)
    self.assertEqual(c7, "flag")

    c8 = pf.CheckExit(v3 , 0.0000005, v1)
    self.assertEqual(c8, "flag")

    c9 = pf.CheckExit(v4, 0.0000005, v1)
    self.assertEqual(c9, 0)

  def testTravelTime(self):
    v = []
    v[0] = 0.005
    v[1] = 0.0051234515
    v[2] = -0.005
    v[3] = 2.4
    v[4] = 1.3
    v[5] = -3.

    A = []
    
    

# types of assertions
# self.assertEqual(x,2)
# self.assertAlmostEqual(x,1./5.)
# self.assertRaises(Exception)
#   only works if there is a function inside that does the following
#     raise Exception ('ya messed up')

if __name__ == "__main__":
  print " \n"
  suite = unittest.TestLoader().loadTestsFromTestCase(TestParticleFunctions)
  unittest.TextTestRunner(verbosity=2).run(suite)

