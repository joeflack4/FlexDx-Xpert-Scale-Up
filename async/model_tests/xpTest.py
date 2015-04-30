import unittest
import sys
import os
#This test file is occupied mainly with the existence
#of the xpert model and the run() function that is used
#to execute it

#This is adding the parent directory to the python
#search path to allow the import below to work with
#the current model code.
sys.path.append (os.path.split(os.getcwd())[0])

try:
    import xpert_bg_inter
except:
    xpert_bg_inter = False
    pass


class XpertModelLoad(unittest.TestCase):

    def test_A_XpertModuleLoaded(self):
        #Did the model load?
        self.assertTrue(xpert_bg_inter)

    def test_B_HasRunMember(self):
        #Does it have a .run member?
        self.assertTrue(hasattr(xpert_bg_inter,'run'))

    def test_C_HasRunFunc(self):
        #Is that .run member a function?
        def f():
            pass

        self.assertEqual( type (f), type (xpert_bg_inter.run) )


if __name__=='__main__':
    unittest.main()
