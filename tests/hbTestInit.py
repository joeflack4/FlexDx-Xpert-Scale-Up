import unittest
import sys
import os
#This test file is occupied mainly with the existence
#of the homebrew model and the run() function that is used
#to execute it

#This is adding the parent directory to the python
#search path to allow the import below to work with
#the current model code.
sys.path.append (os.path.split(os.getcwd())[0])

try:
    from async import homebrew_model
except:
    homebrew_model = False


class HoBrModelInit(unittest.TestCase):

    def test_A_HoBrModuleLoaded(self):
        #Did the model load?
        self.assertTrue(homebrew_model, 'Module Homebrew: not loaded')

    def test_B_HoBrHasRunMember(self):
        #Does it have a .run member?
        self.assertTrue(hasattr(homebrew_model,'run'), 'Module Homebrew: no member run')

    def test_C_HoBrHasRunFunc(self):
        #Is that .run member a function?
        def f():
            pass
        if self.assertTrue(hasattr(homebrew_model,'run'), 'Module Homebrew: no member run'):
            self.assertEqual( type (f), type (homebrew_model.run), 'Module Homebrew: run is not a function' )


if __name__=='__main__':
    unittest.main()
