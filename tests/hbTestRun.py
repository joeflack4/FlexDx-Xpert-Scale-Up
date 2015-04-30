import unittest
import json
import sys
import os
#This test file is an example of running the homebrew
#model using setUp and tearDown functions that run 
#before and after each test.  Tests are identified 
#as functions starting with 'test'

#This is adding the parent directory to the python
#search path to allow the import below to work with
#the current model code.
sys.path.append (os.path.split(os.getcwd())[0])

try:
    from async import homebrew_model
except:
    homebrew_model = False

class HoBrModelRun(unittest.TestCase):

    def setUp(self):
 
        self.data = {"filename" : "HBModelRunTest.json",
                     "ip": "127.0.0.1", 
                     "ud_tests": {}, 
                     "tests": ["Smear"], 
                     "diag": {
                        "hiv-ret": ["Smear", "None", 0.25, 0.05], 
                        "hiv+new": ["Smear", "None", 0.25, 0.05], 
                        "hiv+ret": ["Smear", "None", 0.25, 0.05], 
                        "hiv-new": ["Smear", "None", 0.25, 0.05]
                     }, 
                     "dst": {
                         "hiv-ret": ["DST_None", "DST_None"], 
                         "hiv+new": ["DST_None", "DST_None"], 
                         "hiv+ret": ["DST_None", "DST_None"], 
                         "hiv-new": ["DST_None", "DST_None"]
                     }, 
                     "model_inputs": {
                         "sdgxp_cost": "30", 
                         "hiv": "0.83", 
                         "outpt_cost": "10", 
                         "sm_cost": "2", 
                         "drug3_cost": "5000", 
                         "drug2_cost": "1000", 
                         "drug1_cost": "500", 
                         "mdr": "3.7", 
                         "gxp_cost": "15", 
                         "inc": "250"
                     }, 
                     "ud_strat_name": "User Defined",
                     "homebrew": True}
        with open(self.data['filename'], 'w') as fp:
            json.dump(self.data, fp)    

    def tearDown(self):
        #Maybe delete the json file?  For now nothing
        pass

    def test_A_TestExecution(self):
        #Does the model run?
        homebrew_model.run( self.data )
        with open(self.data['filename'], 'r') as fp:
            data = json.load(fp)

        self.assertEqual(data['progress'], 8, 'Homebrew ran all 9 strategies')


if __name__=='__main__':
    unittest.main()
