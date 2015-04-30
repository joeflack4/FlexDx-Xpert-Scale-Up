import unittest
import json
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
    from async import xpert_bg_inter
except:
    xpert_bg_inter = False
    pass


class XpertModelRun(unittest.TestCase):

    def setUp(self):
        self.data = {'sdgxp_cost': 30.0, 
                     'target_mdr': 3.7, 
                     'target_inc': 250.0, 
                     'outpt_cost': 10.0, 
                     'sm_cost': 2.0, 
                     'drug3_cost': 5000.0, 
                     'drug2_cost': 1000.0, 
                     'drug1_cost': 500.0, 
                     'target_hiv': 0.83, 
                     'homebrew': False, 
                     'int_select': 9, 
                     'gxp_cost': 15.0, 
                     'filename' : 'ModelRunTest.json'}
        with open(self.data['filename'], 'w') as fp:
            json.dump(self.data, fp)    

    def tearDown(self):
        #Maybe delete the json file?  For now nothing
        pass

    def test_A_TestExecution(self):
        #Does the model run?
        xpert_bg_inter.run( json.dumps(self.data) )
        with open(self.data['filename'], 'r') as fp:
            data = json.load(fp)

        self.assertEqual(data['progress'], 8, 'Model ran all 9 strategies')


if __name__=='__main__':
    unittest.main()
