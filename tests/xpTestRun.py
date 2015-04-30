import unittest
import json
import sys
import os
#This test file is an example of running the xpert
#model using setUp and tearDown functions that run 
#before and after each test.  Tests are identified 
#as functions starting with 'test'
#Two Cases: First with run(), second with command
#line



#This is adding the parent directory to the python
#search path to allow the import below to work with
#the current model code.
sys.path.append (os.path.split(os.getcwd())[0])

try:
    from async import xpert_bg_inter
except:
    xpert_bg_inter = False


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
                     'filename' : 'XpertModelRunTest.json'}
        with open(self.data['filename'], 'w') as fp:
            json.dump(self.data, fp)    

    def tearDown(self):
        #Maybe delete the json file?  For now nothing
        pass


class XpertModelRunFunc(XpertModelRun):
    #Does the model run?

    def test_A_asFunction(self):
    
        xpert_bg_inter.run( self.data )
    
        with open(self.data['filename'], 'r') as fp:
            data = json.load(fp)

        self.assertTrue(data, 'No data was loaded from the file as Function')
    
    def test_B_asCommandLine(self):
    
        os.system("../async/xpert_bg_inter.py {}".format(self.data['filename']))
    
        with open(self.data['filename'], 'r') as fp:
            data = json.load(fp)

        self.assertTrue(data, 'No data was loaded from the file as CommandLine')

if __name__=='__main__':
    unittest.main()
