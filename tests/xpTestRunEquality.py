import unittest
import json
import sys
import os
#This test file is an example of running the xpert
#model using setUp and tearDown functions that run 
#before and after each test.  Tests are identified 
#as functions starting with 'test'
#One Case:  Equality of value between run() and
#command line.  Runs the model twice


#This is adding the parent directory to the python
#search path to allow the import below to work with
#the current model code.
sys.path.append (os.path.split(os.getcwd())[0])

try:
    from async import xpert_bg_inter
except:
    xpert_bg_inter = False


class XpertModelRunEqual(unittest.TestCase):

    def setUp(self):
        self.filenames = ['XPEqualTest1.json','XPEqualTest2.json']

        self.data = {'filename' : self.filenames[0],
                     'sdgxp_cost': 30.0, 
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
                     'gxp_cost': 15.0}
        with open(self.data['filename'], 'w') as fp:
            json.dump(self.data, fp)    

        self.data['filename'] = self.filenames[1]

        with open(self.data['filename'], 'w') as fp:
            json.dump(self.data, fp)

    def tearDown(self):
        #Maybe delete the json file?  For now nothing
        pass

    def test_A_Equality(self):

        xpert_bg_inter.run( self.data ) #self.data should still have self.filenames[1].

        os.system("../async/xpert_bg_inter.py {}".format(self.filenames[0]))
        
        with open(self.filenames[0], 'r') as fp:
            data1 = json.load(fp)
            del data1['filename']
        with open(self.filenames[1], 'r') as fp:
            data2 = json.load(fp)
            del data2['filename']

        self.assertTrue(data1 == data2, 'Command line results were not equal to run() results')

if __name__=='__main__':
    unittest.main()
