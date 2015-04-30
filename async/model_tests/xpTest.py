import unittest

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


def suite():
    suite = unittest.TestSuite()
    suite.addTest(XpertModelLoad('test_A_XpertModuleLoaded'))
    suite.addTest(XpertModelLoad('test_B_HasRunMember'))
    suite.addTest(XpertModelLoad('test_C_HasRunFunc'))

    return suite

if __name__=='__main__':
    unittest.main()
