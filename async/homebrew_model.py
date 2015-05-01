#!/usr/bin/env python
 
#############################################################
# "TB Diagnostic Selection" Model
# David Dowdy, Jason Andrews, Peter Dodd, Robert Gilman
#############################################################
#############################################################
# 1. MODEL DESCRIPTION
#############################################################
# This model has the following compartments (see appendix):
# 0. uninfected
# 1-3. latently infected, DS/INH/MDR
# 4-6. "early active", DS/INH/MDR - lasts 9 mos
# 7-12. active, seeking dx:
#       7-8: DS, smear-pos/smear-neg
#       9-10: INH, smear-pos/smear-neg
#       11-12: MDR, smear-pos/smear-neg
# 13-18. active, dx in progress, same categories as 7-12
# 19-24. active, on inappropriate tx, same categories as 7-12
# 25-49. same as above, only HIV-infected
# 50-99. same as all the above, previously treated for TB
############################################################
############################################################
# 2. IMPORT PACKAGES
############################################################
import numpy as np
import scipy as sp
from scipy.integrate import odeint
from scipy.optimize import fsolve
from copy import deepcopy as dcpy
from itertools import *
import json
#JP MOD
import sys #Command Line arguments
import os #Environment Variables
from datetime import datetime #for logging

# ############################################################
# _                          _
# | |__   ___  _ __ ___   ___| |__  _ __ _____      __
# | '_ \ / _ \| '_ ` _ \ / _ \ '_ \| '__/ _ \ \ /\ / /
# | | | | (_) | | | | | |  __/ |_) | | |  __/\ V  V /
# |_| |_|\___/|_| |_| |_|\___|_.__/|_|  \___| \_/\_/

###############################################################
#BEGIN JPMOD
###############################################################


###############################################################
# DEFINE PARAMETERS
###############################################################

class Gl_Vars:
    def __init__(self):
        self.ip_addr = 'Unknown'
        self.ud_strat_name = 'User Defined'
        self.intopt = 0
        self.int_select = 9
        #Inputs for later
        self.target_inc = 0.0#float(jdata['model_inputs']['inc']) # 250.
        self.target_mdr = 0.0#float(jdata['model_inputs']['mdr']) / 100.0 #3.7 / 100
        self.target_hiv = 0.0#float(jdata['model_inputs']['hiv']) # 0.83
        self.drug1_cost = 0.0#float(jdata['model_inputs']['drug1_cost']) #500.
        self.drug2_cost = 0.0#float(jdata['model_inputs']['drug2_cost']) #1000.
        self.drug3_cost = 0.0#float(jdata['model_inputs']['drug3_cost']) #5000.
        self.outpt_cost = 0.0#float(jdata['model_inputs']['outpt_cost']) #10.
        self.sm_cost = 0.0#float(jdata['model_inputs']['sm_cost']) #2.
        self.gxp_cost = 0.0#float(jdata['model_inputs']['gxp_cost']) #15.
        self.sdgxp_cost = 0.0#float(jdata['model_inputs']['sdgxp_cost']) - gxp_cost #30. - gxp_cost

        # Parameters that will later be fit to user inputs
        # The actual values inputted here are NOT used by the program.
        self.beta_s = 8.0             # transmission rate
        self.fit_m = 0.707            # relative fitness of MDR-TB
        self.fit_i = 1.-(1.-self.fit_m)*0.25 # relative fitness of INHr-TB
        self.hiv_inc = 0.00064        # HIV incidence

        # quantities that will be needed for later calculation.
        self.force = 0.               # force of infection, initiate at zero
        self.tb_num = 100.            # number of compartments = 100
        self.mortv = np.zeros(self.tb_num) # mortality vector (to maintain stable pop)

        # Keep the population size at 100,000 for understanding of results.
        self.pop_size = 100000.       # population size = 100,000 (for std rates)

        ###################################################################################
        ### The following parameters can all be modified according to user preferences: ###
        ###################################################################################

        # in emerging-MDR scenario, change "target_mdr" below to the desired
        # MDR prevalence among new cases at the END of 5 years
        #   (input value will define the prevalence at the beginning of 5 years)
        # INPUT AS A PROPORTION (e.g., 0.1 = 10% prevalence)
        # Thus, for example, if you want a scenario that goes from MDR prevalence
        #   of 4.0% to 6.0% over 5 years, input "4.0" when asked in the interface,
        #   and change the line below to read, "target_mdr2 = 0.06".

        self.react = 0.0005           # HIV-negative reactivation rate
        self.react_h = 0.05           # HIV-positive reactivation rate
        self.relbeta_sn = 0.22        # relative infectiousness of smear-neg TB, incl extrapulm
        self.rapid = 0.14             # proportion of infections w/ primary progression
        self.rapid_h = 0.47           # primary progression proportion, HIV+
        self.prot = 0.79              # reduction in rate of reinfection if LTBI & HIV-neg
        self.cure_sn = 0.267          # "self-cure" rate, smear-neg (20% CFR at 3 yrs)
        self.cure_sp = 0.1            # self-cure rate, smear-pos (70% CFR at 3 yrs)
        self.mort = 1.0/55.           # background mortality (life exp = 70)
        self.mort_h = 0.05            # mortality/prevalence of HIV
        self.mort_sn = 0.0667         # mortality rate, smear-negative, 20% CFR at 3 yrs
        self.mort_sp = 0.233          # mortality rate, smear-positive, 70% CFR at 3 yrs
        self.mort_tbh = 1.0           # assume 1-year duration of TB before uniform death if HIV+
        self.prop_inf = 0.63          # proportion of incident cases that are smear-pos
        self.prop_inf_h = 0.50        # proportion smear-pos if HIV-pos

        self.fail_s = 0.04            # probability of failure or relapse, DS-TB
        self.fail_rt = self.fail_s         # probability of failure or relapse, retreatment DS-TB
        self.fail_i1 = 0.21           # prob of failure or relapse, INHr treated with 1st-line
        self.fail_i2 = 0.16           # prob of failure or relapse, INHr treated with streptomycin
        self.fail_m1 = 0.50           # prob of failure or relapse, MDR-TB treated with 1st-line
        self.fail_m2 = 28./93.        # prob of failure or relapse, MDR-TB treated w/ 2nd-line
        self.predx_dur = 0.75         # duration of infectiousness before seeking care, HIV-
        self.predx_dur_h = 1./12      # duration of infectiousness before seeking care, HIV+
        self.acq_s = 0.003/self.fail_s     # proportion of failed treatments acquiring resistance, DS-TB
        self.acq_mdr = 0.33           # proportion of acq resi among DS-TB that is MDR (vs. INH-mono)
        self.acq_i1 = 0.045/self.fail_i1   # proportion of failed treatments acquiring resistance, INHr tx w/ 1st-line
        self.acq_i2 = 0.017/self.fail_i2   # proportion of failed treatments acquiring resistance, INHr tx w/ 2nd-line
        self.defvfail_s = 6./7.       # proportion of inappropriate tx's that are default (vs. failure), DS-TB
        self.defvfail_i = 2./3.       # proportion of inapp tx's that are default, INHr
        self.defvfail_m = 11./25.     # proportion of inapp tx's that are default, MDR

        self.dx_rate = 5.             # diagnostic rate: 5 diagnostic attempts per year
        self.ltfu = 0.15**1              # 15% initial default proportion, smear/GXP
        self.ltfu_both = 0.2          # 25% initial default proportion, smear then GXP
        self.emp_tx = 0.25**6            # 25% of pts treated empirically - SET ZERO
        self.cx_sens = 0.85           # Sensitivity of single culture
        self.gxp_sens = 0.72          # Sensitivity of single GXP
        self.sm_spec = 0.98           # specificity of smear
        self.cx_spec = 0.98           # specificity of culture
        self.gxp_spec = 0.98          # specificity of GXP
        self.cxr_sens = 0.98          # sensitivity of MODS for RIF resistance
        self.cxi_sens = self.cxr_sens      # sensitivity of MODS for INH resistance
        self.cxr_spec = 0.994         # specificity of MODS for RIF
        self.cxi_spec = 0.958         # specificity of MODS for INH
        self.gxr_sens = 0.944         # sensitivity of GXP for RIF
        self.gxr_spec = 0.983         # specificity of GXP for RIF
        self.t_smgxp = 7.0/365        # delay from smear/GXP result to tx: 7 days
        self.t_both = 14.0/365        # delay from smear + GXP: 14 days
        self.t_cxr = 30.0/365         # time to culture-based DST is the same (30 days)

# LPA and chest X-ray parameters
        self.t_lpa = 20.0/365         # time to DST using LPA, per Boehme Lancet
        self.lpa_sens = .97           # from Medscape http://www.medscape.com/viewarticle/740253_9
        self.lpa_spec = .99           # from Medscape http://www.medscape.com/viewarticle/740253_9
        self.lpa_cost = 9.4           # from FIND http://www.finddiagnostics.org/about/what_we_do/successes/find-negotiated-prices/mtbdrplus.html
        self.t_rad = 1.0/365                # time to chest X-ray result
        self.rad_sens = 0.978         # X-ray van t'Hoog's review http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
        self.rad_spec = 0.754         # X-ray van t'Hoog's review http://www.who.int/tb/Review2Accuracyofscreeningtests.pdf
        self.rad_cost = 10            # bit of a guess
        self.t_lam = 1.0/365          #quick
        self.lam_sens = 0.5
        self.lam_spec = .999
        self.lam_cost = 3.5           #http://www.aidsmap.com/Determine-LAM-urine-antigen-TB-test-is-highly-cost-effective-for-use-in-hospitalised-people-living-with-HIV/page/2455261/


        self.dur_fail = 6.            # number of months until failing cases are re-evaluated
        self.prop_cough = 0.01**1        # proportion of people w/ non-TB cough eval'd per yr - SET ZERO

        self.Te = self.t_smgxp             # empirical tx delay -- added pjd
##################################
### End modifiable parameters. ###
##################################

g = Gl_Vars()
#HELPER FUNCTIONS (JPMOD)

FLEXDX_HOME = os.environ.get('FLEXDX_HOME')

if FLEXDX_HOME:
    prod_dev_path = "{}/".format(FLEXDX_HOME)
else:
    prod_dev_path = "{}/".format(os.getcwd())

def json_write( filename, data ):
    with open (filename, 'w+b') as fp:
        json.dump(data, fp, separators=(',',':'))

def hb_log ( message ):
    with open (prod_dev_path + "homebrew.log", 'a') as fp:
        fp.write( '{} - ({}): {}\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),g.ip_addr, message ) )

hb_log("------- Execution Begins -------")

hb_log('Past Constant Definition')


#stderr redirection for debugging
sys.stderr = open(prod_dev_path + 'homebrew_python.err','a')

with open (prod_dev_path + 'homebrew_python.err', 'a') as fp:
    fp.write( '---- New run at : {} ----\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) )

###############################################################

class Diagnostic:
    '''This is a class descibing a single diagnostic, at the moment for inputs that are
    TB-, TB+DR-, TB+HR, TB+RR and with outputs TB-,TB+DR-,TB+DR+. These can be chained together by specifying \'next\'
    and the getTables() method calculates the transition matrix for the subtree rooted at a given
    diagnostic, as well as the average costs and delayed for each transition. '''

    def __init__(self,sens,spec,cost,delay,ltfu,label='',DRsens=0.0,DRspec=1.0):
        '''This constructs an instance of a diagnostic based on data
        on the sensitivity (smr+/- as array) and specificity for Mtb , the cost and delay incurred,
        and (if specified) the spensitivity and specificity for DR given a
        Mtb positive result. DST can be omitted in which case, nothing turns up TB+DR+.'''

        if type(sens)==float:   # duplicate if only 1 given
            sens = [sens,sens]
        self.sens = sens; self.spec = spec; self.cost = cost;self.delay = delay;self.ltfu = ltfu;
        self.DRsens = DRsens;self.DRspec = DRspec
        self.label = label
        self.root = True        # by default, this is a root

        # smr-/+ x an array (TB-,TB+DR-, TB+HR+, TB+RR+) x (-,+n,+y)
        self.transition = np.zeros((2,4,3))
        for i in range(2):      # smr- then smr+
            self.transition[i,:,:] = np.array([[ spec    ,    (1-spec)       , 0                 ],
                                               [1-sens[i], sens[i]*DRspec    , sens[i]*(1-DRspec)],
                                               [1-sens[i], sens[i]*DRspec    , sens[i]*(1-DRspec)],
                                               [1-sens[i], sens[i]*(1-DRsens), sens[i]*DRsens]   ] )

        # 0=no tx, 1=1st line, 2=2nd line = TB-, TB+DRn, TB+DRy
        self.next = [0,1,2]     # default as treatments, can be next diagnostics
        self.cost = cost
        self.delay = delay

    def __str__(self):
        '''For display purposes'''
        print self.transition
        return "Diagnostic instance %s, for specified h x p. Data on characteristics and next test. Use help(instance) to learn more" % self.label

    def setnext(self,k,newdx):
        '''This is for setting the next diagnostic. Same as naively changing vector next, but takes care of root flag at the same time.'''
        if k not in [0,1,2]:
            print 'Error! k must be in 0..2'
        self.next[k] = dcpy(newdx)    # add to tree, as copy
        self.next[k].root = False      # record that this is no longer a root

    def getTables(self):
        '''To get the overall table of (TB-,TB+DR-,TB+H+,TB+R+) x (no tx, 1st, 2nd), and similarly for costs and delays.
        It returns two sets of these: one for smr- entries, one for smr+ entries, meaning that the false-positive calculations are duplicated.
        Computes transitions as T_{ij} = \sum_k A_{ik}B^k_{ij}recursively.'''
        txo = np.zeros((2,4,3))   # transitions table
        cost = np.zeros((2,4,3))  # cost for in state getting to outcome etc.
        delay = np.zeros((2,4,3)) # delay
        if self.root:
            cost += self.cost * self.transition # add on own costs if root
        for k  in [0,1,2]:
            if isinstance(self.next[k],Diagnostic):
                txon, costn, delayn = self.next[k].getTables()
                for j in [0,1,2]:
                    nextbit = self.transition[:,:,k] * txon[:,:,j]
                    txo[:,:,j] += (1-self.ltfu) * nextbit
                    cost[:,:,j] += self.next[k].cost * (1-self.ltfu) * nextbit # new costs lookahead
                    delay[:,:,j] += self.delay * (1-self.ltfu) * nextbit
            elif self.next[k] == k:
                txo[:,:,k] += (1-self.ltfu) * self.transition[:,:,k]
                cost[:,:,k] += 0 # how do we include root costs?
                delay[:,:,k] += self.delay *  (1-self.ltfu) * self.transition[:,:,k]
            else:
                print "Invalid node in next for getTables! Need next[k]==k if not a subsequent test! label=%s, node=%j" % (self.label,j)
        return txo, cost, delay



class Algorithm:
    '''This is for working with a set of 4 trees associated with a Diagnostic - 1 for each HIV x retx possibility'''

    def __init__(self,HnRn,HnRp,HpRn,HpRp,label=''):
        '''Construct an instance from array of Diagnostics (trees) for each of HIV+/-, retx +/-
        '''
        self.trees = [[dcpy(HnRn),dcpy(HnRp)],[dcpy(HpRn),dcpy(HpRp)]] # deep copies made here now
        self.label = label

    def __str__(self):
        '''For display purposes'''
        return "Algorithm instance %s: 4 trees for  h x p.  Use help(instance) to learn more" % self.label

    def getTables(self):
        '''Gather relevant data as for trees and return - HIV x retx x smr x (TB in) x (outcome)'''
        txo = np.zeros((2,2,2,4,3))   # transitions table
        cost = np.ones((2,2,2,4,3))
        delay = np.ones((2,2,2,4,3))
        for i in range(2):
            for j in range(2):
                T,C,D = self.trees[i][j].getTables()
                txo[i,j,:,:] = T
                cost[i,j,:,:] = C
                delay[i,j,:,:] = D

        return txo, cost, delay # add next stage to convert these to necessaries for txrate etc

    def getFuns(self):
        '''This returns functions txrate, invspec, and costdx - wrapping the information in the tables to interface with David\'s code.
        Use as: \"new_txrate, new_invspec, new_costdx = myAlg.getFuns()\", where the new_functions are as their namesakes, but without the intervention argument.
        '''
        #Global Vars used locally
        drug1_cost = g.drug1_cost
        drug2_cost = g.drug2_cost
        drug3_cost = g.drug3_cost
        fail_s = g.fail_s
        fail_rt = g.fail_rt
        fail_i1 = g.fail_i1
        fail_i2 = g.fail_i2
        fail_m1 = g.fail_m1
        fail_m2 = g.fail_m2
        dx_rate = g.dx_rate

        # preliminary reshaping of data
        # treatment costs
        c1 = np.zeros((2,2,3,2))                  # cost of 1st line tx
        c1[:,0,:,:] = drug1_cost; c1[:,1,:,:] = drug2_cost
        c2 = drug3_cost * np.ones((2,2,3,2))                  # cost of 2nd line tx
        # treatment success probabilities
        p1suc = np.zeros((2,2,3,2))               # probabilities of 1st line tx success
        p1suc[:,0,0,:] = 1-fail_s; p1suc[:,1,0,:] = 1-fail_rt;   # DS
        p1suc[:,0,1,:] = 1-fail_i1; p1suc[:,1,1,:] = 1-fail_i2; p1suc[:,:,2,:] = 1-fail_m1;    # DR
        p2suc = np.zeros((2,2,3,2))               # probabilities of 2nd line tx success
        p2suc[:,0,0,:] = 1-fail_s; p2suc[:,1,0,:] = 1-fail_rt   # DS
        p2suc[:,:,1,:] = 1-fail_i2; p2suc[:,:,2,:] = 1-fail_m2   # DR

        # get tree data
        txo, cst, dly = self.getTables()   # calculate the relevant data  HPI (tb3)

        # data needed: hpdi
        p1 = np.zeros((2,2,3,2))     #
        p2 = np.zeros((2,2,3,2))     #
        pe = np.zeros((2,2,3,2))     #
        txdata = np.zeros((2,2,3,2))     #
        txidata = np.zeros((2,2,3,2))     #
        cdata = np.zeros((2,2,3,2))     #
        cidata = np.zeros((2,2,3,2))     #
        dlydata = np.zeros((2,2,3,2))     #
        cdxdata = np.zeros((2,2,3,2))     #
        for k in product(range(2),range(2),range(3),range(2)):   # loop over hpdi - d has different meaning
            h,p,d,i = k                                          # hpdi=(HIV, previous tx, dr, smr)
            dd = d + 1                                           # different d in trees, skips 'no treatment'
            p1[h,p,d,i] = txo[h,p,i,dd,1]                        # prob 1st line tx
            p2[h,p,d,i] = txo[h,p,i,dd,2]                        # prob 2nd line tx
            pe[h,p,d,i] = g.emp_tx * (1-p1[h,p,d,i]-p2[h,p,d,i])   # empiric treatment for those not treated
            dlydata[h,p,d,i] = (dly[h,p,i,dd,1] + dly[h,p,i,dd,2] + g.Te*pe[h,p,d,i])/ (p1[h,p,d,i]+p2[h,p,d,i]+pe[h,p,d,i]+1e-9) # dly include probs, conditioned on treatment
            cdxdata[h,p,d,i] = cst[h,p,i,dd,:].sum()   # cost of dx, sum as assuming ans weighted by probs 1,2 as tx

        # emptx either applied up front with no ltfu or as subsequent test with ltfu
        # for dxcost
        dxcdata = np.zeros((2,2))         # TB-ve diagnosis costs
        for h in range(2):
            for p in range(2):
                dxcdata[h,p] = cst[h,p,0,0,:].sum()   # mean cost of dx - needs smr 0 as data for TB- duplicates here
                dxcdata[h,p] += txo[h,p,0,0,1]*c1[h,p,0,0] + txo[h,p,0,0,2]*c2[h,p,0,0] # mean cost of treatment - see defn of c1/2 above - switch on p
                # 0,0 = smr-ve,TB-ve

        # for invspec
        ispdata = np.zeros((2,2))
        for h in range(2):
            for p in range(2):
                ispdata[h,p] =  1-txo[h,p,0,0,1:2].sum()   # 1-mean prob TB+, assuming smr-,retx-

        # further computations
        txdata = dx_rate * ((p1+pe)*p1suc + p2*p2suc)
        txidata = dx_rate * ((p1+pe)*(1-p1suc) + p2*(1-p2suc))
        cdata =  ((p1+pe)*p1suc*c1 + p2*p2suc*c2)/((p1+pe)*p1suc + p2*p2suc+1e-9)
        cidata = ((p1+pe)*(1-p1suc)*c1 + p2*(1-p2suc)*c2)/((p1+pe)*(1-p1suc) + p2*(1-p2suc)+1e-9)

        # invert
        dlydata = 1.0/(dlydata+1e-9)

        # txrate
        def TXR(smear,dr,retx,hiv):
            # and then calculates: tx, delay, tx_inapp, cost_dx, cost_tx, cost_inapp
            # 0. Rate of successful treatment (dx rate * prob of successful tx)
            # 1. Delay from sending diagnostic test to initiation of treatment -- NB really corresponding rate
            # 2. Rate of treatment that will fail (dx rate * prob of unsuccessful tx)
            # 3. Cost of diagnosis
            # 4. Cost of treatment if successful
            # 5. Cost of treatment if unsuccessful
            tx = txdata[hiv,retx,dr,smear]
            delay = dlydata[hiv,retx,dr,smear]   # already inverted above
            tx_inapp = txidata[hiv,retx,dr,smear]
            cost_dx = cdxdata[hiv,retx,dr,smear]
            cost_tx = cdata[hiv,retx,dr,smear]
            cost_inapp = cidata[hiv,retx,dr,smear]
            return tx, delay, tx_inapp, cost_dx, cost_tx, cost_inapp
        # invspec
        def INVS(hiv,retx):
            # spec of TB- by HIV
            sp = ispdata[hiv,retx]
            return sp
        # costdx
        def CDX(retx, hiv):
            # cost of diagnosing & treating TB-negatives who present with symptoms,
            cost = dxcdata[hiv,retx]
            return cost
        # return these functions
        return TXR, INVS, CDX

def run (jdata): #Comes here as python dict()
#JP MOD removed the inputs, set to model inputs

    g.target_inc = float(jdata['model_inputs']['inc']) # 250.
    g.target_mdr = float(jdata['model_inputs']['mdr']) / 100.0 #3.7 / 100
    g.target_hiv = float(jdata['model_inputs']['hiv']) # 0.83
    g.drug1_cost = float(jdata['model_inputs']['drug1_cost']) #500.
    g.drug2_cost = float(jdata['model_inputs']['drug2_cost']) #1000.
    g.drug3_cost = float(jdata['model_inputs']['drug3_cost']) #5000.
    g.outpt_cost = float(jdata['model_inputs']['outpt_cost']) #10.
    g.sm_cost = float(jdata['model_inputs']['sm_cost']) #2.
    g.gxp_cost = float(jdata['model_inputs']['gxp_cost']) #15.
    g.sdgxp_cost = float(jdata['model_inputs']['sdgxp_cost']) - g.gxp_cost #30. - gxp_cost

    try:
        g.ud_strat_name = jdata['ud_strat_name']
    except:
        pass
    else:
        del jdata['ud_strat_name']

    try:
        g.ip_addr = jdata['ip'];
    except:
        pass
    else:
        del jdata['ip']

    tmp_filename = jdata['filename']

    hb_log("Write Target Filename : {}".format(tmp_filename))

    outdata = {}

    outdata['int_select'] = g.int_select
    outdata['target_inc'] = g.target_inc
    outdata['target_mdr'] = g.target_mdr * 100
    outdata['target_hiv'] = g.target_hiv
    outdata['drug1_cost'] = g.drug1_cost
    outdata['drug2_cost'] = g.drug2_cost
    outdata['drug3_cost'] = g.drug3_cost
    outdata['outpt_cost'] = g.outpt_cost
    outdata['sm_cost'] = g.sm_cost
    outdata['gxp_cost'] = g.gxp_cost
    outdata['sdgxp_cost'] = g.sdgxp_cost + g.gxp_cost# - gxp_cost
    outdata['progress'] = -1

    hb_log("Deleting old json data")

#Remove the contents of the temp file
    open(tmp_filename, 'w').close()

    hb_log ("Writing the new data")
#Switch the JSON over to write back to the page
    json_write( tmp_filename, outdata )
     
###############################################################
#END JPMOD
###############################################################

# In low-incidence scenario, fit initially to target incidence of 50:
    g.target_inc2 = g.target_inc
    if g.target_inc2<50:
        g.target_inc=50.

    g.target_mdr2 = g.target_mdr


#--------------------------------------  this bit defines the interventions ----------------------------------
#                            _                      _
#   _____  ___ __   ___ _ __(_)_ __ ___   ___ _ __ | |_ ___
#  / _ \ \/ / '_ \ / _ \ '__| | '_ ` _ \ / _ \ '_ \| __/ __|
# |  __/>  <| |_) |  __/ |  | | | | | | |  __/ | | | |_\__ \
#  \___/_/\_\ .__/ \___|_|  |_|_| |_| |_|\___|_| |_|\__|___/
#           |_|


# basic tests
    test_smr = Diagnostic([0.,1.0],g.sm_spec,(g.sm_cost+g.outpt_cost),g.t_smgxp,g.ltfu,label='smr',DRsens=0,DRspec=1) # smear as a test
    test_sdgxp = Diagnostic([g.gxp_sens,1.0],g.gxp_spec,(g.gxp_cost+g.outpt_cost+g.sdgxp_cost),1.0/365,0.0,label='sdgxp',DRsens=g.gxr_sens,DRspec=g.gxr_spec) # sd GXP as a test
    test_gxp = Diagnostic([g.gxp_sens,1.0],g.gxp_spec,(g.gxp_cost+g.outpt_cost),g.t_smgxp,g.ltfu,label='gxp',DRsens=g.gxr_sens,DRspec=g.gxr_spec) # GXP as a test
    test_gxpC = Diagnostic([1.0,1.0],0,(g.gxp_cost),1e-9,0,label='gxpc',DRsens=g.gxr_sens,DRspec=g.gxr_spec) # GXP confirmatory

# some test combinations
# GXP for smear positive (e.g. retreatment)
    test_smrXp = dcpy(test_smr)    # make copy of smear
    test_smrXp.setnext(1,test_gxpC)# Xpert if smr +ve

# GXP for smear negative (e.g. hiv+)
    test_smrXn = dcpy(test_smr)    # make copy of smear
    test_smrXn.setnext(0,test_gxp)  # Xpert if smr -ve


# GXP fu regardless of smear status
    test_smrX = dcpy(test_smr)    # make copy of smear
    test_smrX.setnext(0,test_gxp)# Xpert if smr -ve
    test_smrX.setnext(1,test_gxpC)# Xpert if smr +ve, confirmatory

# myAlg = Algorithm(testnn,testnp,testpn,testpp)
    Alg0 = Algorithm(test_smr,test_smrXp,test_smr,test_smrXp) # "0 = baseline"
    Alg1 = Algorithm(test_smrXp,test_smrXp,test_smrXp,test_smrXp) # "1 = Xpert for smear-positive"
    Alg2 = Algorithm(test_smr,test_smrXp,test_gxp,test_gxp) # "2 = Xpert for HIV+"
    Alg3 = Algorithm(test_smr,test_gxp,test_smr,test_gxp) # "3 = Xpert for previously treated"
    Alg4 = Algorithm(test_smr,test_gxp,test_smrXn,test_gxp) # "4 = Xpert for sm-neg HIV+ or prev tx"
    Alg5 = Algorithm(test_smr,test_gxp,test_gxp,test_gxp) # "5 = Xpert for all HIV+ or prev tx"
    Alg6 = Algorithm(test_smrXn,test_smrX,test_smrXn,test_smrX) #  "6 = Xpert for smear-negative"
    Alg7 = Algorithm(test_gxp,test_gxp,test_gxp,test_gxp) # "7 = Xpert for all"
    Alg8 = Algorithm(test_sdgxp,test_sdgxp,test_sdgxp,test_sdgxp) # "8 = Xpert for all, same-day"

# Algz = [Alg0,Alg1,Alg2,Alg3,Alg4,Alg5,Alg6,Alg7,Alg8]
# AllFz = [Algz[i].getFuns() for i in range(9)] # get all the functions for later

    hb_log('Past Algorithm Definitions')

#               _       _                   _
# __      _____| |__   (_)_ __  _ __  _   _| |_
# \ \ /\ / / _ \ '_ \  | | '_ \| '_ \| | | | __|
#  \ V  V /  __/ |_) | | | | | | |_) | |_| | |_
#   \_/\_/ \___|_.__/  |_|_| |_| .__/ \__,_|\__|
#                              |_|


    def algfromdict(J):

        ptypes = ['hiv-new','hiv-ret','hiv+new','hiv+ret']   # note different order
        # check we have the relevant tests defined
        testlist = ['None','DST_None','Smear','GXP','DST_GXP', \
                    'E_hiv+new','E_hiv+ret','E_hiv-new','E_hiv-ret']         # these need default definitions
        # testlist += J['ud_tests'].keys()      #any user-defined tests
        ok = True
        for pt in ptypes:
            for i in [0,1]:
                if (J['diag'][pt][i] not in testlist) and (J['diag'][pt][i] not in J['tests']) and (J['diag'][pt][i] not in J['ud_tests'].keys()):
                    print 'Test type ' + J['diag'][pt][i] + ' not defined for test ' + str(i+1) +' for ' + pt
                    hb_log ('algfromdict: Test type ' + J['diag'][pt][i] + ' not defined for test ' + str(i+1) +' for ' + pt)
                    ok = False
                if (J['dst'][pt][i] not in testlist) and (J['dst'][pt][i] not in J['tests']) and (J['dst'][pt][i] not in J['ud_tests'].keys()):
                    print 'DST test type ' + J['dst'][pt][i] + ' not defined for test ' + str(i+1) +' for ' + pt
                    hb_log ('algfromdict: DST test type ' + J['dst'][pt][i] + ' not defined for test ' + str(i+1) +' for ' + pt)
                    ok = False
        if not ok:
            hb_log('algfromdict: BAILING 1')
            return 0
        hb_log("algfromdict: Before Test definition")
        # define the  tests
        testdict = {}
        # define known tests
        testdict['CXR'] = Diagnostic( [g.rad_sens,g.rad_sens], g.rad_spec, g.rad_cost, g.t_both, 0, label='Chest X-ray', DRsens=0, DRspec=1 ) # chest X-ray
        testdict['Smear'] = Diagnostic( [0,1], g.sm_spec, (g.sm_cost+g.outpt_cost), g.t_smgxp, g.ltfu, label='Smear', DRsens=0, DRspec=1 ) 
        testdict['GXP'] = Diagnostic( [g.gxp_sens,1], g.gxp_spec, (g.gxp_cost+g.outpt_cost), g.t_smgxp, g.ltfu, label='GXP', DRsens=g.gxr_sens, DRspec=g.gxr_spec )
        testdict['LAM'] = Diagnostic( [g.lam_sens,g.lam_sens], g.lam_spec, (g.outpt_cost + g.lam_cost), g.t_lam, 0, label='LAM', DRsens=0, DRspec=1 ) 
        testdict['DST_GXP'] = Diagnostic( [1,1], 0, g.gxp_cost, g.t_smgxp, 0, label='DST_GXP', DRsens=g.gxr_sens, DRspec=g.gxr_spec )
        testdict['DST_LPA'] = Diagnostic( [1,1], 0, g.lpa_cost, g.t_lpa, 0, label='DST_LPA', DRsens=g.lpa_sens, DRspec=g.lpa_spec )
        for test in J['diag'].keys():         #defining the empirical tests
            testdict['E_'+test] = Diagnostic( [J['diag'][test][2]]*2, 1-J['diag'][test][3], 0, 0, 0, label='E_'+test, DRsens=0, DRspec=1 )
        testdict['Culture'] = Diagnostic( [g.cx_sens,1], g.cx_spec, (g.gxp_cost+g.outpt_cost), g.t_cxr, g.ltfu, label='CX', DRsens=0, DRspec=1 )
        testdict['DST_Culture'] = Diagnostic( [1,1], 0, g.gxp_cost, g.t_cxr, g.ltfu, label='DST_CX', DRsens=g.cxr_sens, DRspec=g.cxr_spec )
        hb_log("algfromdict: After Test definition")
        K = J['ud_tests']
        for test in K:
            if test in testlist[2:]:          # check not a naming collision
                print 'Test ' + test + ' already exists! Bailing!'
                hb_log ('algfromdict: Test ' + test + ' already exists! Bailing 2!')
                return 0
            else:                             #test name OK
                # also: NB I haven't put in any checking that the relevant fields have been filled in with valid values
                # assuming that the delay on user defined DST is t_smgxp and 0 ltfu
                if test[0:3]=='DST':          #is it a DST test - NB this will break if someone calls a normal test DST*
                    testdict[test] = Diagnostic( [1,1], 0, K[test]['cost'], g.t_smgxp, 0, label=test, DRsens=K[test]['sens'], DRspec=K[test]['spec'])
                else:
                    testdict[test] = Diagnostic( [K[test]['senssmn'],K[test]['senssmp']], K[test]['spec'], K[test]['cost'], g.t_smgxp, K[test]['spec'], label=test, DRsens=0, DRspec=1)

        hb_log('algfromdict: After ud_tests')
        K = J['diag']
        # linking tests: test2 only done if test1 negative
        testused = [0]*4                      # test used for each patient type
        for i in range(4):
            hb_log('algfromdict: testused' + str(i) + ' Build')
            testused[i] = dcpy(testdict[ K[ptypes[i]][0] ])
            if K[ptypes[i]][1]!='None':            # there is a second test used
                hb_log('algfromdict: testused' + str(i) + ' Second Test')
                negpos = 1 if K[ptypes[i]][0]=='CXR' else 0 # 2nd test if 1st: +ve for CXR; -ve for others
                testused[i].setnext(negpos, testdict[ K[ptypes[i]][1] ] )   # set it (if first test -/+ve)
                testused[i].next[negpos].setnext(0, testdict[ 'E_'+ ptypes[i] ] )   # add on empiric treatment for -ve
                if J['dst'][ ptypes[i] ][1]!='DST_None':               #there is a DST for test 2+ve
                    hb_log('algfromdict: testused' + str(i) + ' a DST for test 2+ve')
                    testused[i].next[negpos].setnext(1, testdict[ J['dst'][ptypes[i]][1] ] )   # DST if test 1-/+, test2  +-
                    testused[i].next[negpos].setnext(2, testdict[ J['dst'][ptypes[i]][1] ] )   # DST if test 1-/+, test2  ++
            else:                                                      #no second test
                hb_log('algfromdict: testused' + str(i) + ' One Test')
                testused[i].setnext(0, testdict[ 'E_'+ ptypes[i] ] )   # add on empiric treatment if test 1 -ve
            # NB next section applies whether 2nd test defined or not. However, it will overwrite in CXR case
            if J['dst'][ ptypes[i] ][0]!='DST_None':               #there is a DST for test 1+ve
                hb_log('algfromdict: testused' + str(i) + ' DST for test 1+ve')
                testused[i].setnext(1, testdict[ J['dst'][ptypes[i]][0] ] )   # DST if test 1 +-
                testused[i].setnext(2, testdict[ J['dst'][ptypes[i]][0] ] )   # DST if test 1 ++
        
        hb_log('algfromdict: After test linkage')
        # construct return value
        newalg = Algorithm( testused[0], testused[1], testused[2], testused[3] ) # the algorithm specified
        return newalg

###############################################################
#BEGIN JPMOD
###############################################################


    hb_log('jdata (in JSON): ' + json.dumps(jdata))

#LOAD THE HOMEBREW JPMOD
    Alg9 = algfromdict(jdata)                 # convert json data into algorithm; passed as python dict

    hb_log('Past the algfromdict() call')

    Algz = [Alg0,Alg1,Alg2,Alg3,Alg4,Alg5,Alg6,Alg7,Alg9]   # same as usual, but with last algorithm swapped out with new version
    AllFz = [Algz[i].getFuns() for i in range(9)] # get all the functions for later

    hb_log('Past All Homebrew Definitions')

###############################################################
#END JPMOD
###############################################################


##########################################################
# 4. DEFINE INITIAL POPULATION
##########################################################
# The values in each of the initial compartments were manually determined
#   to provide a reasonable starting population -- needed for the equation solver to
#   converge on appropriate roots later.

    I = np.zeros(100)           # initial population vector
    inc_fact = g.target_inc/100.  # factor based on inputted incidence
    lat_prob = 0.5*inc_fact/(1+0.5*inc_fact) # odds of latency relate to this

    mdr_prob = g.target_mdr
    inh_prob = g.target_mdr*2.3

    hiv_prob = g.target_hiv*0.01
    retx_prob = 0.2

# I[0] through I[99] are the initial "guesses" at each compartment size.
    I[1] = 68900*lat_prob*(1-inh_prob-mdr_prob)*(1-hiv_prob)
    I[2] = 100103*lat_prob*(inh_prob)*(1-hiv_prob)
    I[3] = 72300*lat_prob*(mdr_prob)*(1-hiv_prob)
    I[4:7] = I[1:4]*0.0026
    smpos = np.array([7,9,11])
    smneg = np.array([8,10,12])
    I[smpos]=I[1:4]*0.00046
    I[smneg]=I[1:4]*0.00089
    I[smpos+6] = I[1:4]*0.00003
    I[smneg+6] = I[1:4]*0.000015
    I[19]=I[1]*0.000035
    I[20]=I[1]*0.000019
    I[21]=I[2]*0.00017
    I[22]=I[2]*0.0001
    I[23]=I[3]*0.0006
    I[24]=I[3]*0.0003
    I[25] = 103000*(1-lat_prob+inh_prob+mdr_prob)*hiv_prob
    I[26:29] = I[1:4]* hiv_prob*0.97/((1-hiv_prob)*0.996)*0.7
    I[29:32] = I[26:29]*0.0053
    I[smpos+25] = I[26:29]*0.00576
    I[smneg+25] = I[26:29]*0.0135
    I[smpos+31] = I[26:29]*0.00035
    I[smneg+31] = I[26:29]*0.00025
    I[19+25]=I[1+25]*0.0003
    I[20+25]=I[1+25]*0.0002
    I[21+25]=I[2+25]*0.0015
    I[22+25]=I[2+25]*0.001
    I[23+25]=I[3+25]*0.0052
    I[24+25]=I[3+25]*0.0034
    I[50] = (1-lat_prob)*(g.prop_cough)*(1-hiv_prob)*100000.
    I[51:54] = I[1:4]* retx_prob/(1.-retx_prob)*0.7
    I[54:57] = I[51:54]*0.00064
    I[57] = I[51]*0.00018
    I[58] = I[51]*0.00034
    I[59] = I[52]*0.00042
    I[60] = I[52]*0.00077
    I[61] = I[53]*0.0011
    I[62] = I[53]*0.002
    I[63] = I[51]*0.000015
    I[64] = I[51]*0.000008
    I[65] = I[52]*0.00003
    I[66] = I[52]*0.000015
    I[67] = I[53]*0.00012
    I[68] = I[53]*0.000023
    I[19+50]=I[1+50]*0.000015
    I[20+50]=I[1+50]*0.000008
    I[21+50]=I[2+50]*0.0002
    I[22+50]=I[2+50]*0.0001
    I[23+50]=I[3+50]*0.0021
    I[24+50]=I[3+50]*0.0012
    I[75] = (1-lat_prob)*(g.prop_cough)*(hiv_prob)*100000.
    I[76:79] = I[26:29]* retx_prob/(1.-retx_prob)/0.7/0.7
    I[79:82] = I[76:79]*0.0043
    I[82] = I[76]*0.0048
    I[83] = I[76]*0.011
    I[84] = I[77]*0.0055
    I[85] = I[77]*0.012
    I[86] = I[78]*0.0075
    I[87] = I[78]*0.015
    I[88] = I[76]*0.00038
    I[89] = I[76]*0.00025
    I[90] = I[77]*0.00039
    I[91] = I[77]*0.00024
    I[92] = I[78]*0.0008
    I[93] = I[78]*0.00018
    I[94]=I[76]*0.00025
    I[95]=I[76]*0.00016
    I[96]=I[77]*0.0014
    I[97]=I[77]*0.00086
    I[98]=I[78]*0.0081
    I[99]=I[78]*0.0052

    I[0] = 100000.-np.sum(I[:])

    hb_log('Past Initial Population Definition')
##########################################################
# 5. DEFINE CALCULATION ROUTINES
##########################################################

# This function encodes the decision tree.
# It takes every combination of:
# intervention, smear/infectious status, HIV status, drug resistance status, and new/retreatment status,
#
# and then calculates:
# 0. Rate of successful treatment (dx rate * prob of successful tx)
# 1. Delay from sending diagnostic test to initiation of treatment
# 2. Rate of treatment that will fail (dx rate * prob of unsuccessful tx)
# 3. Cost of diagnosis
# 4. Cost of treatment if successful
# 5. Cost of treatment if unsuccessful

    def txrate(interv,smear,dr,retx,hiv):
        return AllFz[int(interv)][0](smear,dr,retx,hiv)

# calculate specificity of testing for people without TB, according to HIV status
    def invspec(interv, hiv, retx):
        return AllFz[int(interv)][1](hiv,retx)

# cost of diagnosing & treating TB-negatives who present with symptoms,
#    according to selected intervention, prior TB treatment status, and HIV status
    def dxcost(interv, retx, hiv):
        return AllFz[int(interv)][2](retx,hiv)


# simple routine to define proportion of rapid progression, according to HIV status
    def rap(hiv):
        if hiv<25 or (hiv>=50 and hiv<75):
            rp = rapid
        if (hiv>=25 and hiv<50) or (hiv>=75):
            rp = rapid_h
        return rp

# calculate mortality rate according to TB and HIV status
    def mortality(tb,hiv):
        m = g.mort
        if hiv==1:
            m += g.mort_h
        if tb==1:
            if hiv==0:
                m += g.mort_sp
            elif hiv==1:
                m += g.mort_tbh
        if tb==2:
            if hiv==0:
                m += g.mort_sn
            elif hiv==1:
                m += g.mort_tbh
        return m


##########################################################
    hb_log('Past Calculation Routines')

##########################################################
# 6. DIFFERENTIAL EQUATION FUNCTION
##########################################################

# This function defines the differential equations that are later solved.

    def diffeq(population,time):
        #Global values used locally
        tb_num = g.tb_num
        beta_s = g.beta_s
        intopt = g.intopt
        relbeta_sn = g.relbeta_sn
        fit_i = g.fit_i
        fit_m = g.fit_m
        rapid = g.rapid
        rapid_h = g.rapid_h
        react = g.react
        react_h = g.react_h
        prot = g.prot
        predx_dur = g.predx_dur
        predx_dur_h = g.predx_dur_h
        cure_sn = g.cure_sn
        prop_inf = g.prop_inf
        prop_inf_h = g.prop_inf_h
        acq_s = g.acq_s
        acq_i1 = g.acq_i1
        acq_i2 = g.acq_i2
        acq_mdr = g.acq_mdr
        dur_fail = g.dur_fail
        defvfail_s = g.defvfail_s
        defvfail_i = g.defvfail_i
        defvfail_m = g.defvfail_m
        fail_rt = g.fail_rt
        fail_i2 = g.fail_i2
        fail_m2 = g.fail_m2
        prop_cough = g.prop_cough
        cure_sp = g.cure_sp
     
        dxdt = np.zeros(tb_num)  # create the vector of ODEs
        X = population

        # The following equations add up all infectious compartments to calculate the force of infection
        active_s = np.array([7,13])              # smear-pos, DS-TB compartments
        active_s1 = active_s + 25
        active_s = np.concatenate((active_s,active_s1))
        active_s1 = active_s + 50
        active_s = np.concatenate((active_s,active_s1))

        smneg_s = np.array([4,8,14,19,20])       # smear-neg, DS-TB compartments
        smneg_s1 = smneg_s + 25
        smneg_s = np.concatenate((smneg_s,smneg_s1))
        smneg_s1 = smneg_s + 50
        smneg_s = np.concatenate((smneg_s,smneg_s1))

        active_i = np.array([9,15])              # smear-pos, INH-m compartments
        active_i1 = active_i + 25
        active_i = np.concatenate((active_i,active_i1))
        active_i1 = active_i + 50
        active_i = np.concatenate((active_i,active_i1))

        smneg_i = np.array([5,10,16,21,22])       # smear-neg, INH-m compartments
        smneg_i1 = smneg_i + 25
        smneg_i = np.concatenate((smneg_i,smneg_i1))
        smneg_i1 = smneg_i + 50
        smneg_i = np.concatenate((smneg_i,smneg_i1))

        active_m = np.array([11,17])              # smear-pos, MDR-TB compartments
        active_m1 = active_m + 25
        active_m = np.concatenate((active_m,active_m1))
        active_m1 = active_m + 50
        active_m = np.concatenate((active_m,active_m1))

        smneg_m = np.array([6,12,18,23,24])       # smear-neg, MDR-TB compartments
        smneg_m1 = smneg_m + 25
        smneg_m = np.concatenate((smneg_m,smneg_m1))
        smneg_m1 = smneg_m + 50
        smneg_m = np.concatenate((smneg_m,smneg_m1))

        # use the quantities above to calculate forces of infection:
        force_s = beta_s*(X[active_s].sum() + relbeta_sn*X[smneg_s].sum())/population.sum()
        force_i = beta_s*fit_i*(X[active_i].sum() + relbeta_sn*X[smneg_i].sum())/population.sum()
        force_m = beta_s*fit_m*(X[active_m].sum() + relbeta_sn*X[smneg_m].sum())/population.sum()
        force_t = force_s + force_i + force_m

        # create a 100x100 matrix of flow rates:
        # each element is the flow from row to column
        # thus, chg_arr[0,1] is the flow from compartment 0 (uninfected) to
        # compartment 1 (DS-latent)
        chg_arr = np.zeros((tb_num,tb_num))
        i1 = np.array([0,50]) # HIV-neg only
        i2 = np.array([25,75]) # HIV-pos only
        i3 = np.array([0,25,50,75]) # all
        i4 = np.array([0,25]) # new dx
        i5 = np.array([50,75]) # retx
        chg_arr[i1,i1+1] = force_s*(1-rapid)*X[i1] # infection that becomes latent
        chg_arr[i1,i1+2] = force_i*(1-rapid)*X[i1]
        chg_arr[i1,i1+3] = force_m*(1-rapid)*X[i1]
        chg_arr[i2,i2+1] = force_s*(1-rapid_h)*X[i2]
        chg_arr[i2,i2+2] = force_i*(1-rapid_h)*X[i2]
        chg_arr[i2,i2+3] = force_m*(1-rapid_h)*X[i2]
        chg_arr[i1,i1+4] = force_s*rapid*X[i1]     # primary progressive TB
        chg_arr[i1,i1+5] = force_i*rapid*X[i1]
        chg_arr[i1,i1+6] = force_m*rapid*X[i1]
        chg_arr[i2,i2+4] = force_s*rapid_h*X[i2]
        chg_arr[i2,i2+5] = force_i*rapid_h*X[i2]
        chg_arr[i2,i2+6] = force_m*rapid_h*X[i2]
        chg_arr[i1+1,i1+4] += react*X[i1+1]         # reactivation
        chg_arr[i1+2,i1+5] += react*X[i1+2]
        chg_arr[i1+3,i1+6] += react*X[i1+3]
        chg_arr[i2+1,i2+4] += react_h*X[i2+1]
        chg_arr[i2+2,i2+5] += react_h*X[i2+2]
        chg_arr[i2+3,i2+6] += react_h*X[i2+3]
        chg_arr[i1+1,i1+4] += force_s*rapid*(1-prot)*X[i1+1]     # reinfection to active
        chg_arr[i1+1,i1+5] += force_i*rapid*(1-prot)*X[i1+1]
        chg_arr[i1+1,i1+6] += force_m*rapid*(1-prot)*X[i1+1]
        chg_arr[i1+2,i1+4] += force_s*rapid*(1-prot)*X[i1+2]
        chg_arr[i1+2,i1+5] += force_i*rapid*(1-prot)*X[i1+2]
        chg_arr[i1+2,i1+6] += force_m*rapid*(1-prot)*X[i1+2]
        chg_arr[i1+3,i1+4] += force_s*rapid*(1-prot)*X[i1+3]
        chg_arr[i1+3,i1+5] += force_i*rapid*(1-prot)*X[i1+3]
        chg_arr[i1+3,i1+6] += force_m*rapid*(1-prot)*X[i1+3]
        chg_arr[i1+1,i1+2] += force_i*(1-rapid*(1-prot))*X[i1+1] # reinfection to latent
        chg_arr[i1+1,i1+3] += force_m*(1-rapid*(1-prot))*X[i1+1]
        chg_arr[i1+2,i1+1] += force_s*(1-rapid*(1-prot))*X[i1+2]
        chg_arr[i1+2,i1+3] += force_m*(1-rapid*(1-prot))*X[i1+2]
        chg_arr[i1+3,i1+1] += force_s*(1-rapid*(1-prot))*X[i1+3]
        chg_arr[i1+3,i1+2] += force_i*(1-rapid*(1-prot))*X[i1+3]
        chg_arr[i2+1,i2+4] += force_s*rapid_h*X[i2+1]     # reinfection to active (HIV)
        chg_arr[i2+1,i2+5] += force_i*rapid_h*X[i2+1]
        chg_arr[i2+1,i2+6] += force_m*rapid_h*X[i2+1]
        chg_arr[i2+2,i2+4] += force_s*rapid_h*X[i2+2]
        chg_arr[i2+2,i2+5] += force_i*rapid_h*X[i2+2]
        chg_arr[i2+2,i2+6] += force_m*rapid_h*X[i2+2]
        chg_arr[i2+3,i2+4] += force_s*rapid_h*X[i2+3]
        chg_arr[i2+3,i2+5] += force_i*rapid_h*X[i2+3]
        chg_arr[i2+3,i2+6] += force_m*rapid_h*X[i2+3]
        chg_arr[i1+4,i1+5] += force_i*(1-rapid_h)*X[i1+4]  # reinfection to latent (HIV)
        chg_arr[i1+4,i1+6] += force_m*(1-rapid_h)*X[i1+4]
        chg_arr[i1+5,i1+4] += force_s*(1-rapid_h)*X[i1+5]
        chg_arr[i1+5,i1+6] += force_m*(1-rapid_h)*X[i1+5]
        chg_arr[i1+6,i1+4] += force_s*(1-rapid_h)*X[i1+6]
        chg_arr[i1+6,i1+5] += force_i*(1-rapid_h)*X[i1+6]
        chg_arr[i1+4,i1+7] += (1/predx_dur-cure_sn)*prop_inf*X[i1+4]  # progression to dx-seeking
        chg_arr[i1+4,i1+8] += (1/predx_dur-cure_sn)*(1-prop_inf)*X[i1+4]
        chg_arr[i1+5,i1+9] += (1/predx_dur-cure_sn)*prop_inf*X[i1+5]
        chg_arr[i1+5,i1+10] += (1/predx_dur-cure_sn)*(1-prop_inf)*X[i1+5]
        chg_arr[i1+6,i1+11] += (1/predx_dur-cure_sn)*prop_inf*X[i1+6]
        chg_arr[i1+6,i1+12] += (1/predx_dur-cure_sn)*(1-prop_inf)*X[i1+6]
        chg_arr[i2+4,i2+7] += (1/predx_dur_h)*prop_inf_h*X[i2+4]  # progress to dx-seek, HIV
        chg_arr[i2+4,i2+8] += (1/predx_dur_h)*(1-prop_inf_h)*X[i2+4]
        chg_arr[i2+5,i2+9] += (1/predx_dur_h)*prop_inf_h*X[i2+5]
        chg_arr[i2+5,i2+10] += (1/predx_dur_h)*(1-prop_inf_h)*X[i2+5]
        chg_arr[i2+6,i2+11] += (1/predx_dur_h)*prop_inf_h*X[i2+6]
        chg_arr[i2+6,i2+12] += (1/predx_dur_h)*(1-prop_inf_h)*X[i2+6]

        # The equations below use the output from the decision-tree
        # "txrate" function.
        # Arguments for this function are: intervention, smear, DR, retreat, hiv
        # Outputs are: tx rate, dx delay, tx_inapp, cost_dx, cost_tx, cost_inapp
        chg_arr[7,13] += txrate(intopt,1,0,0,0)[0]*X[7] # progress to "dx in progress"
        chg_arr[8,14] += txrate(intopt,0,0,0,0)[0]*X[8] # new
        chg_arr[9,15] += txrate(intopt,1,1,0,0)[0]*X[9]
        chg_arr[10,16] += txrate(intopt,0,1,0,0)[0]*X[10]
        chg_arr[11,17] += txrate(intopt,1,2,0,0)[0]*X[11]
        chg_arr[12,18] += txrate(intopt,0,2,0,0)[0]*X[12]
        chg_arr[57,63] += txrate(intopt,1,0,1,0)[0]*X[57] # progress to "dx in progress"
        chg_arr[58,64] += txrate(intopt,0,0,1,0)[0]*X[58] # retreat
        chg_arr[59,65] += txrate(intopt,1,1,1,0)[0]*X[59]
        chg_arr[60,66] += txrate(intopt,0,1,1,0)[0]*X[60]
        chg_arr[61,67] += txrate(intopt,1,2,1,0)[0]*X[61]
        chg_arr[62,68] += txrate(intopt,0,2,1,0)[0]*X[62]
        chg_arr[7,19] += txrate(intopt,1,0,0,0)[2]*(1-acq_s)*X[7] # progress to "inapp tx"
        chg_arr[8,20] += txrate(intopt,0,0,0,0)[2]*(1-acq_s)*X[8] # new
        chg_arr[9,21] += txrate(intopt,1,1,0,0)[2]*(1-acq_i1)*X[9]
        chg_arr[10,22] += txrate(intopt,0,1,0,0)[2]*(1-acq_i1)*X[10]
        chg_arr[11,23] += txrate(intopt,1,2,0,0)[2]*X[11]
        chg_arr[12,24] += txrate(intopt,0,2,0,0)[2]*X[12]
        chg_arr[7,21] += txrate(intopt,1,0,0,0)[2]*(acq_s)*(1-acq_mdr)*X[7] # new resistance
        chg_arr[8,22] += txrate(intopt,0,0,0,0)[2]*(acq_s)*(1-acq_mdr)*X[8]
        chg_arr[7,23] += txrate(intopt,1,0,0,0)[2]*(acq_s)*(acq_mdr)*X[7] # new resistance
        chg_arr[8,24] += txrate(intopt,0,0,0,0)[2]*(acq_s)*(acq_mdr)*X[8]
        chg_arr[9,23] += txrate(intopt,1,1,0,0)[2]*(acq_i1)*X[9]
        chg_arr[10,24] += txrate(intopt,0,1,0,0)[2]*(acq_i1)*X[10]
        chg_arr[57,69] += txrate(intopt,1,0,1,0)[2]*(1-acq_s)*X[57] # progress to "inapp tx"
        chg_arr[58,70] += txrate(intopt,0,0,1,0)[2]*(1-acq_s)*X[58] # retreat
        chg_arr[59,71] += txrate(intopt,1,1,1,0)[2]*(1-acq_i2)*X[59]
        chg_arr[60,72] += txrate(intopt,0,1,1,0)[2]*(1-acq_i2)*X[60]
        chg_arr[61,73] += txrate(intopt,1,2,1,0)[2]*X[61]
        chg_arr[62,74] += txrate(intopt,0,2,1,0)[2]*X[62]
        chg_arr[57,71] += txrate(intopt,1,0,1,0)[2]*(acq_s)*(1-acq_mdr)*X[57] # new resistance
        chg_arr[58,72] += txrate(intopt,0,0,1,0)[2]*(acq_s)*(1-acq_mdr)*X[58]
        chg_arr[57,73] += txrate(intopt,1,0,1,0)[2]*(acq_s)*(acq_mdr)*X[57] # new resistance
        chg_arr[58,74] += txrate(intopt,0,0,1,0)[2]*(acq_s)*(acq_mdr)*X[58]
        chg_arr[59,73] += txrate(intopt,1,1,1,0)[2]*(acq_i2)*X[59]
        chg_arr[60,74] += txrate(intopt,0,1,1,0)[2]*(acq_i2)*X[60]
        chg_arr[13,51] += txrate(intopt,1,0,0,0)[1]*X[13] # progress to treated
        chg_arr[14,51] += txrate(intopt,0,0,0,0)[1]*X[14] # new
        chg_arr[15,52] += txrate(intopt,1,1,0,0)[1]*X[15]
        chg_arr[16,52] += txrate(intopt,0,1,0,0)[1]*X[16]
        chg_arr[17,53] += txrate(intopt,1,2,0,0)[1]*X[17]
        chg_arr[18,53] += txrate(intopt,0,2,0,0)[1]*X[18]
        chg_arr[63,51] += txrate(intopt,1,0,1,0)[1]*X[63] # progress to treated
        chg_arr[64,51] += txrate(intopt,0,0,1,0)[1]*X[64] # retreat
        chg_arr[65,52] += txrate(intopt,1,1,1,0)[1]*X[65]
        chg_arr[66,52] += txrate(intopt,0,1,1,0)[1]*X[66]
        chg_arr[67,53] += txrate(intopt,1,2,1,0)[1]*X[67]
        chg_arr[68,53] += txrate(intopt,0,2,1,0)[1]*X[68]

        chg_arr[7+25,13+25] += txrate(intopt,1,0,0,1)[0]*X[7+25] # progress to "dx in progress"
        chg_arr[8+25,14+25] += txrate(intopt,0,0,0,1)[0]*X[8+25] # new
        chg_arr[9+25,15+25] += txrate(intopt,1,1,0,1)[0]*X[9+25]
        chg_arr[10+25,16+25] += txrate(intopt,0,1,0,1)[0]*X[10+25]
        chg_arr[11+25,17+25] += txrate(intopt,1,2,0,1)[0]*X[11+25]
        chg_arr[12+25,18+25] += txrate(intopt,0,2,0,1)[0]*X[12+25]
        chg_arr[57+25,63+25] += txrate(intopt,1,0,1,1)[0]*X[57+25] # progress to "dx in progress"
        chg_arr[58+25,64+25] += txrate(intopt,0,0,1,1)[0]*X[58+25] # retreat
        chg_arr[59+25,65+25] += txrate(intopt,1,1,1,1)[0]*X[59+25]
        chg_arr[60+25,66+25] += txrate(intopt,0,1,1,1)[0]*X[60+25]
        chg_arr[61+25,67+25] += txrate(intopt,1,2,1,1)[0]*X[61+25]
        chg_arr[62+25,68+25] += txrate(intopt,0,2,1,1)[0]*X[62+25]
        chg_arr[7+25,19+25] += txrate(intopt,1,0,0,1)[2]*(1-acq_s)*X[7+25] # progress to "inapp tx"
        chg_arr[8+25,20+25] += txrate(intopt,0,0,0,1)[2]*(1-acq_s)*X[8+25] # new
        chg_arr[9+25,21+25] += txrate(intopt,1,1,0,1)[2]*(1-acq_i1)*X[9+25]
        chg_arr[10+25,22+25] += txrate(intopt,0,1,0,1)[2]*(1-acq_i1)*X[10+25]
        chg_arr[11+25,23+25] += txrate(intopt,1,2,0,1)[2]*X[11+25]
        chg_arr[12+25,24+25] += txrate(intopt,0,2,0,1)[2]*X[12+25]
        chg_arr[7+25,21+25] += txrate(intopt,1,0,0,1)[2]*(acq_s)*(1-acq_mdr)*X[7+25] # new resistance
        chg_arr[8+25,22+25] += txrate(intopt,0,0,0,1)[2]*(acq_s)*(1-acq_mdr)*X[8+25]
        chg_arr[7+25,23+25] += txrate(intopt,1,0,0,1)[2]*(acq_s)*(acq_mdr)*X[7+25] # new resistance
        chg_arr[8+25,24+25] += txrate(intopt,0,0,0,1)[2]*(acq_s)*(acq_mdr)*X[8+25]
        chg_arr[9+25,23+25] += txrate(intopt,1,1,0,1)[2]*(acq_i1)*X[9+25]
        chg_arr[10+25,24+25] += txrate(intopt,0,1,0,1)[2]*(acq_i1)*X[10+25]
        chg_arr[57+25,69+25] += txrate(intopt,1,0,1,1)[2]*(1-acq_s)*X[57+25] # progress to "inapp tx"
        chg_arr[58+25,70+25] += txrate(intopt,0,0,1,1)[2]*(1-acq_s)*X[58+25] # retreat
        chg_arr[59+25,71+25] += txrate(intopt,1,1,1,1)[2]*(1-acq_i2)*X[59+25]
        chg_arr[60+25,72+25] += txrate(intopt,0,1,1,1)[2]*(1-acq_i2)*X[60+25]
        chg_arr[61+25,73+25] += txrate(intopt,1,2,1,1)[2]*X[61+25]
        chg_arr[62+25,74+25] += txrate(intopt,0,2,1,1)[2]*X[62+25]
        chg_arr[57+25,71+25] += txrate(intopt,1,0,1,1)[2]*(acq_s)*(1-acq_mdr)*X[57+25] # new resistance
        chg_arr[58+25,72+25] += txrate(intopt,0,0,1,1)[2]*(acq_s)*(1-acq_mdr)*X[58+25]
        chg_arr[57+25,73+25] += txrate(intopt,1,0,1,1)[2]*(acq_s)*(acq_mdr)*X[57+25] # new resistance
        chg_arr[58+25,74+25] += txrate(intopt,0,0,1,1)[2]*(acq_s)*(acq_mdr)*X[58+25]
        chg_arr[59+25,73+25] += txrate(intopt,1,1,1,1)[2]*(acq_i2)*X[59+25]
        chg_arr[60+25,74+25] += txrate(intopt,0,1,1,1)[2]*(acq_i2)*X[60+25]
        chg_arr[13+25,51+25] += txrate(intopt,1,0,0,1)[1]*X[13+25] # progress to treated
        chg_arr[14+25,51+25] += txrate(intopt,0,0,0,1)[1]*X[14+25] # new
        chg_arr[15+25,52+25] += txrate(intopt,1,1,0,1)[1]*X[15+25]
        chg_arr[16+25,52+25] += txrate(intopt,0,1,0,1)[1]*X[16+25]
        chg_arr[17+25,53+25] += txrate(intopt,1,2,0,1)[1]*X[17+25]
        chg_arr[18+25,53+25] += txrate(intopt,0,2,0,1)[1]*X[18+25]
        chg_arr[63+25,51+25] += txrate(intopt,1,0,1,1)[1]*X[63+25] # progress to treated
        chg_arr[64+25,51+25] += txrate(intopt,0,0,1,1)[1]*X[64+25] # retreat
        chg_arr[65+25,52+25] += txrate(intopt,1,1,1,1)[1]*X[65+25]
        chg_arr[66+25,52+25] += txrate(intopt,0,1,1,1)[1]*X[66+25]
        chg_arr[67+25,53+25] += txrate(intopt,1,2,1,1)[1]*X[67+25]
        chg_arr[68+25,53+25] += txrate(intopt,0,2,1,1)[1]*X[68+25]

        chg_arr[i4+19,i4+57] += (12./dur_fail)*defvfail_s*X[i4+19] # default from "inapp tx" (new cases)
        chg_arr[i4+20,i4+58] += (12./dur_fail)*defvfail_s*X[i4+20] # assume zero cost b/c 1/2 of default is cured
        chg_arr[i4+21,i4+59] += (12./dur_fail)*defvfail_i*X[i4+21]
        chg_arr[i4+22,i4+60] += (12./dur_fail)*defvfail_i*X[i4+22]
        chg_arr[i4+23,i4+61] += (12./dur_fail)*defvfail_m*X[i4+23]
        chg_arr[i4+24,i4+62] += (12./dur_fail)*defvfail_m*X[i4+24]
        chg_arr[i5+19,i5+7] += (12./dur_fail)*defvfail_s*X[i5+19] # default from "inapp tx" (retx cases)
        chg_arr[i5+20,i5+8] += (12./dur_fail)*defvfail_s*X[i5+20] # assume zero cost b/c 1/2 of default is cured
        chg_arr[i5+21,i5+9] += (12./dur_fail)*defvfail_i*X[i5+21]
        chg_arr[i5+22,i5+10] += (12./dur_fail)*defvfail_i*X[i5+22]
        chg_arr[i5+23,i5+11] += (12./dur_fail)*defvfail_m*X[i5+23]
        chg_arr[i5+24,i5+12] += (12./dur_fail)*defvfail_m*X[i5+24]
        chg_arr[i4+19,i4+51] += (12./dur_fail)*(1-defvfail_s)*X[i4+19]*(1-fail_rt) # failure from "inapp tx" (new TB cases)
        chg_arr[i4+20,i4+51] += (12./dur_fail)*(1-defvfail_s)*X[i4+20]*(1-fail_rt) # these are the pts that are then successfully treated
        chg_arr[i4+21,i4+52] += (12./dur_fail)*(1-defvfail_i)*X[i4+21]*(1-fail_i2) # the ones below re-fail
        chg_arr[i4+22,i4+52] += (12./dur_fail)*(1-defvfail_i)*X[i4+22]*(1-fail_i2)
        chg_arr[i4+23,i4+53] += (12./dur_fail)*(1-defvfail_m)*X[i4+23]*(1-fail_m2)
        chg_arr[i4+24,i4+53] += (12./dur_fail)*(1-defvfail_m)*X[i4+24]*(1-fail_m2)
        chg_arr[i4+19,i4+69] += (12./dur_fail)*(1-defvfail_s)*X[i4+19]*(fail_rt) # re-failure from "inapp tx" -> now classified as retx
        chg_arr[i4+20,i4+70] += (12./dur_fail)*(1-defvfail_s)*X[i4+20]*(fail_rt)
        chg_arr[i4+21,i4+71] += (12./dur_fail)*(1-defvfail_i)*X[i4+21]*(fail_i2)
        chg_arr[i4+22,i4+72] += (12./dur_fail)*(1-defvfail_i)*X[i4+22]*(fail_i2)
        chg_arr[i4+23,i4+73] += (12./dur_fail)*(1-defvfail_m)*X[i4+23]*(fail_m2)
        chg_arr[i4+24,i4+74] += (12./dur_fail)*(1-defvfail_m)*X[i4+24]*(fail_m2)
        chg_arr[i5+19,i5+1] += (12./dur_fail)*(1-defvfail_s)*X[i5+19]*(1-fail_rt) # failure from "inapp tx"
        chg_arr[i5+20,i5+1] += (12./dur_fail)*(1-defvfail_s)*X[i5+20]*(1-fail_rt) # retx pts, so "re-fail" stays in the same box
        chg_arr[i5+21,i5+2] += (12./dur_fail)*(1-defvfail_i)*X[i5+21]*(1-fail_i2)
        chg_arr[i5+22,i5+2] += (12./dur_fail)*(1-defvfail_i)*X[i5+22]*(1-fail_i2)
        chg_arr[i5+23,i5+3] += (12./dur_fail)*(1-defvfail_m)*X[i5+23]*(1-fail_m2)
        chg_arr[i5+24,i5+3] += (12./dur_fail)*(1-defvfail_m)*X[i5+24]*(1-fail_m2)
        i6 = np.arange(25)                          # HIV infection
        chg_arr[i6,i6+25] += g.hiv_inc*X[i6]
        chg_arr[i6+50,i6+75] += g.hiv_inc*X[i6+50]
        i7 = np.arange(4)                           # inappropriate tx for TB
        chg_arr[i7,i7+50] += prop_cough*(1-invspec(intopt,0,0))*X[i7]
        chg_arr[i7+25,i7+75] += prop_cough*(1-invspec(intopt,1,0))*X[i7+25]
        i8 = np.arange(50)                         # mortality
        chg_arr[i6,0] += mortality(0,0)*X[i6]       # background mortality
        chg_arr[i6+50,0] += mortality(0,0)*X[i6+50]
        chg_arr[i6+25,0] += mortality(0,1)*X[i6+25] # HIV deaths
        chg_arr[i6+75,0] += mortality(0,1)*X[i6+75]
        i9 = np.array([7,9,11,13,15,17])  # smear-pos
        i10 = np.array([4,5,6,8,10,12,14,16,18,19,20,21,22,23,24])  # smear-neg
        chg_arr[i9,0] += mortality(1,0)*X[i9]              # smear-positive TB deaths
        chg_arr[i9+25,0] += mortality(1,1)*X[i9+25]
        chg_arr[i9+50,0] += mortality(1,0)*X[i9+50]
        chg_arr[i9+75,0] += mortality(1,1)*X[i9+75]
        chg_arr[i10,0] += mortality(2,0)*X[i10]              # smear-negative TB deaths
        chg_arr[i10+25,0] += mortality(2,1)*X[i10+25]
        chg_arr[i10+50,0] += mortality(2,0)*X[i10+50]
        chg_arr[i10+75,0] += mortality(2,1)*X[i10+75]
        chg_arr[i1+7,i1+1] += cure_sp*X[i1+7]              # smear-positive self-cure (no HIV)
        chg_arr[i1+9,i1+2] += cure_sp*X[i1+9]
        chg_arr[i1+11,i1+3] += cure_sp*X[i1+11]
        chg_arr[i1+13,i1+1] += cure_sp*X[i1+13]
        chg_arr[i1+15,i1+2] += cure_sp*X[i1+15]
        chg_arr[i1+17,i1+3] += cure_sp*X[i1+17]
        chg_arr[i1+8,i1+1] += cure_sn*X[i1+8]              # smear-negative self-cure (no HIV)
        chg_arr[i1+10,i1+2] += cure_sn*X[i1+10]
        chg_arr[i1+12,i1+3] += cure_sn*X[i1+12]
        chg_arr[i1+14,i1+1] += cure_sn*X[i1+14]
        chg_arr[i1+16,i1+2] += cure_sn*X[i1+16]
        chg_arr[i1+18,i1+3] += cure_sn*X[i1+18]
        chg_arr[i1+4,i1+1] += cure_sn*X[i1+4]
        chg_arr[i1+5,i1+2] += cure_sn*X[i1+5]
        chg_arr[i1+6,i1+3] += cure_sn*X[i1+6]
        chg_arr[i1+19,i1+1] += cure_sn*X[i1+19]
        chg_arr[i1+21,i1+2] += cure_sn*X[i1+21]
        chg_arr[i1+23,i1+3] += cure_sn*X[i1+23]
        chg_arr[i1+20,i1+1] += cure_sn*X[i1+20]
        chg_arr[i1+22,i1+2] += cure_sn*X[i1+22]
        chg_arr[i1+24,i1+3] += cure_sn*X[i1+24]

        # differential equations that are the output of the function:
        # each equation is "sum of all inputs minus sum of all outputs"
        dxdt[0] = np.sum(chg_arr[:,0])-np.sum(chg_arr[0,:])
        dxdt[1] = np.sum(chg_arr[:,1])-np.sum(chg_arr[1,:])
        dxdt[2] = np.sum(chg_arr[:,2])-np.sum(chg_arr[2,:])
        dxdt[3] = np.sum(chg_arr[:,3])-np.sum(chg_arr[3,:])
        dxdt[4] = np.sum(chg_arr[:,4])-np.sum(chg_arr[4,:])
        dxdt[5] = np.sum(chg_arr[:,5])-np.sum(chg_arr[5,:])
        dxdt[6] = np.sum(chg_arr[:,6])-np.sum(chg_arr[6,:])
        dxdt[7] = np.sum(chg_arr[:,7])-np.sum(chg_arr[7,:])
        dxdt[8] = np.sum(chg_arr[:,8])-np.sum(chg_arr[8,:])
        dxdt[9] = np.sum(chg_arr[:,9])-np.sum(chg_arr[9,:])
        dxdt[10] = np.sum(chg_arr[:,10])-np.sum(chg_arr[10,:])
        dxdt[11] = np.sum(chg_arr[:,11])-np.sum(chg_arr[11,:])
        dxdt[12] = np.sum(chg_arr[:,12])-np.sum(chg_arr[12,:])
        dxdt[13] = np.sum(chg_arr[:,13])-np.sum(chg_arr[13,:])
        dxdt[14] = np.sum(chg_arr[:,14])-np.sum(chg_arr[14,:])
        dxdt[15] = np.sum(chg_arr[:,15])-np.sum(chg_arr[15,:])
        dxdt[16] = np.sum(chg_arr[:,16])-np.sum(chg_arr[16,:])
        dxdt[17] = np.sum(chg_arr[:,17])-np.sum(chg_arr[17,:])
        dxdt[18] = np.sum(chg_arr[:,18])-np.sum(chg_arr[18,:])
        dxdt[19] = np.sum(chg_arr[:,19])-np.sum(chg_arr[19,:])
        dxdt[20] = np.sum(chg_arr[:,20])-np.sum(chg_arr[20,:])
        dxdt[21] = np.sum(chg_arr[:,21])-np.sum(chg_arr[21,:])
        dxdt[22] = np.sum(chg_arr[:,22])-np.sum(chg_arr[22,:])
        dxdt[23] = np.sum(chg_arr[:,23])-np.sum(chg_arr[23,:])
        dxdt[24] = np.sum(chg_arr[:,24])-np.sum(chg_arr[24,:])
        dxdt[25] = np.sum(chg_arr[:,25])-np.sum(chg_arr[25,:])
        dxdt[26] = np.sum(chg_arr[:,26])-np.sum(chg_arr[26,:])
        dxdt[27] = np.sum(chg_arr[:,27])-np.sum(chg_arr[27,:])
        dxdt[28] = np.sum(chg_arr[:,28])-np.sum(chg_arr[28,:])
        dxdt[29] = np.sum(chg_arr[:,29])-np.sum(chg_arr[29,:])
        dxdt[30] = np.sum(chg_arr[:,30])-np.sum(chg_arr[30,:])
        dxdt[31] = np.sum(chg_arr[:,31])-np.sum(chg_arr[31,:])
        dxdt[32] = np.sum(chg_arr[:,32])-np.sum(chg_arr[32,:])
        dxdt[33] = np.sum(chg_arr[:,33])-np.sum(chg_arr[33,:])
        dxdt[34] = np.sum(chg_arr[:,34])-np.sum(chg_arr[34,:])
        dxdt[35] = np.sum(chg_arr[:,35])-np.sum(chg_arr[35,:])
        dxdt[36] = np.sum(chg_arr[:,36])-np.sum(chg_arr[36,:])
        dxdt[37] = np.sum(chg_arr[:,37])-np.sum(chg_arr[37,:])
        dxdt[38] = np.sum(chg_arr[:,38])-np.sum(chg_arr[38,:])
        dxdt[39] = np.sum(chg_arr[:,39])-np.sum(chg_arr[39,:])
        dxdt[40] = np.sum(chg_arr[:,40])-np.sum(chg_arr[40,:])
        dxdt[41] = np.sum(chg_arr[:,41])-np.sum(chg_arr[41,:])
        dxdt[42] = np.sum(chg_arr[:,42])-np.sum(chg_arr[42,:])
        dxdt[43] = np.sum(chg_arr[:,43])-np.sum(chg_arr[43,:])
        dxdt[44] = np.sum(chg_arr[:,44])-np.sum(chg_arr[44,:])
        dxdt[45] = np.sum(chg_arr[:,45])-np.sum(chg_arr[45,:])
        dxdt[46] = np.sum(chg_arr[:,46])-np.sum(chg_arr[46,:])
        dxdt[47] = np.sum(chg_arr[:,47])-np.sum(chg_arr[47,:])
        dxdt[48] = np.sum(chg_arr[:,48])-np.sum(chg_arr[48,:])
        dxdt[49] = np.sum(chg_arr[:,49])-np.sum(chg_arr[49,:])
        dxdt[50] = np.sum(chg_arr[:,50])-np.sum(chg_arr[50,:])
        dxdt[51] = np.sum(chg_arr[:,51])-np.sum(chg_arr[51,:])
        dxdt[52] = np.sum(chg_arr[:,52])-np.sum(chg_arr[52,:])
        dxdt[53] = np.sum(chg_arr[:,53])-np.sum(chg_arr[53,:])
        dxdt[54] = np.sum(chg_arr[:,54])-np.sum(chg_arr[54,:])
        dxdt[55] = np.sum(chg_arr[:,55])-np.sum(chg_arr[55,:])
        dxdt[56] = np.sum(chg_arr[:,56])-np.sum(chg_arr[56,:])
        dxdt[57] = np.sum(chg_arr[:,57])-np.sum(chg_arr[57,:])
        dxdt[58] = np.sum(chg_arr[:,58])-np.sum(chg_arr[58,:])
        dxdt[59] = np.sum(chg_arr[:,59])-np.sum(chg_arr[59,:])
        dxdt[60] = np.sum(chg_arr[:,60])-np.sum(chg_arr[60,:])
        dxdt[61] = np.sum(chg_arr[:,61])-np.sum(chg_arr[61,:])
        dxdt[62] = np.sum(chg_arr[:,62])-np.sum(chg_arr[62,:])
        dxdt[63] = np.sum(chg_arr[:,63])-np.sum(chg_arr[63,:])
        dxdt[64] = np.sum(chg_arr[:,64])-np.sum(chg_arr[64,:])
        dxdt[65] = np.sum(chg_arr[:,65])-np.sum(chg_arr[65,:])
        dxdt[66] = np.sum(chg_arr[:,66])-np.sum(chg_arr[66,:])
        dxdt[67] = np.sum(chg_arr[:,67])-np.sum(chg_arr[67,:])
        dxdt[68] = np.sum(chg_arr[:,68])-np.sum(chg_arr[68,:])
        dxdt[69] = np.sum(chg_arr[:,69])-np.sum(chg_arr[69,:])
        dxdt[70] = np.sum(chg_arr[:,70])-np.sum(chg_arr[70,:])
        dxdt[71] = np.sum(chg_arr[:,71])-np.sum(chg_arr[71,:])
        dxdt[72] = np.sum(chg_arr[:,72])-np.sum(chg_arr[72,:])
        dxdt[73] = np.sum(chg_arr[:,73])-np.sum(chg_arr[73,:])
        dxdt[74] = np.sum(chg_arr[:,74])-np.sum(chg_arr[74,:])
        dxdt[75] = np.sum(chg_arr[:,75])-np.sum(chg_arr[75,:])
        dxdt[76] = np.sum(chg_arr[:,76])-np.sum(chg_arr[76,:])
        dxdt[77] = np.sum(chg_arr[:,77])-np.sum(chg_arr[77,:])
        dxdt[78] = np.sum(chg_arr[:,78])-np.sum(chg_arr[78,:])
        dxdt[79] = np.sum(chg_arr[:,79])-np.sum(chg_arr[79,:])
        dxdt[80] = np.sum(chg_arr[:,80])-np.sum(chg_arr[80,:])
        dxdt[81] = np.sum(chg_arr[:,81])-np.sum(chg_arr[81,:])
        dxdt[82] = np.sum(chg_arr[:,82])-np.sum(chg_arr[82,:])
        dxdt[83] = np.sum(chg_arr[:,83])-np.sum(chg_arr[83,:])
        dxdt[84] = np.sum(chg_arr[:,84])-np.sum(chg_arr[84,:])
        dxdt[85] = np.sum(chg_arr[:,85])-np.sum(chg_arr[85,:])
        dxdt[86] = np.sum(chg_arr[:,86])-np.sum(chg_arr[86,:])
        dxdt[87] = np.sum(chg_arr[:,87])-np.sum(chg_arr[87,:])
        dxdt[88] = np.sum(chg_arr[:,88])-np.sum(chg_arr[88,:])
        dxdt[89] = np.sum(chg_arr[:,89])-np.sum(chg_arr[89,:])
        dxdt[90] = np.sum(chg_arr[:,90])-np.sum(chg_arr[90,:])
        dxdt[91] = np.sum(chg_arr[:,91])-np.sum(chg_arr[91,:])
        dxdt[92] = np.sum(chg_arr[:,92])-np.sum(chg_arr[92,:])
        dxdt[93] = np.sum(chg_arr[:,93])-np.sum(chg_arr[93,:])
        dxdt[94] = np.sum(chg_arr[:,94])-np.sum(chg_arr[94,:])
        dxdt[95] = np.sum(chg_arr[:,95])-np.sum(chg_arr[95,:])
        dxdt[96] = np.sum(chg_arr[:,96])-np.sum(chg_arr[96,:])
        dxdt[97] = np.sum(chg_arr[:,97])-np.sum(chg_arr[97,:])
        dxdt[98] = np.sum(chg_arr[:,98])-np.sum(chg_arr[98,:])
        dxdt[99] = np.sum(chg_arr[:,99])-np.sum(chg_arr[99,:])
        return dxdt
###########################################################

    hb_log('Past Differential Equation Calculator')

##########################################################
# 7. SOLVE FOR EQUILIBRIUM
##########################################################

# The following equation replicates the function in
#      6 above (differential equation function).
# However, it substitutes 3 quantities: rec25, rec75, and rec50.
# These quantities are used as placeholders for the given populations,
#      X[25], X[50], and X[75].
# Thus, rec25 denotes the population of the 25th compartment,
#      rec50 denotes the population of the 50th compartment,
#      rec75 denotes the population of the 75th compartment.
# The "initial guess" for these populations is that they are equal in size.
# Concomitantly, X[25] denotes hiv incidence,
#       X[50] denotes the transmission parameter beta,
#       X[75] denotes the reduction in transmission rate for MDR-TB.
# The rates of change in X[25], X[50], and X[75] are now set to
#       the difference between calculated and target HIV prevalence,
#       between calculated and target TB incidence, and
#       between calculated and target MDR-TB incidence, respectively.
# Thus, at the roots of these equations (dx/dt = 0), the population
#      has achieved the target HIV prevalence, TB incidence, and MDR-TB incidence.
# Equations 100 and 101 are then used to maintain a constant population,
#      ensuring that rec25 and rec50 are also not changing at steady-state.

    def solvebeta(population):
        dxdt = np.zeros(g.tb_num+2)  # create the vector of ODEs
        X = population
        rec1 = 100000 - np.sum(X) + X[25] + X[50] + X[75] + X[100] + X[101]
        rec25 = (rec1+X[100]+X[101])/3.
        rec75 = rec25-X[100] # X[100] is the difference b/w X[25] & X[75]
        rec50 = rec25-X[101] # X[101] is the difference b/w X[50] & X[25]

        active_s = np.array([7,13])              # smear-pos, DS-TB compartments
        active_s1 = active_s + 25
        active_s = np.concatenate((active_s,active_s1))
        active_s1 = active_s + 50
        active_s = np.concatenate((active_s,active_s1))

        smneg_s = np.array([4,8,14,19,20])       # smear-neg, DS-TB compartments
        smneg_s1 = smneg_s + 25
        smneg_s = np.concatenate((smneg_s,smneg_s1))
        smneg_s1 = smneg_s + 50
        smneg_s = np.concatenate((smneg_s,smneg_s1))

        active_i = np.array([9,15])              # smear-pos, INH-m compartments
        active_i1 = active_i + 25
        active_i = np.concatenate((active_i,active_i1))
        active_i1 = active_i + 50
        active_i = np.concatenate((active_i,active_i1))

        smneg_i = np.array([5,10,16,21,22])       # smear-neg, INH-m compartments
        smneg_i1 = smneg_i + 25
        smneg_i = np.concatenate((smneg_i,smneg_i1))
        smneg_i1 = smneg_i + 50
        smneg_i = np.concatenate((smneg_i,smneg_i1))

        active_m = np.array([11,17])              # smear-pos, MDR-TB compartments
        active_m1 = active_m + 25
        active_m = np.concatenate((active_m,active_m1))
        active_m1 = active_m + 50
        active_m = np.concatenate((active_m,active_m1))

        smneg_m = np.array([6,12,18,23,24])       # smear-neg, MDR-TB compartments
        smneg_m1 = smneg_m + 25
        smneg_m = np.concatenate((smneg_m,smneg_m1))
        smneg_m1 = smneg_m + 50
        smneg_m = np.concatenate((smneg_m,smneg_m1))

        # forces of infection:
        force_s = X[50]*(X[active_s].sum() + g.relbeta_sn*X[smneg_s].sum())/100000
        force_i = X[50]*(1.-X[75]/4.)*(X[active_i].sum() + g.relbeta_sn*X[smneg_i].sum())/100000
        force_m = X[50]*(1.-X[75])*(X[active_m].sum() + g.relbeta_sn*X[smneg_m].sum())/100000
        force_t = force_s + force_i + force_m

        # 100x100 matrix of flow rates:
        chg_arr = np.zeros((g.tb_num,g.tb_num))
        i1 = np.array([0,50]) # HIV-neg only
        i2 = np.array([25,75]) # HIV-pos only
        i3 = np.array([0,25,50,75]) # all
        i4 = np.array([0,25]) # new dx
        i5 = np.array([50,75]) # retx
        chg_arr[0,1] = force_s*(1-g.rapid)*X[0] # infection that becomes latent
        chg_arr[0,2] = force_i*(1-g.rapid)*X[0]
        chg_arr[0,3] = force_m*(1-g.rapid)*X[0]
        chg_arr[50,51] = force_s*(1-g.rapid)*rec50 # infection that becomes latent
        chg_arr[50,52] = force_i*(1-g.rapid)*rec50
        chg_arr[50,53] = force_m*(1-g.rapid)*rec50
        chg_arr[25,26] = force_s*(1-g.rapid_h)*rec25
        chg_arr[25,27] = force_i*(1-g.rapid_h)*rec25
        chg_arr[25,28] = force_m*(1-g.rapid_h)*rec25
        chg_arr[75,76] = force_s*(1-g.rapid_h)*rec75
        chg_arr[75,77] = force_i*(1-g.rapid_h)*rec75
        chg_arr[75,78] = force_m*(1-g.rapid_h)*rec75
        chg_arr[0,4] = force_s*g.rapid*X[0]     # primary progressive TB
        chg_arr[0,5] = force_i*g.rapid*X[0]
        chg_arr[0,6] = force_m*g.rapid*X[0]
        chg_arr[50,54] = force_s*g.rapid*rec50     # primary progressive TB
        chg_arr[50,55] = force_i*g.rapid*rec50
        chg_arr[50,56] = force_m*g.rapid*rec50
        chg_arr[25,29] = force_s*g.rapid_h*rec25
        chg_arr[25,30] = force_i*g.rapid_h*rec25
        chg_arr[25,31] = force_m*g.rapid_h*rec25
        chg_arr[75,79] = force_s*g.rapid_h*rec75
        chg_arr[75,80] = force_i*g.rapid_h*rec75
        chg_arr[75,81] = force_m*g.rapid_h*rec75
        chg_arr[i1+1,i1+4] += g.react*X[i1+1]         # reactivation
        chg_arr[i1+2,i1+5] += g.react*X[i1+2]
        chg_arr[i1+3,i1+6] += g.react*X[i1+3]
        chg_arr[i2+1,i2+4] += g.react_h*X[i2+1]
        chg_arr[i2+2,i2+5] += g.react_h*X[i2+2]
        chg_arr[i2+3,i2+6] += g.react_h*X[i2+3]
        chg_arr[i1+1,i1+4] += force_s*g.rapid*(1-g.prot)*X[i1+1]     # reinfection to active
        chg_arr[i1+1,i1+5] += force_i*g.rapid*(1-g.prot)*X[i1+1]
        chg_arr[i1+1,i1+6] += force_m*g.rapid*(1-g.prot)*X[i1+1]
        chg_arr[i1+2,i1+4] += force_s*g.rapid*(1-g.prot)*X[i1+2]
        chg_arr[i1+2,i1+5] += force_i*g.rapid*(1-g.prot)*X[i1+2]
        chg_arr[i1+2,i1+6] += force_m*g.rapid*(1-g.prot)*X[i1+2]
        chg_arr[i1+3,i1+4] += force_s*g.rapid*(1-g.prot)*X[i1+3]
        chg_arr[i1+3,i1+5] += force_i*g.rapid*(1-g.prot)*X[i1+3]
        chg_arr[i1+3,i1+6] += force_m*g.rapid*(1-g.prot)*X[i1+3]
        chg_arr[i1+1,i1+2] += force_i*(1-g.rapid*(1-g.prot))*X[i1+1] # reinfection to latent
        chg_arr[i1+1,i1+3] += force_m*(1-g.rapid*(1-g.prot))*X[i1+1]
        chg_arr[i1+2,i1+1] += force_s*(1-g.rapid*(1-g.prot))*X[i1+2]
        chg_arr[i1+2,i1+3] += force_m*(1-g.rapid*(1-g.prot))*X[i1+2]
        chg_arr[i1+3,i1+1] += force_s*(1-g.rapid*(1-g.prot))*X[i1+3]
        chg_arr[i1+3,i1+2] += force_i*(1-g.rapid*(1-g.prot))*X[i1+3]
        chg_arr[i2+1,i2+4] += force_s*g.rapid_h*X[i2+1]     # reinfection to active (HIV)
        chg_arr[i2+1,i2+5] += force_i*g.rapid_h*X[i2+1]
        chg_arr[i2+1,i2+6] += force_m*g.rapid_h*X[i2+1]
        chg_arr[i2+2,i2+4] += force_s*g.rapid_h*X[i2+2]
        chg_arr[i2+2,i2+5] += force_i*g.rapid_h*X[i2+2]
        chg_arr[i2+2,i2+6] += force_m*g.rapid_h*X[i2+2]
        chg_arr[i2+3,i2+4] += force_s*g.rapid_h*X[i2+3]
        chg_arr[i2+3,i2+5] += force_i*g.rapid_h*X[i2+3]
        chg_arr[i2+3,i2+6] += force_m*g.rapid_h*X[i2+3]
        chg_arr[i1+4,i1+5] += force_i*(1-g.rapid_h)*X[i1+4]  # reinfection to latent (HIV)
        chg_arr[i1+4,i1+6] += force_m*(1-g.rapid_h)*X[i1+4]
        chg_arr[i1+5,i1+4] += force_s*(1-g.rapid_h)*X[i1+5]
        chg_arr[i1+5,i1+6] += force_m*(1-g.rapid_h)*X[i1+5]
        chg_arr[i1+6,i1+4] += force_s*(1-g.rapid_h)*X[i1+6]
        chg_arr[i1+6,i1+5] += force_i*(1-g.rapid_h)*X[i1+6]
        chg_arr[i1+4,i1+7] += (1/g.predx_dur-g.cure_sn)*g.prop_inf*X[i1+4]  # progression to dx-seeking
        chg_arr[i1+4,i1+8] += (1/g.predx_dur-g.cure_sn)*(1-g.prop_inf)*X[i1+4]
        chg_arr[i1+5,i1+9] += (1/g.predx_dur-g.cure_sn)*g.prop_inf*X[i1+5]
        chg_arr[i1+5,i1+10] += (1/g.predx_dur-g.cure_sn)*(1-g.prop_inf)*X[i1+5]
        chg_arr[i1+6,i1+11] += (1/g.predx_dur-g.cure_sn)*g.prop_inf*X[i1+6]
        chg_arr[i1+6,i1+12] += (1/g.predx_dur-g.cure_sn)*(1-g.prop_inf)*X[i1+6]
        chg_arr[i2+4,i2+7] += (1/g.predx_dur_h)*g.prop_inf_h*X[i2+4]  # progress to dx-seek, HIV
        chg_arr[i2+4,i2+8] += (1/g.predx_dur_h)*(1-g.prop_inf_h)*X[i2+4]
        chg_arr[i2+5,i2+9] += (1/g.predx_dur_h)*g.prop_inf_h*X[i2+5]
        chg_arr[i2+5,i2+10] += (1/g.predx_dur_h)*(1-g.prop_inf_h)*X[i2+5]
        chg_arr[i2+6,i2+11] += (1/g.predx_dur_h)*g.prop_inf_h*X[i2+6]
        chg_arr[i2+6,i2+12] += (1/g.predx_dur_h)*(1-g.prop_inf_h)*X[i2+6]

        chg_arr[7,13] += txrate(g.intopt,1,0,0,0)[0]*X[7] # progress to "dx in progress"
        chg_arr[8,14] += txrate(g.intopt,0,0,0,0)[0]*X[8] # new
        chg_arr[9,15] += txrate(g.intopt,1,1,0,0)[0]*X[9]
        chg_arr[10,16] += txrate(g.intopt,0,1,0,0)[0]*X[10]
        chg_arr[11,17] += txrate(g.intopt,1,2,0,0)[0]*X[11]
        chg_arr[12,18] += txrate(g.intopt,0,2,0,0)[0]*X[12]
        chg_arr[57,63] += txrate(g.intopt,1,0,1,0)[0]*X[57] # progress to "dx in progress"
        chg_arr[58,64] += txrate(g.intopt,0,0,1,0)[0]*X[58] # retreat
        chg_arr[59,65] += txrate(g.intopt,1,1,1,0)[0]*X[59]
        chg_arr[60,66] += txrate(g.intopt,0,1,1,0)[0]*X[60]
        chg_arr[61,67] += txrate(g.intopt,1,2,1,0)[0]*X[61]
        chg_arr[62,68] += txrate(g.intopt,0,2,1,0)[0]*X[62]
        chg_arr[7,19] += txrate(g.intopt,1,0,0,0)[2]*(1-g.acq_s)*X[7] # progress to "inapp tx"
        chg_arr[8,20] += txrate(g.intopt,0,0,0,0)[2]*(1-g.acq_s)*X[8] # new
        chg_arr[9,21] += txrate(g.intopt,1,1,0,0)[2]*(1-g.acq_i1)*X[9]
        chg_arr[10,22] += txrate(g.intopt,0,1,0,0)[2]*(1-g.acq_i1)*X[10]
        chg_arr[11,23] += txrate(g.intopt,1,2,0,0)[2]*X[11]
        chg_arr[12,24] += txrate(g.intopt,0,2,0,0)[2]*X[12]
        chg_arr[7,21] += txrate(g.intopt,1,0,0,0)[2]*(g.acq_s)*(1-g.acq_mdr)*X[7] # new resistance
        chg_arr[8,22] += txrate(g.intopt,0,0,0,0)[2]*(g.acq_s)*(1-g.acq_mdr)*X[8]
        chg_arr[7,23] += txrate(g.intopt,1,0,0,0)[2]*(g.acq_s)*(g.acq_mdr)*X[7] # new resistance
        chg_arr[8,24] += txrate(g.intopt,0,0,0,0)[2]*(g.acq_s)*(g.acq_mdr)*X[8]
        chg_arr[9,23] += txrate(g.intopt,1,1,0,0)[2]*(g.acq_i1)*X[9]
        chg_arr[10,24] += txrate(g.intopt,0,1,0,0)[2]*(g.acq_i1)*X[10]
        chg_arr[57,69] += txrate(g.intopt,1,0,1,0)[2]*(1-g.acq_s)*X[57] # progress to "inapp tx"
        chg_arr[58,70] += txrate(g.intopt,0,0,1,0)[2]*(1-g.acq_s)*X[58] # retreat
        chg_arr[59,71] += txrate(g.intopt,1,1,1,0)[2]*(1-g.acq_i2)*X[59]
        chg_arr[60,72] += txrate(g.intopt,0,1,1,0)[2]*(1-g.acq_i2)*X[60]
        chg_arr[61,73] += txrate(g.intopt,1,2,1,0)[2]*X[61]
        chg_arr[62,74] += txrate(g.intopt,0,2,1,0)[2]*X[62]
        chg_arr[57,71] += txrate(g.intopt,1,0,1,0)[2]*(g.acq_s)*(1-g.acq_mdr)*X[57] # new resistance
        chg_arr[58,72] += txrate(g.intopt,0,0,1,0)[2]*(g.acq_s)*(1-g.acq_mdr)*X[58]
        chg_arr[57,73] += txrate(g.intopt,1,0,1,0)[2]*(g.acq_s)*(g.acq_mdr)*X[57] # new resistance
        chg_arr[58,74] += txrate(g.intopt,0,0,1,0)[2]*(g.acq_s)*(g.acq_mdr)*X[58]
        chg_arr[59,73] += txrate(g.intopt,1,1,1,0)[2]*(g.acq_i2)*X[59]
        chg_arr[60,74] += txrate(g.intopt,0,1,1,0)[2]*(g.acq_i2)*X[60]
        chg_arr[13,51] += txrate(g.intopt,1,0,0,0)[1]*X[13] # progress to treated
        chg_arr[14,51] += txrate(g.intopt,0,0,0,0)[1]*X[14] # new
        chg_arr[15,52] += txrate(g.intopt,1,1,0,0)[1]*X[15]
        chg_arr[16,52] += txrate(g.intopt,0,1,0,0)[1]*X[16]
        chg_arr[17,53] += txrate(g.intopt,1,2,0,0)[1]*X[17]
        chg_arr[18,53] += txrate(g.intopt,0,2,0,0)[1]*X[18]
        chg_arr[63,51] += txrate(g.intopt,1,0,1,0)[1]*X[63] # progress to treated
        chg_arr[64,51] += txrate(g.intopt,0,0,1,0)[1]*X[64] # retreat
        chg_arr[65,52] += txrate(g.intopt,1,1,1,0)[1]*X[65]
        chg_arr[66,52] += txrate(g.intopt,0,1,1,0)[1]*X[66]
        chg_arr[67,53] += txrate(g.intopt,1,2,1,0)[1]*X[67]
        chg_arr[68,53] += txrate(g.intopt,0,2,1,0)[1]*X[68]

        chg_arr[7+25,13+25] += txrate(g.intopt,1,0,0,1)[0]*X[7+25] # progress to "dx in progress"
        chg_arr[8+25,14+25] += txrate(g.intopt,0,0,0,1)[0]*X[8+25] # new
        chg_arr[9+25,15+25] += txrate(g.intopt,1,1,0,1)[0]*X[9+25]
        chg_arr[10+25,16+25] += txrate(g.intopt,0,1,0,1)[0]*X[10+25]
        chg_arr[11+25,17+25] += txrate(g.intopt,1,2,0,1)[0]*X[11+25]
        chg_arr[12+25,18+25] += txrate(g.intopt,0,2,0,1)[0]*X[12+25]
        chg_arr[57+25,63+25] += txrate(g.intopt,1,0,1,1)[0]*X[57+25] # progress to "dx in progress"
        chg_arr[58+25,64+25] += txrate(g.intopt,0,0,1,1)[0]*X[58+25] # retreat
        chg_arr[59+25,65+25] += txrate(g.intopt,1,1,1,1)[0]*X[59+25]
        chg_arr[60+25,66+25] += txrate(g.intopt,0,1,1,1)[0]*X[60+25]
        chg_arr[61+25,67+25] += txrate(g.intopt,1,2,1,1)[0]*X[61+25]
        chg_arr[62+25,68+25] += txrate(g.intopt,0,2,1,1)[0]*X[62+25]
        chg_arr[7+25,19+25] += txrate(g.intopt,1,0,0,1)[2]*(1-g.acq_s)*X[7+25] # progress to "inapp tx"
        chg_arr[8+25,20+25] += txrate(g.intopt,0,0,0,1)[2]*(1-g.acq_s)*X[8+25] # new
        chg_arr[9+25,21+25] += txrate(g.intopt,1,1,0,1)[2]*(1-g.acq_i1)*X[9+25]
        chg_arr[10+25,22+25] += txrate(g.intopt,0,1,0,1)[2]*(1-g.acq_i1)*X[10+25]
        chg_arr[11+25,23+25] += txrate(g.intopt,1,2,0,1)[2]*X[11+25]
        chg_arr[12+25,24+25] += txrate(g.intopt,0,2,0,1)[2]*X[12+25]
        chg_arr[7+25,21+25] += txrate(g.intopt,1,0,0,1)[2]*(g.acq_s)*(1-g.acq_mdr)*X[7+25] # new resistance
        chg_arr[8+25,22+25] += txrate(g.intopt,0,0,0,1)[2]*(g.acq_s)*(1-g.acq_mdr)*X[8+25]
        chg_arr[7+25,23+25] += txrate(g.intopt,1,0,0,1)[2]*(g.acq_s)*(g.acq_mdr)*X[7+25] # new resistance
        chg_arr[8+25,24+25] += txrate(g.intopt,0,0,0,1)[2]*(g.acq_s)*(g.acq_mdr)*X[8+25]
        chg_arr[9+25,23+25] += txrate(g.intopt,1,1,0,1)[2]*(g.acq_i1)*X[9+25]
        chg_arr[10+25,24+25] += txrate(g.intopt,0,1,0,1)[2]*(g.acq_i1)*X[10+25]
        chg_arr[57+25,69+25] += txrate(g.intopt,1,0,1,1)[2]*(1-g.acq_s)*X[57+25] # progress to "inapp tx"
        chg_arr[58+25,70+25] += txrate(g.intopt,0,0,1,1)[2]*(1-g.acq_s)*X[58+25] # retreat
        chg_arr[59+25,71+25] += txrate(g.intopt,1,1,1,1)[2]*(1-g.acq_i2)*X[59+25]
        chg_arr[60+25,72+25] += txrate(g.intopt,0,1,1,1)[2]*(1-g.acq_i2)*X[60+25]
        chg_arr[61+25,73+25] += txrate(g.intopt,1,2,1,1)[2]*X[61+25]
        chg_arr[62+25,74+25] += txrate(g.intopt,0,2,1,1)[2]*X[62+25]
        chg_arr[57+25,71+25] += txrate(g.intopt,1,0,1,1)[2]*(g.acq_s)*(1-g.acq_mdr)*X[57+25] # new resistance
        chg_arr[58+25,72+25] += txrate(g.intopt,0,0,1,1)[2]*(g.acq_s)*(1-g.acq_mdr)*X[58+25]
        chg_arr[57+25,73+25] += txrate(g.intopt,1,0,1,1)[2]*(g.acq_s)*(g.acq_mdr)*X[57+25] # new resistance
        chg_arr[58+25,74+25] += txrate(g.intopt,0,0,1,1)[2]*(g.acq_s)*(g.acq_mdr)*X[58+25]
        chg_arr[59+25,73+25] += txrate(g.intopt,1,1,1,1)[2]*(g.acq_i2)*X[59+25]
        chg_arr[60+25,74+25] += txrate(g.intopt,0,1,1,1)[2]*(g.acq_i2)*X[60+25]
        chg_arr[13+25,51+25] += txrate(g.intopt,1,0,0,1)[1]*X[13+25] # progress to treated
        chg_arr[14+25,51+25] += txrate(g.intopt,0,0,0,1)[1]*X[14+25] # new
        chg_arr[15+25,52+25] += txrate(g.intopt,1,1,0,1)[1]*X[15+25]
        chg_arr[16+25,52+25] += txrate(g.intopt,0,1,0,1)[1]*X[16+25]
        chg_arr[17+25,53+25] += txrate(g.intopt,1,2,0,1)[1]*X[17+25]
        chg_arr[18+25,53+25] += txrate(g.intopt,0,2,0,1)[1]*X[18+25]
        chg_arr[63+25,51+25] += txrate(g.intopt,1,0,1,1)[1]*X[63+25] # progress to treated
        chg_arr[64+25,51+25] += txrate(g.intopt,0,0,1,1)[1]*X[64+25] # retreat
        chg_arr[65+25,52+25] += txrate(g.intopt,1,1,1,1)[1]*X[65+25]
        chg_arr[66+25,52+25] += txrate(g.intopt,0,1,1,1)[1]*X[66+25]
        chg_arr[67+25,53+25] += txrate(g.intopt,1,2,1,1)[1]*X[67+25]
        chg_arr[68+25,53+25] += txrate(g.intopt,0,2,1,1)[1]*X[68+25]

        chg_arr[i4+19,i4+57] += (12./g.dur_fail)*g.defvfail_s*X[i4+19] # default from "inapp tx" (new cases)
        chg_arr[i4+20,i4+58] += (12./g.dur_fail)*g.defvfail_s*X[i4+20] # assume zero cost b/c 1/2 of default is cured
        chg_arr[i4+21,i4+59] += (12./g.dur_fail)*g.defvfail_i*X[i4+21]
        chg_arr[i4+22,i4+60] += (12./g.dur_fail)*g.defvfail_i*X[i4+22]
        chg_arr[i4+23,i4+61] += (12./g.dur_fail)*g.defvfail_m*X[i4+23]
        chg_arr[i4+24,i4+62] += (12./g.dur_fail)*g.defvfail_m*X[i4+24]
        chg_arr[i5+19,i5+7] += (12./g.dur_fail)*g.defvfail_s*X[i5+19] # default from "inapp tx" (retx cases)
        chg_arr[i5+20,i5+8] += (12./g.dur_fail)*g.defvfail_s*X[i5+20] # assume zero cost b/c 1/2 of default is cured
        chg_arr[i5+21,i5+9] += (12./g.dur_fail)*g.defvfail_i*X[i5+21]
        chg_arr[i5+22,i5+10] += (12./g.dur_fail)*g.defvfail_i*X[i5+22]
        chg_arr[i5+23,i5+11] += (12./g.dur_fail)*g.defvfail_m*X[i5+23]
        chg_arr[i5+24,i5+12] += (12./g.dur_fail)*g.defvfail_m*X[i5+24]
        chg_arr[i4+19,i4+51] += (12./g.dur_fail)*(1-g.defvfail_s)*X[i4+19]*(1-g.fail_rt) # failure from "inapp tx" (new TB cases)
        chg_arr[i4+20,i4+51] += (12./g.dur_fail)*(1-g.defvfail_s)*X[i4+20]*(1-g.fail_rt) # these are the pts that are then successfully treated
        chg_arr[i4+21,i4+52] += (12./g.dur_fail)*(1-g.defvfail_i)*X[i4+21]*(1-g.fail_i2) # the ones below re-fail
        chg_arr[i4+22,i4+52] += (12./g.dur_fail)*(1-g.defvfail_i)*X[i4+22]*(1-g.fail_i2)
        chg_arr[i4+23,i4+53] += (12./g.dur_fail)*(1-g.defvfail_m)*X[i4+23]*(1-g.fail_m2)
        chg_arr[i4+24,i4+53] += (12./g.dur_fail)*(1-g.defvfail_m)*X[i4+24]*(1-g.fail_m2)
        chg_arr[i4+19,i4+69] += (12./g.dur_fail)*(1-g.defvfail_s)*X[i4+19]*(g.fail_rt) # re-failure from "inapp tx" -> now classified as retx
        chg_arr[i4+20,i4+70] += (12./g.dur_fail)*(1-g.defvfail_s)*X[i4+20]*(g.fail_rt)
        chg_arr[i4+21,i4+71] += (12./g.dur_fail)*(1-g.defvfail_i)*X[i4+21]*(g.fail_i2)
        chg_arr[i4+22,i4+72] += (12./g.dur_fail)*(1-g.defvfail_i)*X[i4+22]*(g.fail_i2)
        chg_arr[i4+23,i4+73] += (12./g.dur_fail)*(1-g.defvfail_m)*X[i4+23]*(g.fail_m2)
        chg_arr[i4+24,i4+74] += (12./g.dur_fail)*(1-g.defvfail_m)*X[i4+24]*(g.fail_m2)
        chg_arr[i5+19,i5+1] += (12./g.dur_fail)*(1-g.defvfail_s)*X[i5+19]*(1-g.fail_rt) # failure from "inapp tx"
        chg_arr[i5+20,i5+1] += (12./g.dur_fail)*(1-g.defvfail_s)*X[i5+20]*(1-g.fail_rt) # retx pts, so "re-fail" stays in the same box
        chg_arr[i5+21,i5+2] += (12./g.dur_fail)*(1-g.defvfail_i)*X[i5+21]*(1-g.fail_i2)
        chg_arr[i5+22,i5+2] += (12./g.dur_fail)*(1-g.defvfail_i)*X[i5+22]*(1-g.fail_i2)
        chg_arr[i5+23,i5+3] += (12./g.dur_fail)*(1-g.defvfail_m)*X[i5+23]*(1-g.fail_m2)
        chg_arr[i5+24,i5+3] += (12./g.dur_fail)*(1-g.defvfail_m)*X[i5+24]*(1-g.fail_m2)
        i6 = np.arange(25)                          # HIV infection
        chg_arr[i6,i6+25] += X[25]*X[i6]
        chg_arr[i6+50,i6+75] += X[25]*X[i6+50]
        chg_arr[50,75]-=X[25]*X[50]
        chg_arr[50,75]+=X[25]*rec50
        i7 = np.arange(4)                           # inappropriate tx for TB
        chg_arr[i7,i7+50] += g.prop_cough*(1-invspec(g.intopt,0,0))*X[i7]
        chg_arr[i7+25,i7+75] += g.prop_cough*(1-invspec(g.intopt,1,0))*X[i7+25]
        chg_arr[25,75] -= g.prop_cough*(1-invspec(g.intopt,1,0))*X[25]
        chg_arr[25,75] += g.prop_cough*(1-invspec(g.intopt,1,0))*rec25
        i8 = np.arange(50)                         # mortality
        chg_arr[i6,0] += mortality(0,0)*X[i6]       # background mortality
        chg_arr[i6+50,0] += mortality(0,0)*X[i6+50]
        chg_arr[50,0]-=mortality(0,0)*X[50]
        chg_arr[50,0]+=mortality(0,0)*rec50
        chg_arr[i6+25,0] += mortality(0,1)*X[i6+25] # HIV deaths
        chg_arr[i6+75,0] += mortality(0,1)*X[i6+75]
        chg_arr[75,0]-=mortality(0,1)*X[75]
        chg_arr[75,0]+=mortality(0,1)*rec75
        chg_arr[25,0]-=mortality(0,1)*X[25]
        chg_arr[25,0]+=mortality(0,1)*rec25
        i9 = np.array([7,9,11,13,15,17])  # smear-pos
        i10 = np.array([4,5,6,8,10,12,14,16,18,19,20,21,22,23,24])  # smear-neg
        chg_arr[i9,0] += mortality(1,0)*X[i9]              # smear-positive TB deaths
        chg_arr[i9+25,0] += mortality(1,1)*X[i9+25]
        chg_arr[i9+50,0] += mortality(1,0)*X[i9+50]
        chg_arr[i9+75,0] += mortality(1,1)*X[i9+75]
        chg_arr[i10,0] += mortality(2,0)*X[i10]              # smear-negative TB deaths
        chg_arr[i10+25,0] += mortality(2,1)*X[i10+25]
        chg_arr[i10+50,0] += mortality(2,0)*X[i10+50]
        chg_arr[i10+75,0] += mortality(2,1)*X[i10+75]
        chg_arr[i1+7,i1+1] += g.cure_sp*X[i1+7]              # smear-positive self-cure (no HIV)
        chg_arr[i1+9,i1+2] += g.cure_sp*X[i1+9]
        chg_arr[i1+11,i1+3] += g.cure_sp*X[i1+11]
        chg_arr[i1+13,i1+1] += g.cure_sp*X[i1+13]
        chg_arr[i1+15,i1+2] += g.cure_sp*X[i1+15]
        chg_arr[i1+17,i1+3] += g.cure_sp*X[i1+17]
        chg_arr[i1+8,i1+1] += g.cure_sn*X[i1+8]              # smear-negative self-cure (no HIV)
        chg_arr[i1+10,i1+2] += g.cure_sn*X[i1+10]
        chg_arr[i1+12,i1+3] += g.cure_sn*X[i1+12]
        chg_arr[i1+14,i1+1] += g.cure_sn*X[i1+14]
        chg_arr[i1+16,i1+2] += g.cure_sn*X[i1+16]
        chg_arr[i1+18,i1+3] += g.cure_sn*X[i1+18]
        chg_arr[i1+4,i1+1] += g.cure_sn*X[i1+4]
        chg_arr[i1+5,i1+2] += g.cure_sn*X[i1+5]
        chg_arr[i1+6,i1+3] += g.cure_sn*X[i1+6]
        chg_arr[i1+19,i1+1] += g.cure_sn*X[i1+19]
        chg_arr[i1+21,i1+2] += g.cure_sn*X[i1+21]
        chg_arr[i1+23,i1+3] += g.cure_sn*X[i1+23]
        chg_arr[i1+20,i1+1] += g.cure_sn*X[i1+20]
        chg_arr[i1+22,i1+2] += g.cure_sn*X[i1+22]
        chg_arr[i1+24,i1+3] += g.cure_sn*X[i1+24]

        # calculate incidence (to ensure that it reaches the target)
        incidence = np.sum(chg_arr[i4,i4+4]) + np.sum(chg_arr[i4,i4+5]) + np.sum(chg_arr[i4,i4+6]) + \
                 np.sum(chg_arr[i4+1,i4+4]) + np.sum(chg_arr[i4+1,i4+5]) + np.sum(chg_arr[i4+1,i4+6]) + \
                 np.sum(chg_arr[i4+2,i4+4]) + np.sum(chg_arr[i4+2,i4+5]) + np.sum(chg_arr[i4+2,i4+6]) + \
                 np.sum(chg_arr[i4+3,i4+4]) + np.sum(chg_arr[i4+3,i4+5]) + np.sum(chg_arr[i4+3,i4+6]) + \
                 np.sum(chg_arr[i5,i5+4]) + np.sum(chg_arr[i5,i5+5]) + np.sum(chg_arr[i5,i5+6]) + \
                 np.sum(chg_arr[i5+1,i5+4]) + np.sum(chg_arr[i5+1,i5+5]) + np.sum(chg_arr[i5+1,i5+6]) + \
                 np.sum(chg_arr[i5+2,i5+4]) + np.sum(chg_arr[i5+2,i5+5]) + np.sum(chg_arr[i5+2,i5+6]) + \
                 np.sum(chg_arr[i5+3,i5+4]) + np.sum(chg_arr[i5+3,i5+5]) + np.sum(chg_arr[i5+3,i5+6]) + \
                 np.sum(chg_arr[i4+19,i4+57]) + np.sum(chg_arr[i4+20,i4+58]) + np.sum(chg_arr[i4+21,i4+59]) + \
                 np.sum(chg_arr[i4+22,i4+60]) + np.sum(chg_arr[i4+23,i4+61]) + np.sum(chg_arr[i4+24,i4+62]) + \
                 np.sum(chg_arr[i5+19,i5+7]) + np.sum(chg_arr[i5+20,i5+8]) + np.sum(chg_arr[i5+21,i5+9]) + \
                 np.sum(chg_arr[i5+22,i5+10]) + np.sum(chg_arr[i5+23,i5+11]) + np.sum(chg_arr[i5+24,i5+12]) + \
                 np.sum((12./g.dur_fail)*(1-g.defvfail_s)*X[i3+19]) + np.sum((12./g.dur_fail)*(1-g.defvfail_s)*X[i3+20]) + \
                 np.sum((12./g.dur_fail)*(1-g.defvfail_i)*X[i3+21]) + np.sum((12./g.dur_fail)*(1-g.defvfail_i)*X[i3+22]) + \
                 np.sum((12./g.dur_fail)*(1-g.defvfail_m)*X[i3+23]) + np.sum((12./g.dur_fail)*(1-g.defvfail_m)*X[i3+24])


        # calculate incidence of non-retreatment TB (to match MDR-TB incidence to target)
        incnew = np.sum(chg_arr[i4,i4+4]) + np.sum(chg_arr[i4,i4+5]) + np.sum(chg_arr[i4,i4+6]) + \
                 np.sum(chg_arr[i4+1,i4+4]) + np.sum(chg_arr[i4+1,i4+5]) + np.sum(chg_arr[i4+1,i4+6]) + \
                 np.sum(chg_arr[i4+2,i4+4]) + np.sum(chg_arr[i4+2,i4+5]) + np.sum(chg_arr[i4+2,i4+6]) + \
                 np.sum(chg_arr[i4+3,i4+4]) + np.sum(chg_arr[i4+3,i4+5]) + np.sum(chg_arr[i4+3,i4+6])

        # calculate incidence of new MDR-TB
        incmdrnew = np.sum(chg_arr[i4,i4+6]) + np.sum(chg_arr[i4+1,i4+6]) + np.sum(chg_arr[i4+2,i4+6]) + np.sum(chg_arr[i4+3,i4+6])

        # calculate HIV prevalence (to match to target)
        i12 = np.concatenate((np.arange(25,50), np.arange(75,100)))
        HIVprev = (sum(X[i12])-X[25]-X[75]+rec25+rec75)/1000

        # differential equations:
        dxdt[0] = np.sum(chg_arr[:,0])-np.sum(chg_arr[0,:])
        dxdt[1] = np.sum(chg_arr[:,1])-np.sum(chg_arr[1,:])
        dxdt[2] = np.sum(chg_arr[:,2])-np.sum(chg_arr[2,:])
        dxdt[3] = np.sum(chg_arr[:,3])-np.sum(chg_arr[3,:])
        dxdt[4] = np.sum(chg_arr[:,4])-np.sum(chg_arr[4,:])
        dxdt[5] = np.sum(chg_arr[:,5])-np.sum(chg_arr[5,:])
        dxdt[6] = np.sum(chg_arr[:,6])-np.sum(chg_arr[6,:])
        dxdt[7] = np.sum(chg_arr[:,7])-np.sum(chg_arr[7,:])
        dxdt[8] = np.sum(chg_arr[:,8])-np.sum(chg_arr[8,:])
        dxdt[9] = np.sum(chg_arr[:,9])-np.sum(chg_arr[9,:])
        dxdt[10] = np.sum(chg_arr[:,10])-np.sum(chg_arr[10,:])
        dxdt[11] = np.sum(chg_arr[:,11])-np.sum(chg_arr[11,:])
        dxdt[12] = np.sum(chg_arr[:,12])-np.sum(chg_arr[12,:])
        dxdt[13] = np.sum(chg_arr[:,13])-np.sum(chg_arr[13,:])
        dxdt[14] = np.sum(chg_arr[:,14])-np.sum(chg_arr[14,:])
        dxdt[15] = np.sum(chg_arr[:,15])-np.sum(chg_arr[15,:])
        dxdt[16] = np.sum(chg_arr[:,16])-np.sum(chg_arr[16,:])
        dxdt[17] = np.sum(chg_arr[:,17])-np.sum(chg_arr[17,:])
        dxdt[18] = np.sum(chg_arr[:,18])-np.sum(chg_arr[18,:])
        dxdt[19] = np.sum(chg_arr[:,19])-np.sum(chg_arr[19,:])
        dxdt[20] = np.sum(chg_arr[:,20])-np.sum(chg_arr[20,:])
        dxdt[21] = np.sum(chg_arr[:,21])-np.sum(chg_arr[21,:])
        dxdt[22] = np.sum(chg_arr[:,22])-np.sum(chg_arr[22,:])
        dxdt[23] = np.sum(chg_arr[:,23])-np.sum(chg_arr[23,:])
        dxdt[24] = np.sum(chg_arr[:,24])-np.sum(chg_arr[24,:])
        dxdt[25] = g.target_hiv - HIVprev # equals zero when HIV prevalence is at target
        dxdt[26] = np.sum(chg_arr[:,26])-np.sum(chg_arr[26,:])
        dxdt[27] = np.sum(chg_arr[:,27])-np.sum(chg_arr[27,:])
        dxdt[28] = np.sum(chg_arr[:,28])-np.sum(chg_arr[28,:])
        dxdt[29] = np.sum(chg_arr[:,29])-np.sum(chg_arr[29,:])
        dxdt[30] = np.sum(chg_arr[:,30])-np.sum(chg_arr[30,:])
        dxdt[31] = np.sum(chg_arr[:,31])-np.sum(chg_arr[31,:])
        dxdt[32] = np.sum(chg_arr[:,32])-np.sum(chg_arr[32,:])
        dxdt[33] = np.sum(chg_arr[:,33])-np.sum(chg_arr[33,:])
        dxdt[34] = np.sum(chg_arr[:,34])-np.sum(chg_arr[34,:])
        dxdt[35] = np.sum(chg_arr[:,35])-np.sum(chg_arr[35,:])
        dxdt[36] = np.sum(chg_arr[:,36])-np.sum(chg_arr[36,:])
        dxdt[37] = np.sum(chg_arr[:,37])-np.sum(chg_arr[37,:])
        dxdt[38] = np.sum(chg_arr[:,38])-np.sum(chg_arr[38,:])
        dxdt[39] = np.sum(chg_arr[:,39])-np.sum(chg_arr[39,:])
        dxdt[40] = np.sum(chg_arr[:,40])-np.sum(chg_arr[40,:])
        dxdt[41] = np.sum(chg_arr[:,41])-np.sum(chg_arr[41,:])
        dxdt[42] = np.sum(chg_arr[:,42])-np.sum(chg_arr[42,:])
        dxdt[43] = np.sum(chg_arr[:,43])-np.sum(chg_arr[43,:])
        dxdt[44] = np.sum(chg_arr[:,44])-np.sum(chg_arr[44,:])
        dxdt[45] = np.sum(chg_arr[:,45])-np.sum(chg_arr[45,:])
        dxdt[46] = np.sum(chg_arr[:,46])-np.sum(chg_arr[46,:])
        dxdt[47] = np.sum(chg_arr[:,47])-np.sum(chg_arr[47,:])
        dxdt[48] = np.sum(chg_arr[:,48])-np.sum(chg_arr[48,:])
        dxdt[49] = np.sum(chg_arr[:,49])-np.sum(chg_arr[49,:])
        dxdt[50] = g.target_inc-incidence # equals zero when TB incidence is at target
        dxdt[51] = np.sum(chg_arr[:,51])-np.sum(chg_arr[51,:])
        dxdt[52] = np.sum(chg_arr[:,52])-np.sum(chg_arr[52,:])
        dxdt[53] = np.sum(chg_arr[:,53])-np.sum(chg_arr[53,:])
        dxdt[54] = np.sum(chg_arr[:,54])-np.sum(chg_arr[54,:])
        dxdt[55] = np.sum(chg_arr[:,55])-np.sum(chg_arr[55,:])
        dxdt[56] = np.sum(chg_arr[:,56])-np.sum(chg_arr[56,:])
        dxdt[57] = np.sum(chg_arr[:,57])-np.sum(chg_arr[57,:])
        dxdt[58] = np.sum(chg_arr[:,58])-np.sum(chg_arr[58,:])
        dxdt[59] = np.sum(chg_arr[:,59])-np.sum(chg_arr[59,:])
        dxdt[60] = np.sum(chg_arr[:,60])-np.sum(chg_arr[60,:])
        dxdt[61] = np.sum(chg_arr[:,61])-np.sum(chg_arr[61,:])
        dxdt[62] = np.sum(chg_arr[:,62])-np.sum(chg_arr[62,:])
        dxdt[63] = np.sum(chg_arr[:,63])-np.sum(chg_arr[63,:])
        dxdt[64] = np.sum(chg_arr[:,64])-np.sum(chg_arr[64,:])
        dxdt[65] = np.sum(chg_arr[:,65])-np.sum(chg_arr[65,:])
        dxdt[66] = np.sum(chg_arr[:,66])-np.sum(chg_arr[66,:])
        dxdt[67] = np.sum(chg_arr[:,67])-np.sum(chg_arr[67,:])
        dxdt[68] = np.sum(chg_arr[:,68])-np.sum(chg_arr[68,:])
        dxdt[69] = np.sum(chg_arr[:,69])-np.sum(chg_arr[69,:])
        dxdt[70] = np.sum(chg_arr[:,70])-np.sum(chg_arr[70,:])
        dxdt[71] = np.sum(chg_arr[:,71])-np.sum(chg_arr[71,:])
        dxdt[72] = np.sum(chg_arr[:,72])-np.sum(chg_arr[72,:])
        dxdt[73] = np.sum(chg_arr[:,73])-np.sum(chg_arr[73,:])
        dxdt[74] = np.sum(chg_arr[:,74])-np.sum(chg_arr[74,:])
        dxdt[75] = g.target_mdr-incmdrnew/incnew # equals zero when new MDR-TB prevalence is at target
        dxdt[76] = np.sum(chg_arr[:,76])-np.sum(chg_arr[76,:])
        dxdt[77] = np.sum(chg_arr[:,77])-np.sum(chg_arr[77,:])
        dxdt[78] = np.sum(chg_arr[:,78])-np.sum(chg_arr[78,:])
        dxdt[79] = np.sum(chg_arr[:,79])-np.sum(chg_arr[79,:])
        dxdt[80] = np.sum(chg_arr[:,80])-np.sum(chg_arr[80,:])
        dxdt[81] = np.sum(chg_arr[:,81])-np.sum(chg_arr[81,:])
        dxdt[82] = np.sum(chg_arr[:,82])-np.sum(chg_arr[82,:])
        dxdt[83] = np.sum(chg_arr[:,83])-np.sum(chg_arr[83,:])
        dxdt[84] = np.sum(chg_arr[:,84])-np.sum(chg_arr[84,:])
        dxdt[85] = np.sum(chg_arr[:,85])-np.sum(chg_arr[85,:])
        dxdt[86] = np.sum(chg_arr[:,86])-np.sum(chg_arr[86,:])
        dxdt[87] = np.sum(chg_arr[:,87])-np.sum(chg_arr[87,:])
        dxdt[88] = np.sum(chg_arr[:,88])-np.sum(chg_arr[88,:])
        dxdt[89] = np.sum(chg_arr[:,89])-np.sum(chg_arr[89,:])
        dxdt[90] = np.sum(chg_arr[:,90])-np.sum(chg_arr[90,:])
        dxdt[91] = np.sum(chg_arr[:,91])-np.sum(chg_arr[91,:])
        dxdt[92] = np.sum(chg_arr[:,92])-np.sum(chg_arr[92,:])
        dxdt[93] = np.sum(chg_arr[:,93])-np.sum(chg_arr[93,:])
        dxdt[94] = np.sum(chg_arr[:,94])-np.sum(chg_arr[94,:])
        dxdt[95] = np.sum(chg_arr[:,95])-np.sum(chg_arr[95,:])
        dxdt[96] = np.sum(chg_arr[:,96])-np.sum(chg_arr[96,:])
        dxdt[97] = np.sum(chg_arr[:,97])-np.sum(chg_arr[97,:])
        dxdt[98] = np.sum(chg_arr[:,98])-np.sum(chg_arr[98,:])
        dxdt[99] = np.sum(chg_arr[:,99])-np.sum(chg_arr[99,:])
        dxdt[100] = np.sum(chg_arr[:,50])-np.sum(chg_arr[50,:])
        dxdt[101] = np.sum(chg_arr[:,25])-np.sum(chg_arr[25,:])
        return dxdt

# The initial population for the "solvebeta" function requires initial
#    estimates for the populations of compartments 100 and 101.
    I2 = np.append(I,[I[25]-I[75]])
    I3 = np.append(I2,[I[25]-I[50]])
# As above, substitute compartments 25, 50, and 75 with the quantities
#     of interest
    I3[25] = g.hiv_inc
    I3[50] = g.beta_s
    I3[75] = g.fit_m

# solve for the roots of the system of equations.
# compartment 25 is now hiv_inc, 50 beta_s, and 75 fit_m at steady-state.
    equipop_pre = fsolve(solvebeta,I3)
    g.beta_s = equipop_pre[50]
    g.fit_m = 1.-equipop_pre[75]
    g.fit_i = 1.-equipop_pre[75]*0.25
    g.hiv_inc = equipop_pre[25]

# now, add back in the actual population sizes to compartments 25, 50, and 75.
    equipop_pre[25] = (100000-sum(equipop_pre)+equipop_pre[25]+equipop_pre[50] + \
                             equipop_pre[75]+2.*equipop_pre[100]+2.*equipop_pre[101])/3.
    equipop_pre[75] = equipop_pre[25]-equipop_pre[100]
    equipop_pre[50] = equipop_pre[25]-equipop_pre[101]

# "equipop" is now the equilibrium population that fits the user-specified values
# of TB incidence, MDR-TB incidence, and HIV prevalence.
    equipop = np.zeros(g.tb_num)
    equipop[:] = equipop_pre[0:100]

###########################################################

# Check to make sure no negative cells - will stop running if they are detected:
    for xyz in range(100):
        if equipop[xyz]<0:
            print "negative population - recalibrating..."
            hb_log("negative population - recalibrating...")

            I3[4:25] +=10.
            I3[29:50] +=10.
            I3[54:75]+=10.
            I3[79:100]+=10.
            I3[100] = I[25]-I[75]
            I3[101] = I[25]-I[50]
            I3[0] = 100000 - np.sum(I3[1:100])
            equipop_pre = fsolve(solvebeta,I3)
            # print equipop
            print "beta", equipop_pre[50]
            print "MDR fitness", 1.-equipop_pre[75]
            print "INH fitness", 1.-equipop_pre[75]*0.25
            print "HIV incidence", equipop_pre[25]
            hb_log( "beta: {}".format(equipop_pre[50]))
            hb_log( "MDR fitness: {}".format( 1.-equipop_pre[75]))
            hb_log( "INH fitness: {}".format( 1.-equipop_pre[75]*0.25))
            hb_log( "HIV incidence: {}".format(equipop_pre[25]))

            beta_s = equipop_pre[50]
            fit_m = 1.-equipop_pre[75]
            fit_i = 1.-equipop_pre[75]*0.25
            hiv_inc = equipop_pre[25]
            # i10 = np.array([17,18,42,43,68,93])
            # equipop_pre[i10] = 0.
            equipop_pre[25] = (100000-np.sum(equipop_pre)+equipop_pre[25]+equipop_pre[50] + \
                                     equipop_pre[75]+2.*equipop_pre[100]+2.*equipop_pre[101])/3.
            equipop_pre[75] = equipop_pre[25]-equipop_pre[100]
            equipop_pre[50] = equipop_pre[25]-equipop_pre[101]
            equipop = np.zeros(tb_num)
            equipop[:] = equipop_pre[0:100]
            print "total pop (check for 100,000)", sum(equipop[:])
            hb_log("total pop (check for 100,000): {}".format(sum(equipop[:])))
            print " "

    for xyz2 in range(100):
        if equipop[xyz2]<0:
            print "negative population - failed run"
            hb_log("negative population - failed run")
            exit(0)
            

###########################################################
# 8. CALCULATE INCIDENCE, PREVALENCE, MORTALITY, COST
###########################################################
# This function takes the same equations as in 6 and 7 above
# and uses them to calculate incidence (according to treatment & HIV status),
# as well as TB and HIV prevalence, TB mortality, and total costs.

    def incprevmort(pop1):

        #Global vars used here:
        intopt = g.intopt
        hiv_inc = g.hiv_inc
        tb_num = g.tb_num
        beta_s = g.beta_s
        relbeta_sn = g.relbeta_sn
        fit_i = g.fit_i
        fit_m = g.fit_m
        rapid = g.rapid
        rapid_h = g.rapid_h
        react = g.react
        react_h = g.react_h
        prot = g.prot
        predx_dur = g.predx_dur
        predx_dur_h = g.predx_dur_h
        cure_sn = g.cure_sn
        prop_inf = g.prop_inf
        prop_inf_h = g.prop_inf_h
        acq_s = g.acq_s
        acq_i1 = g.acq_i1
        acq_i2 = g.acq_i2
        acq_mdr = g.acq_mdr
        dur_fail = g.dur_fail
        defvfail_s = g.defvfail_s
        defvfail_i = g.defvfail_i
        defvfail_m = g.defvfail_m
        fail_rt = g.fail_rt
        fail_i2 = g.fail_i2
        fail_m2 = g.fail_m2
        prop_cough = g.prop_cough
        cure_sp = g.cure_sp
        dx_rate = g.dx_rate



        dxdt = np.zeros(tb_num)  # create the vector of ODEs
        X = pop1

        active_s = np.array([7,13])              # smear-pos, DS-TB compartments
        active_s1 = active_s + 25
        active_s = np.concatenate((active_s,active_s1))
        active_s1 = active_s + 50
        active_s = np.concatenate((active_s,active_s1))

        smneg_s = np.array([4,8,14,19,20])       # smear-neg, DS-TB compartments
        smneg_s1 = smneg_s + 25
        smneg_s = np.concatenate((smneg_s,smneg_s1))
        smneg_s1 = smneg_s + 50
        smneg_s = np.concatenate((smneg_s,smneg_s1))

        active_i = np.array([9,15])              # smear-pos, INH-m compartments
        active_i1 = active_i + 25
        active_i = np.concatenate((active_i,active_i1))
        active_i1 = active_i + 50
        active_i = np.concatenate((active_i,active_i1))

        smneg_i = np.array([5,10,16,21,22])       # smear-neg, INH-m compartments
        smneg_i1 = smneg_i + 25
        smneg_i = np.concatenate((smneg_i,smneg_i1))
        smneg_i1 = smneg_i + 50
        smneg_i = np.concatenate((smneg_i,smneg_i1))

        active_m = np.array([11,17])              # smear-pos, MDR-TB compartments
        active_m1 = active_m + 25
        active_m = np.concatenate((active_m,active_m1))
        active_m1 = active_m + 50
        active_m = np.concatenate((active_m,active_m1))

        smneg_m = np.array([6,12,18,23,24])       # smear-neg, MDR-TB compartments
        smneg_m1 = smneg_m + 25
        smneg_m = np.concatenate((smneg_m,smneg_m1))
        smneg_m1 = smneg_m + 50
        smneg_m = np.concatenate((smneg_m,smneg_m1))

        # forces of infection:
        force_s = beta_s*(X[active_s].sum() + relbeta_sn*X[smneg_s].sum())/pop1.sum()
        force_i = beta_s*fit_i*(X[active_i].sum() + relbeta_sn*X[smneg_i].sum())/pop1.sum()
        force_m = beta_s*fit_m*(X[active_m].sum() + relbeta_sn*X[smneg_m].sum())/pop1.sum()

        # 100x100 matrix of flow rates (cut-and-pasted from above):
        chg_arr = np.zeros((tb_num,tb_num))
        i1 = np.array([0,50]) # HIV-neg only
        i2 = np.array([25,75]) # HIV-pos only
        i3 = np.array([0,25,50,75]) # all
        i4 = np.array([0,25]) # new dx
        i5 = np.array([50,75]) # retx
        chg_arr[i1,i1+1] = force_s*(1-rapid)*X[i1] # infection that becomes latent
        chg_arr[i1,i1+2] = force_i*(1-rapid)*X[i1]
        chg_arr[i1,i1+3] = force_m*(1-rapid)*X[i1]
        chg_arr[i2,i2+1] = force_s*(1-rapid_h)*X[i2]
        chg_arr[i2,i2+2] = force_i*(1-rapid_h)*X[i2]
        chg_arr[i2,i2+3] = force_m*(1-rapid_h)*X[i2]
        chg_arr[i1,i1+4] = force_s*rapid*X[i1]     # primary progressive TB
        chg_arr[i1,i1+5] = force_i*rapid*X[i1]
        chg_arr[i1,i1+6] = force_m*rapid*X[i1]
        chg_arr[i2,i2+4] = force_s*rapid_h*X[i2]
        chg_arr[i2,i2+5] = force_i*rapid_h*X[i2]
        chg_arr[i2,i2+6] = force_m*rapid_h*X[i2]
        chg_arr[i1+1,i1+4] += react*X[i1+1]         # reactivation
        chg_arr[i1+2,i1+5] += react*X[i1+2]
        chg_arr[i1+3,i1+6] += react*X[i1+3]
        chg_arr[i2+1,i2+4] += react_h*X[i2+1]
        chg_arr[i2+2,i2+5] += react_h*X[i2+2]
        chg_arr[i2+3,i2+6] += react_h*X[i2+3]
        chg_arr[i1+1,i1+4] += force_s*rapid*(1-prot)*X[i1+1]     # reinfection to active
        chg_arr[i1+1,i1+5] += force_i*rapid*(1-prot)*X[i1+1]
        chg_arr[i1+1,i1+6] += force_m*rapid*(1-prot)*X[i1+1]
        chg_arr[i1+2,i1+4] += force_s*rapid*(1-prot)*X[i1+2]
        chg_arr[i1+2,i1+5] += force_i*rapid*(1-prot)*X[i1+2]
        chg_arr[i1+2,i1+6] += force_m*rapid*(1-prot)*X[i1+2]
        chg_arr[i1+3,i1+4] += force_s*rapid*(1-prot)*X[i1+3]
        chg_arr[i1+3,i1+5] += force_i*rapid*(1-prot)*X[i1+3]
        chg_arr[i1+3,i1+6] += force_m*rapid*(1-prot)*X[i1+3]
        chg_arr[i1+1,i1+2] += force_i*(1-rapid*(1-prot))*X[i1+1] # reinfection to latent
        chg_arr[i1+1,i1+3] += force_m*(1-rapid*(1-prot))*X[i1+1]
        chg_arr[i1+2,i1+1] += force_s*(1-rapid*(1-prot))*X[i1+2]
        chg_arr[i1+2,i1+3] += force_m*(1-rapid*(1-prot))*X[i1+2]
        chg_arr[i1+3,i1+1] += force_s*(1-rapid*(1-prot))*X[i1+3]
        chg_arr[i1+3,i1+2] += force_i*(1-rapid*(1-prot))*X[i1+3]
        chg_arr[i2+1,i2+4] += force_s*rapid_h*X[i2+1]     # reinfection to active (HIV)
        chg_arr[i2+1,i2+5] += force_i*rapid_h*X[i2+1]
        chg_arr[i2+1,i2+6] += force_m*rapid_h*X[i2+1]
        chg_arr[i2+2,i2+4] += force_s*rapid_h*X[i2+2]
        chg_arr[i2+2,i2+5] += force_i*rapid_h*X[i2+2]
        chg_arr[i2+2,i2+6] += force_m*rapid_h*X[i2+2]
        chg_arr[i2+3,i2+4] += force_s*rapid_h*X[i2+3]
        chg_arr[i2+3,i2+5] += force_i*rapid_h*X[i2+3]
        chg_arr[i2+3,i2+6] += force_m*rapid_h*X[i2+3]
        chg_arr[i1+4,i1+5] += force_i*(1-rapid_h)*X[i1+4]  # reinfection to latent (HIV)
        chg_arr[i1+4,i1+6] += force_m*(1-rapid_h)*X[i1+4]
        chg_arr[i1+5,i1+4] += force_s*(1-rapid_h)*X[i1+5]
        chg_arr[i1+5,i1+6] += force_m*(1-rapid_h)*X[i1+5]
        chg_arr[i1+6,i1+4] += force_s*(1-rapid_h)*X[i1+6]
        chg_arr[i1+6,i1+5] += force_i*(1-rapid_h)*X[i1+6]
        chg_arr[i1+4,i1+7] += (1/predx_dur-cure_sn)*prop_inf*X[i1+4]  # progression to dx-seeking
        chg_arr[i1+4,i1+8] += (1/predx_dur-cure_sn)*(1-prop_inf)*X[i1+4]
        chg_arr[i1+5,i1+9] += (1/predx_dur-cure_sn)*prop_inf*X[i1+5]
        chg_arr[i1+5,i1+10] += (1/predx_dur-cure_sn)*(1-prop_inf)*X[i1+5]
        chg_arr[i1+6,i1+11] += (1/predx_dur-cure_sn)*prop_inf*X[i1+6]
        chg_arr[i1+6,i1+12] += (1/predx_dur-cure_sn)*(1-prop_inf)*X[i1+6]
        chg_arr[i2+4,i2+7] += (1/predx_dur_h)*prop_inf_h*X[i2+4]  # progress to dx-seek, HIV
        chg_arr[i2+4,i2+8] += (1/predx_dur_h)*(1-prop_inf_h)*X[i2+4]
        chg_arr[i2+5,i2+9] += (1/predx_dur_h)*prop_inf_h*X[i2+5]
        chg_arr[i2+5,i2+10] += (1/predx_dur_h)*(1-prop_inf_h)*X[i2+5]
        chg_arr[i2+6,i2+11] += (1/predx_dur_h)*prop_inf_h*X[i2+6]
        chg_arr[i2+6,i2+12] += (1/predx_dur_h)*(1-prop_inf_h)*X[i2+6]

        chg_arr[7,13] += txrate(intopt,1,0,0,0)[0]*X[7] # progress to "dx in progress"
        chg_arr[8,14] += txrate(intopt,0,0,0,0)[0]*X[8] # new
        chg_arr[9,15] += txrate(intopt,1,1,0,0)[0]*X[9]
        chg_arr[10,16] += txrate(intopt,0,1,0,0)[0]*X[10]
        chg_arr[11,17] += txrate(intopt,1,2,0,0)[0]*X[11]
        chg_arr[12,18] += txrate(intopt,0,2,0,0)[0]*X[12]
        chg_arr[57,63] += txrate(intopt,1,0,1,0)[0]*X[57] # progress to "dx in progress"
        chg_arr[58,64] += txrate(intopt,0,0,1,0)[0]*X[58] # retreat
        chg_arr[59,65] += txrate(intopt,1,1,1,0)[0]*X[59]
        chg_arr[60,66] += txrate(intopt,0,1,1,0)[0]*X[60]
        chg_arr[61,67] += txrate(intopt,1,2,1,0)[0]*X[61]
        chg_arr[62,68] += txrate(intopt,0,2,1,0)[0]*X[62]
        chg_arr[7,19] += txrate(intopt,1,0,0,0)[2]*(1-acq_s)*X[7] # progress to "inapp tx"
        chg_arr[8,20] += txrate(intopt,0,0,0,0)[2]*(1-acq_s)*X[8] # new
        chg_arr[9,21] += txrate(intopt,1,1,0,0)[2]*(1-acq_i1)*X[9]
        chg_arr[10,22] += txrate(intopt,0,1,0,0)[2]*(1-acq_i1)*X[10]
        chg_arr[11,23] += txrate(intopt,1,2,0,0)[2]*X[11]
        chg_arr[12,24] += txrate(intopt,0,2,0,0)[2]*X[12]
        chg_arr[7,21] += txrate(intopt,1,0,0,0)[2]*(acq_s)*(1-acq_mdr)*X[7] # new resistance
        chg_arr[8,22] += txrate(intopt,0,0,0,0)[2]*(acq_s)*(1-acq_mdr)*X[8]
        chg_arr[7,23] += txrate(intopt,1,0,0,0)[2]*(acq_s)*(acq_mdr)*X[7] # new resistance
        chg_arr[8,24] += txrate(intopt,0,0,0,0)[2]*(acq_s)*(acq_mdr)*X[8]
        chg_arr[9,23] += txrate(intopt,1,1,0,0)[2]*(acq_i1)*X[9]
        chg_arr[10,24] += txrate(intopt,0,1,0,0)[2]*(acq_i1)*X[10]
        chg_arr[57,69] += txrate(intopt,1,0,1,0)[2]*(1-acq_s)*X[57] # progress to "inapp tx"
        chg_arr[58,70] += txrate(intopt,0,0,1,0)[2]*(1-acq_s)*X[58] # retreat
        chg_arr[59,71] += txrate(intopt,1,1,1,0)[2]*(1-acq_i2)*X[59]
        chg_arr[60,72] += txrate(intopt,0,1,1,0)[2]*(1-acq_i2)*X[60]
        chg_arr[61,73] += txrate(intopt,1,2,1,0)[2]*X[61]
        chg_arr[62,74] += txrate(intopt,0,2,1,0)[2]*X[62]
        chg_arr[57,71] += txrate(intopt,1,0,1,0)[2]*(acq_s)*(1-acq_mdr)*X[57] # new resistance
        chg_arr[58,72] += txrate(intopt,0,0,1,0)[2]*(acq_s)*(1-acq_mdr)*X[58]
        chg_arr[57,73] += txrate(intopt,1,0,1,0)[2]*(acq_s)*(acq_mdr)*X[57] # new resistance
        chg_arr[58,74] += txrate(intopt,0,0,1,0)[2]*(acq_s)*(acq_mdr)*X[58]
        chg_arr[59,73] += txrate(intopt,1,1,1,0)[2]*(acq_i2)*X[59]
        chg_arr[60,74] += txrate(intopt,0,1,1,0)[2]*(acq_i2)*X[60]
        chg_arr[13,51] += txrate(intopt,1,0,0,0)[1]*X[13] # progress to treated
        chg_arr[14,51] += txrate(intopt,0,0,0,0)[1]*X[14] # new
        chg_arr[15,52] += txrate(intopt,1,1,0,0)[1]*X[15]
        chg_arr[16,52] += txrate(intopt,0,1,0,0)[1]*X[16]
        chg_arr[17,53] += txrate(intopt,1,2,0,0)[1]*X[17]
        chg_arr[18,53] += txrate(intopt,0,2,0,0)[1]*X[18]
        chg_arr[63,51] += txrate(intopt,1,0,1,0)[1]*X[63] # progress to treated
        chg_arr[64,51] += txrate(intopt,0,0,1,0)[1]*X[64] # retreat
        chg_arr[65,52] += txrate(intopt,1,1,1,0)[1]*X[65]
        chg_arr[66,52] += txrate(intopt,0,1,1,0)[1]*X[66]
        chg_arr[67,53] += txrate(intopt,1,2,1,0)[1]*X[67]
        chg_arr[68,53] += txrate(intopt,0,2,1,0)[1]*X[68]

        chg_arr[7+25,13+25] += txrate(intopt,1,0,0,1)[0]*X[7+25] # progress to "dx in progress"
        chg_arr[8+25,14+25] += txrate(intopt,0,0,0,1)[0]*X[8+25] # new
        chg_arr[9+25,15+25] += txrate(intopt,1,1,0,1)[0]*X[9+25]
        chg_arr[10+25,16+25] += txrate(intopt,0,1,0,1)[0]*X[10+25]
        chg_arr[11+25,17+25] += txrate(intopt,1,2,0,1)[0]*X[11+25]
        chg_arr[12+25,18+25] += txrate(intopt,0,2,0,1)[0]*X[12+25]
        chg_arr[57+25,63+25] += txrate(intopt,1,0,1,1)[0]*X[57+25] # progress to "dx in progress"
        chg_arr[58+25,64+25] += txrate(intopt,0,0,1,1)[0]*X[58+25] # retreat
        chg_arr[59+25,65+25] += txrate(intopt,1,1,1,1)[0]*X[59+25]
        chg_arr[60+25,66+25] += txrate(intopt,0,1,1,1)[0]*X[60+25]
        chg_arr[61+25,67+25] += txrate(intopt,1,2,1,1)[0]*X[61+25]
        chg_arr[62+25,68+25] += txrate(intopt,0,2,1,1)[0]*X[62+25]
        chg_arr[7+25,19+25] += txrate(intopt,1,0,0,1)[2]*(1-acq_s)*X[7+25] # progress to "inapp tx"
        chg_arr[8+25,20+25] += txrate(intopt,0,0,0,1)[2]*(1-acq_s)*X[8+25] # new
        chg_arr[9+25,21+25] += txrate(intopt,1,1,0,1)[2]*(1-acq_i1)*X[9+25]
        chg_arr[10+25,22+25] += txrate(intopt,0,1,0,1)[2]*(1-acq_i1)*X[10+25]
        chg_arr[11+25,23+25] += txrate(intopt,1,2,0,1)[2]*X[11+25]
        chg_arr[12+25,24+25] += txrate(intopt,0,2,0,1)[2]*X[12+25]
        chg_arr[7+25,21+25] += txrate(intopt,1,0,0,1)[2]*(acq_s)*(1-acq_mdr)*X[7+25] # new resistance
        chg_arr[8+25,22+25] += txrate(intopt,0,0,0,1)[2]*(acq_s)*(1-acq_mdr)*X[8+25]
        chg_arr[7+25,23+25] += txrate(intopt,1,0,0,1)[2]*(acq_s)*(acq_mdr)*X[7+25] # new resistance
        chg_arr[8+25,24+25] += txrate(intopt,0,0,0,1)[2]*(acq_s)*(acq_mdr)*X[8+25]
        chg_arr[9+25,23+25] += txrate(intopt,1,1,0,1)[2]*(acq_i1)*X[9+25]
        chg_arr[10+25,24+25] += txrate(intopt,0,1,0,1)[2]*(acq_i1)*X[10+25]
        chg_arr[57+25,69+25] += txrate(intopt,1,0,1,1)[2]*(1-acq_s)*X[57+25] # progress to "inapp tx"
        chg_arr[58+25,70+25] += txrate(intopt,0,0,1,1)[2]*(1-acq_s)*X[58+25] # retreat
        chg_arr[59+25,71+25] += txrate(intopt,1,1,1,1)[2]*(1-acq_i2)*X[59+25]
        chg_arr[60+25,72+25] += txrate(intopt,0,1,1,1)[2]*(1-acq_i2)*X[60+25]
        chg_arr[61+25,73+25] += txrate(intopt,1,2,1,1)[2]*X[61+25]
        chg_arr[62+25,74+25] += txrate(intopt,0,2,1,1)[2]*X[62+25]
        chg_arr[57+25,71+25] += txrate(intopt,1,0,1,1)[2]*(acq_s)*(1-acq_mdr)*X[57+25] # new resistance
        chg_arr[58+25,72+25] += txrate(intopt,0,0,1,1)[2]*(acq_s)*(1-acq_mdr)*X[58+25]
        chg_arr[57+25,73+25] += txrate(intopt,1,0,1,1)[2]*(acq_s)*(acq_mdr)*X[57+25] # new resistance
        chg_arr[58+25,74+25] += txrate(intopt,0,0,1,1)[2]*(acq_s)*(acq_mdr)*X[58+25]
        chg_arr[59+25,73+25] += txrate(intopt,1,1,1,1)[2]*(acq_i2)*X[59+25]
        chg_arr[60+25,74+25] += txrate(intopt,0,1,1,1)[2]*(acq_i2)*X[60+25]
        chg_arr[13+25,51+25] += txrate(intopt,1,0,0,1)[1]*X[13+25] # progress to treated
        chg_arr[14+25,51+25] += txrate(intopt,0,0,0,1)[1]*X[14+25] # new
        chg_arr[15+25,52+25] += txrate(intopt,1,1,0,1)[1]*X[15+25]
        chg_arr[16+25,52+25] += txrate(intopt,0,1,0,1)[1]*X[16+25]
        chg_arr[17+25,53+25] += txrate(intopt,1,2,0,1)[1]*X[17+25]
        chg_arr[18+25,53+25] += txrate(intopt,0,2,0,1)[1]*X[18+25]
        chg_arr[63+25,51+25] += txrate(intopt,1,0,1,1)[1]*X[63+25] # progress to treated
        chg_arr[64+25,51+25] += txrate(intopt,0,0,1,1)[1]*X[64+25] # retreat
        chg_arr[65+25,52+25] += txrate(intopt,1,1,1,1)[1]*X[65+25]
        chg_arr[66+25,52+25] += txrate(intopt,0,1,1,1)[1]*X[66+25]
        chg_arr[67+25,53+25] += txrate(intopt,1,2,1,1)[1]*X[67+25]
        chg_arr[68+25,53+25] += txrate(intopt,0,2,1,1)[1]*X[68+25]

        chg_arr[i4+19,i4+57] += (12./dur_fail)*defvfail_s*X[i4+19] # default from "inapp tx" (new cases)
        chg_arr[i4+20,i4+58] += (12./dur_fail)*defvfail_s*X[i4+20] # assume zero cost b/c 1/2 of default is cured
        chg_arr[i4+21,i4+59] += (12./dur_fail)*defvfail_i*X[i4+21]
        chg_arr[i4+22,i4+60] += (12./dur_fail)*defvfail_i*X[i4+22]
        chg_arr[i4+23,i4+61] += (12./dur_fail)*defvfail_m*X[i4+23]
        chg_arr[i4+24,i4+62] += (12./dur_fail)*defvfail_m*X[i4+24]
        chg_arr[i5+19,i5+7] += (12./dur_fail)*defvfail_s*X[i5+19] # default from "inapp tx" (retx cases)
        chg_arr[i5+20,i5+8] += (12./dur_fail)*defvfail_s*X[i5+20] # assume zero cost b/c 1/2 of default is cured
        chg_arr[i5+21,i5+9] += (12./dur_fail)*defvfail_i*X[i5+21]
        chg_arr[i5+22,i5+10] += (12./dur_fail)*defvfail_i*X[i5+22]
        chg_arr[i5+23,i5+11] += (12./dur_fail)*defvfail_m*X[i5+23]
        chg_arr[i5+24,i5+12] += (12./dur_fail)*defvfail_m*X[i5+24]
        chg_arr[i4+19,i4+51] += (12./dur_fail)*(1-defvfail_s)*X[i4+19]*(1-fail_rt) # failure from "inapp tx" (new TB cases)
        chg_arr[i4+20,i4+51] += (12./dur_fail)*(1-defvfail_s)*X[i4+20]*(1-fail_rt) # these are the pts that are then successfully treated
        chg_arr[i4+21,i4+52] += (12./dur_fail)*(1-defvfail_i)*X[i4+21]*(1-fail_i2) # the ones below re-fail
        chg_arr[i4+22,i4+52] += (12./dur_fail)*(1-defvfail_i)*X[i4+22]*(1-fail_i2)
        chg_arr[i4+23,i4+53] += (12./dur_fail)*(1-defvfail_m)*X[i4+23]*(1-fail_m2)
        chg_arr[i4+24,i4+53] += (12./dur_fail)*(1-defvfail_m)*X[i4+24]*(1-fail_m2)
        chg_arr[i4+19,i4+69] += (12./dur_fail)*(1-defvfail_s)*X[i4+19]*(fail_rt) # re-failure from "inapp tx" -> now classified as retx
        chg_arr[i4+20,i4+70] += (12./dur_fail)*(1-defvfail_s)*X[i4+20]*(fail_rt)
        chg_arr[i4+21,i4+71] += (12./dur_fail)*(1-defvfail_i)*X[i4+21]*(fail_i2)
        chg_arr[i4+22,i4+72] += (12./dur_fail)*(1-defvfail_i)*X[i4+22]*(fail_i2)
        chg_arr[i4+23,i4+73] += (12./dur_fail)*(1-defvfail_m)*X[i4+23]*(fail_m2)
        chg_arr[i4+24,i4+74] += (12./dur_fail)*(1-defvfail_m)*X[i4+24]*(fail_m2)
        chg_arr[i5+19,i5+1] += (12./dur_fail)*(1-defvfail_s)*X[i5+19]*(1-fail_rt) # failure from "inapp tx"
        chg_arr[i5+20,i5+1] += (12./dur_fail)*(1-defvfail_s)*X[i5+20]*(1-fail_rt) # retx pts, so "re-fail" stays in the same box
        chg_arr[i5+21,i5+2] += (12./dur_fail)*(1-defvfail_i)*X[i5+21]*(1-fail_i2)
        chg_arr[i5+22,i5+2] += (12./dur_fail)*(1-defvfail_i)*X[i5+22]*(1-fail_i2)
        chg_arr[i5+23,i5+3] += (12./dur_fail)*(1-defvfail_m)*X[i5+23]*(1-fail_m2)
        chg_arr[i5+24,i5+3] += (12./dur_fail)*(1-defvfail_m)*X[i5+24]*(1-fail_m2)
        i6 = np.arange(25)                          # HIV infection
        chg_arr[i6,i6+25] += hiv_inc*X[i6]
        chg_arr[i6+50,i6+75] += hiv_inc*X[i6+50]
        i7 = np.arange(4)                           # inappropriate tx for TB
        chg_arr[i7,i7+50] += prop_cough*(1-invspec(intopt,0,0))*X[i7]
        chg_arr[i7+25,i7+75] += prop_cough*(1-invspec(intopt,1,0))*X[i7+25]
        i8 = np.arange(50)                         # mortality
        chg_arr[i6,0] += mortality(0,0)*X[i6]       # background mortality
        chg_arr[i6+50,0] += mortality(0,0)*X[i6+50]
        chg_arr[i6+25,0] += mortality(0,1)*X[i6+25] # HIV deaths
        chg_arr[i6+75,0] += mortality(0,1)*X[i6+75]
        i9 = np.array([7,9,11,13,15,17])  # smear-pos
        i10 = np.array([4,5,6,8,10,12,14,16,18,19,20,21,22,23,24])  # smear-neg
        chg_arr[i9,0] += mortality(1,0)*X[i9]              # smear-positive TB deaths
        chg_arr[i9+25,0] += mortality(1,1)*X[i9+25]
        chg_arr[i9+50,0] += mortality(1,0)*X[i9+50]
        chg_arr[i9+75,0] += mortality(1,1)*X[i9+75]
        chg_arr[i10,0] += mortality(2,0)*X[i10]              # smear-negative TB deaths
        chg_arr[i10+25,0] += mortality(2,1)*X[i10+25]
        chg_arr[i10+50,0] += mortality(2,0)*X[i10+50]
        chg_arr[i10+75,0] += mortality(2,1)*X[i10+75]
        chg_arr[i1+7,i1+1] += cure_sp*X[i1+7]              # smear-positive self-cure (no HIV)
        chg_arr[i1+9,i1+2] += cure_sp*X[i1+9]
        chg_arr[i1+11,i1+3] += cure_sp*X[i1+11]
        chg_arr[i1+13,i1+1] += cure_sp*X[i1+13]
        chg_arr[i1+15,i1+2] += cure_sp*X[i1+15]
        chg_arr[i1+17,i1+3] += cure_sp*X[i1+17]
        chg_arr[i1+8,i1+1] += cure_sn*X[i1+8]              # smear-negative self-cure (no HIV)
        chg_arr[i1+10,i1+2] += cure_sn*X[i1+10]
        chg_arr[i1+12,i1+3] += cure_sn*X[i1+12]
        chg_arr[i1+14,i1+1] += cure_sn*X[i1+14]
        chg_arr[i1+16,i1+2] += cure_sn*X[i1+16]
        chg_arr[i1+18,i1+3] += cure_sn*X[i1+18]
        chg_arr[i1+4,i1+1] += cure_sn*X[i1+4]
        chg_arr[i1+5,i1+2] += cure_sn*X[i1+5]
        chg_arr[i1+6,i1+3] += cure_sn*X[i1+6]
        chg_arr[i1+19,i1+1] += cure_sn*X[i1+19]
        chg_arr[i1+21,i1+2] += cure_sn*X[i1+21]
        chg_arr[i1+23,i1+3] += cure_sn*X[i1+23]
        chg_arr[i1+20,i1+1] += cure_sn*X[i1+20]
        chg_arr[i1+22,i1+2] += cure_sn*X[i1+22]
        chg_arr[i1+24,i1+3] += cure_sn*X[i1+24]

        # 0 = incidence of new (non-retreatment) TB:
        incnew = np.sum(chg_arr[i4,i4+4]) + np.sum(chg_arr[i4,i4+5]) + np.sum(chg_arr[i4,i4+6]) + \
                 np.sum(chg_arr[i4+1,i4+4]) + np.sum(chg_arr[i4+1,i4+5]) + np.sum(chg_arr[i4+1,i4+6]) + \
                 np.sum(chg_arr[i4+2,i4+4]) + np.sum(chg_arr[i4+2,i4+5]) + np.sum(chg_arr[i4+2,i4+6]) + \
                 np.sum(chg_arr[i4+3,i4+4]) + np.sum(chg_arr[i4+3,i4+5]) + np.sum(chg_arr[i4+3,i4+6])
        # 1 = incidence of retreatment TB:
        incretx = np.sum(chg_arr[i5,i5+4]) + np.sum(chg_arr[i5,i5+5]) + np.sum(chg_arr[i5,i5+6]) + \
                 np.sum(chg_arr[i5+1,i5+4]) + np.sum(chg_arr[i5+1,i5+5]) + np.sum(chg_arr[i5+1,i5+6]) + \
                 np.sum(chg_arr[i5+2,i5+4]) + np.sum(chg_arr[i5+2,i5+5]) + np.sum(chg_arr[i5+2,i5+6]) + \
                 np.sum(chg_arr[i5+3,i5+4]) + np.sum(chg_arr[i5+3,i5+5]) + np.sum(chg_arr[i5+3,i5+6]) + \
                 np.sum(chg_arr[i4+19,i4+57]) + np.sum(chg_arr[i4+20,i4+58]) + np.sum(chg_arr[i4+21,i4+59]) + \
                 np.sum(chg_arr[i4+22,i4+60]) + np.sum(chg_arr[i4+23,i4+61]) + np.sum(chg_arr[i4+24,i4+62]) + \
                 np.sum(chg_arr[i5+19,i5+7]) + np.sum(chg_arr[i5+20,i5+8]) + np.sum(chg_arr[i5+21,i5+9]) + \
                 np.sum(chg_arr[i5+22,i5+10]) + np.sum(chg_arr[i5+23,i5+11]) + np.sum(chg_arr[i5+24,i5+12]) + \
                 np.sum((12./dur_fail)*(1-defvfail_s)*X[i3+19]) + np.sum((12./dur_fail)*(1-defvfail_s)*X[i3+20]) + \
                 np.sum((12./dur_fail)*(1-defvfail_i)*X[i3+21]) + np.sum((12./dur_fail)*(1-defvfail_i)*X[i3+22]) + \
                 np.sum((12./dur_fail)*(1-defvfail_m)*X[i3+23]) + np.sum((12./dur_fail)*(1-defvfail_m)*X[i3+24])
        # 2 = incidence of new, INH-resistant TB:
        incinhnew = np.sum(chg_arr[i4,i4+5]) + np.sum(chg_arr[i4+1,i4+5])+ + np.sum(chg_arr[i4+2,i4+5]) + \
                    np.sum(chg_arr[i4+3,i4+5])
        # 3 = incidence of retreatment, INH-resistant TB:
        incinhretx = np.sum(chg_arr[i5,i5+5]) + np.sum(chg_arr[i5+1,i5+5]) + np.sum(chg_arr[i5+2,i5+5]) + \
                     np.sum(chg_arr[i5+3,i5+5]) + np.sum(chg_arr[i4+21,i4+59]) + \
                     np.sum(chg_arr[i4+22,i4+60]) + np.sum(chg_arr[i5+21,i5+9]) + np.sum(chg_arr[i5+22,i5+10]) + \
                     np.sum((12./dur_fail)*(1-defvfail_i)*X[i3+21]) + np.sum((12./dur_fail)*(1-defvfail_i)*X[i3+22])
        # 4 = incidence of new, MDR-TB:
        incmdrnew = np.sum(chg_arr[i4,i4+6]) + np.sum(chg_arr[i4+1,i4+6])+ + np.sum(chg_arr[i4+2,i4+6]) + \
                    np.sum(chg_arr[i4+3,i4+6])
        # 5 = incidence of retreatment, MDR-TB:
        incmdrretx = np.sum(chg_arr[i5,i5+6]) + np.sum(chg_arr[i5+1,i5+6]) + np.sum(chg_arr[i5+2,i5+6]) + \
                     np.sum(chg_arr[i5+3,i5+6]) + np.sum(chg_arr[i4+23,i4+61]) + \
                     np.sum(chg_arr[i4+24,i4+62]) + np.sum(chg_arr[i5+23,i5+11]) + np.sum(chg_arr[i5+24,i5+12]) + \
                     np.sum((12./dur_fail)*(1-defvfail_m)*X[i3+23]) + np.sum((12./dur_fail)*(1-defvfail_m)*X[i3+24])
        # 6 = incidence of HIV-associated TB:
        inctbhiv = np.sum(chg_arr[i2,i2+4]) + np.sum(chg_arr[i2,i2+5]) + np.sum(chg_arr[i2,i2+6]) + \
          np.sum(chg_arr[i2+1,i2+4]) + np.sum(chg_arr[i2+1,i2+5]) + np.sum(chg_arr[i2+1,i2+6]) +\
           np.sum(chg_arr[i2+2,i2+4]) + np.sum(chg_arr[i2+2,i2+5]) + np.sum(chg_arr[i2+2,i2+6]) +\
            np.sum(chg_arr[i2+3,i2+4]) + np.sum(chg_arr[i2+3,i2+5]) + np.sum(chg_arr[i2+3,i2+6]) +\
             np.sum(chg_arr[25+19,25+57]) + np.sum(chg_arr[25+20,25+58]) + np.sum(chg_arr[25+21,25+59]) +\
              np.sum(chg_arr[25+22,25+60]) + np.sum(chg_arr[25+23,25+61]) + np.sum(chg_arr[25+24,25+62]) +\
               np.sum(chg_arr[75+19,75+7]) + np.sum(chg_arr[75+20,75+8]) + np.sum(chg_arr[75+21,75+9]) + np.sum(chg_arr[75+22,75+10]) +\
                np.sum(chg_arr[75+23,75+11]) + np.sum(chg_arr[75+24,75+12])  +\
                 np.sum((12./dur_fail)*(1-defvfail_s)*X[i2+19]) + np.sum((12./dur_fail)*(1-defvfail_s)*X[i2+20]) +\
                 np.sum((12./dur_fail)*(1-defvfail_i)*X[i2+21]) + np.sum((12./dur_fail)*(1-defvfail_i)*X[i2+22]) +\
                      np.sum((12./dur_fail)*(1-defvfail_m)*X[i2+23]) + np.sum((12./dur_fail)*(1-defvfail_m)*X[i2+24])
        # 7 = TB mortality
        TBmort = np.sum(mortality(1,0)*X[i9]) + np.sum(mortality(1,1)*X[i9+25]) + np.sum(mortality(1,0)*X[i9+50]) + \
                 np.sum(mortality(1,1)*X[i9+75]) + np.sum(mortality(2,0)*X[i10]) + np.sum(mortality(2,1)*X[i10+25]) + \
                 np.sum(mortality(2,0)*X[i10+50]) + np.sum(mortality(2,1)*X[i10+75])

        i11 = np.arange(4,25)
        i11b = i11 + 25
        i11 = np.concatenate((i11,i11b))
        i11b = i11 + 50
        i11 = np.concatenate((i11,i11b))
        # 8 = TB prevalence
        # sum(X[i11]) includes all active cases
        # the next indent is people on 1st-line therapy (for 6 months)
        # the next indent is people on cat2 regimen (for 8 months)
        # and the final indent is people on MDR/2nd-line therapy (for 20 months)
        TBprev = np.sum(X[i11])+ 0.5*(np.sum(txrate(intopt,1,0,0,0)[1]*X[13]) + np.sum(txrate(intopt,0,0,0,0)[1]*X[14]) + \
                                   np.sum(txrate(intopt,1,1,0,0)[1]*X[15]) + np.sum(txrate(intopt,0,1,0,0)[1]*X[16]) + \
                                   np.sum(txrate(intopt,1,2,0,0)[1]*X[17]) + np.sum(txrate(intopt,0,2,0,0)[1]*X[18]) + \
                                   np.sum(txrate(intopt,1,0,0,1)[1]*X[25+13]) + np.sum(txrate(intopt,0,0,0,1)[1]*X[25+14]) + \
                                   np.sum(txrate(intopt,1,1,0,1)[1]*X[25+15]) + np.sum(txrate(intopt,0,1,0,1)[1]*X[25+16]) + \
                                   np.sum(txrate(intopt,1,2,0,1)[1]*X[25+17]) + np.sum(txrate(intopt,0,2,0,1)[1]*X[25+18])) + \
                                   (8./12.)*(np.sum(txrate(intopt,1,0,1,0)[1]*X[63]) + np.sum(txrate(intopt,0,0,1,0)[1]*X[64]) + \
                                             np.sum(txrate(intopt,1,1,1,0)[1]*X[65]) + np.sum(txrate(intopt,0,1,1,0)[1]*X[66]) + \
                                             np.sum(txrate(intopt,1,0,1,1)[1]*X[25+63]) + np.sum(txrate(intopt,0,0,1,1)[1]*X[25+64]) + \
                                             np.sum(txrate(intopt,1,1,1,1)[1]*X[25+65]) + np.sum(txrate(intopt,0,1,1,1)[1]*X[25+66]) + \
                                             np.sum((12./dur_fail)*0.5*X[i4+19]*(1-fail_rt)) + \
                                             np.sum((12./dur_fail)*0.5*X[i4+20]*(1-fail_rt)) + \
                                             np.sum((12./dur_fail)*0.5*X[i4+21]*(1-fail_i2)) + \
                                             np.sum((12./dur_fail)*0.5*X[i4+22]*(1-fail_i2)) + \
                                             np.sum((12./dur_fail)*0.5*X[i5+19]*(1-fail_rt)) + \
                                             np.sum((12./dur_fail)*0.5*X[i5+20]*(1-fail_rt)) + \
                                             np.sum((12./dur_fail)*0.5*X[i5+20]*(1-fail_rt)) + \
                                             np.sum((12./dur_fail)*0.5*X[i5+21]*(1-fail_i2)) + \
                                             np.sum((12./dur_fail)*0.5*X[i5+22]*(1-fail_i2))) + \
                                             (20./12.)*(np.sum(txrate(intopt,1,2,1,0)[1]*X[67]) + \
                                                        np.sum(txrate(intopt,0,2,1,0)[1]*X[68]) + \
                                                        np.sum(txrate(intopt,1,2,1,1)[1]*X[25+67]) + \
                                                        np.sum(txrate(intopt,0,2,1,1)[1]*X[25+68]) + \
                                                        np.sum((12./dur_fail)*0.5*X[i4+23]*(1-fail_m2)) + \
                                                        np.sum((12./dur_fail)*0.5*X[i4+24]*(1-fail_m2)) + \
                                                        np.sum((12./dur_fail)*0.5*X[i5+23]*(1-fail_m2)) + \
                                                        np.sum((12./dur_fail)*0.5*X[i5+24]*(1-fail_m2)))
        i12 = np.concatenate((np.arange(25,50), np.arange(75,100)))
        # 9 = HIV prevalence
        HIVprev = np.sum(X[i12])

        # 10 = Cost
        # First indent is diagnostic costs for people without TB.
        # Second indent is diagnostic costs.
        # Third indent is costs of inappropriate treatment (e.g., 1st-line for MDR-TB) after diagnosis
        # Fourth indent is costs of appropriate treatment
        # Fifth indent is inappropriate treatment that is made empirically
        # Sixth indent is failure: inappropriate treatment followed by appropriate treatment
        # Seventh through tenth indents are indents 3-6 but extended to people with HIV
        i13a = np.arange(4)
        i13b = i13a + 25
        i13 = np.concatenate((i13a,i13b),axis=1)
        i14a = np.arange(50,54)
        i14b = i14a + 25
        i14 = np.concatenate((i14a,i14b),axis=1)
        cost = np.sum(prop_cough*X[i13a]*dxcost(intopt,0,0))+np.sum(prop_cough*X[i14a]*dxcost(intopt,1,0)) + \
               np.sum(prop_cough*X[i13b]*dxcost(intopt,0,1))+np.sum(prop_cough*X[i14b]*dxcost(intopt,1,1)) + \
               dx_rate*(np.sum(txrate(intopt,1,0,0,0)[3]*X[7]) + np.sum(txrate(intopt,0,0,0,0)[3]*X[8]) + \
                        np.sum(txrate(intopt,1,1,0,0)[3]*X[9]) + np.sum(txrate(intopt,0,1,0,0)[3]*X[10]) + \
                        np.sum(txrate(intopt,1,2,0,0)[3]*X[11]) + np.sum(txrate(intopt,0,2,0,0)[3]*X[12]) + \
                        np.sum(txrate(intopt,1,0,1,0)[3]*X[57]) + np.sum(txrate(intopt,0,0,1,0)[3]*X[58]) + \
                        np.sum(txrate(intopt,1,1,1,0)[3]*X[59]) + np.sum(txrate(intopt,0,1,1,0)[3]*X[60]) + \
                        np.sum(txrate(intopt,1,2,1,0)[3]*X[61]) + np.sum(txrate(intopt,0,2,1,0)[3]*X[62]) + \
                        np.sum(txrate(intopt,1,0,0,1)[3]*X[25+7]) + np.sum(txrate(intopt,0,0,0,1)[3]*X[25+8]) + \
                        np.sum(txrate(intopt,1,1,0,1)[3]*X[25+9]) + np.sum(txrate(intopt,0,1,0,1)[3]*X[25+10]) + \
                        np.sum(txrate(intopt,1,2,0,1)[3]*X[25+11]) + np.sum(txrate(intopt,0,2,0,1)[3]*X[25+12]) + \
                        np.sum(txrate(intopt,1,0,1,1)[3]*X[25+57]) + np.sum(txrate(intopt,0,0,1,1)[3]*X[25+58]) + \
                        np.sum(txrate(intopt,1,1,1,1)[3]*X[25+59]) + np.sum(txrate(intopt,0,1,1,1)[3]*X[25+60]) + \
                        np.sum(txrate(intopt,1,2,1,1)[3]*X[25+61]) + np.sum(txrate(intopt,0,2,1,1)[3]*X[25+62])) + \
                    np.sum(txrate(intopt,1,0,0,0)[2]*txrate(intopt,1,0,0,0)[5]*X[7]) + \
                    np.sum(txrate(intopt,0,0,0,0)[2]*txrate(intopt,0,0,0,0)[5]*X[8]) + \
                    np.sum(txrate(intopt,1,1,0,0)[2]*txrate(intopt,1,1,0,0)[5]*X[9]) + \
                    np.sum(txrate(intopt,0,1,0,0)[2]*txrate(intopt,0,1,0,0)[5]*X[10]) + \
                    np.sum(txrate(intopt,1,2,0,0)[2]*txrate(intopt,1,2,0,0)[5]*X[11]) + \
                    np.sum(txrate(intopt,0,2,0,0)[2]*txrate(intopt,0,2,0,0)[5]*X[12]) + \
                    np.sum(txrate(intopt,1,0,1,0)[2]*txrate(intopt,1,0,1,0)[5]*X[57]) + \
                    np.sum(txrate(intopt,0,0,1,0)[2]*txrate(intopt,0,0,1,0)[5]*X[58]) + \
                    np.sum(txrate(intopt,1,1,1,0)[2]*txrate(intopt,1,1,1,0)[5]*X[59]) + \
                    np.sum(txrate(intopt,0,1,1,0)[2]*txrate(intopt,0,1,1,0)[5]*X[60]) + \
                    np.sum(txrate(intopt,1,2,1,0)[2]*txrate(intopt,1,2,1,0)[5]*X[61]) + \
                    np.sum(txrate(intopt,0,2,1,0)[2]*txrate(intopt,0,2,1,0)[5]*X[62]) + \
                        np.sum(txrate(intopt,1,0,0,0)[0]*txrate(intopt,1,0,0,0)[4]*X[7]) + \
                        np.sum(txrate(intopt,0,0,0,0)[0]*txrate(intopt,0,0,0,0)[4]*X[8]) + \
                        np.sum(txrate(intopt,1,1,0,0)[0]*txrate(intopt,1,1,0,0)[4]*X[9]) + \
                        np.sum(txrate(intopt,0,1,0,0)[0]*txrate(intopt,0,1,0,0)[4]*X[10]) + \
                        np.sum(txrate(intopt,1,2,0,0)[0]*txrate(intopt,1,2,0,0)[4]*X[11]) + \
                        np.sum(txrate(intopt,0,2,0,0)[0]*txrate(intopt,0,2,0,0)[4]*X[12]) + \
                        np.sum(txrate(intopt,1,0,1,0)[0]*txrate(intopt,1,0,1,0)[4]*X[57]) + \
                        np.sum(txrate(intopt,0,0,1,0)[0]*txrate(intopt,0,0,1,0)[4]*X[58]) + \
                        np.sum(txrate(intopt,1,1,1,0)[0]*txrate(intopt,1,1,1,0)[4]*X[59]) + \
                        np.sum(txrate(intopt,0,1,1,0)[0]*txrate(intopt,0,1,1,0)[4]*X[60]) + \
                        np.sum(txrate(intopt,1,2,1,0)[0]*txrate(intopt,1,2,1,0)[4]*X[61]) + \
                        np.sum(txrate(intopt,0,2,1,0)[0]*txrate(intopt,0,2,1,0)[4]*X[62]) - \
                    np.sum(chg_arr[19,57]*txrate(intopt,1,0,0,0)[5]) - \
                    np.sum(chg_arr[20,58]*txrate(intopt,0,0,0,0)[5]) - \
                    np.sum(chg_arr[21,59]*txrate(intopt,1,1,0,0)[5]) - \
                    np.sum(chg_arr[22,60]*txrate(intopt,0,1,0,0)[5]) - \
                    np.sum(chg_arr[23,61]*txrate(intopt,1,2,0,0)[5]) - \
                    np.sum(chg_arr[24,62]*txrate(intopt,0,2,0,0)[5]) - \
                    np.sum(chg_arr[69,57]*txrate(intopt,1,0,1,0)[5]) - \
                    np.sum(chg_arr[70,58]*txrate(intopt,0,0,1,0)[5]) - \
                    np.sum(chg_arr[71,59]*txrate(intopt,1,1,1,0)[5]) - \
                    np.sum(chg_arr[72,60]*txrate(intopt,0,1,1,0)[5]) - \
                    np.sum(chg_arr[73,61]*txrate(intopt,1,2,1,0)[5]) - \
                    np.sum(chg_arr[74,62]*txrate(intopt,0,2,1,0)[5]) + \
                        np.sum((txrate(intopt,1,0,1,0)[4]+txrate(intopt,1,2,1,0)[3])*(12./dur_fail)*(1-defvfail_s)*X[i1+19]) + \
                        np.sum((txrate(intopt,0,0,1,0)[4]+txrate(intopt,1,2,1,0)[3])*(12./dur_fail)*(1-defvfail_s)*X[i1+20]) + \
                        np.sum((txrate(intopt,1,1,1,0)[4]+txrate(intopt,1,2,1,0)[3])*(12./dur_fail)*(1-defvfail_s)*X[i1+21]) + \
                        np.sum((txrate(intopt,0,1,1,0)[4]+txrate(intopt,1,2,1,0)[3])*(12./dur_fail)*(1-defvfail_s)*X[i1+22]) + \
                        np.sum((txrate(intopt,1,2,1,0)[4]+txrate(intopt,1,2,1,0)[3])*(12./dur_fail)*(1-defvfail_s)*X[i1+23]) + \
                        np.sum((txrate(intopt,0,2,1,0)[4]+txrate(intopt,1,2,1,0)[3])*(12./dur_fail)*(1-defvfail_s)*X[i1+24]) + \
                    np.sum(txrate(intopt,1,0,0,1)[2]*txrate(intopt,1,0,0,1)[5]*X[7+25]) + \
                    np.sum(txrate(intopt,0,0,0,1)[2]*txrate(intopt,0,0,0,1)[5]*X[8+25]) + \
                    np.sum(txrate(intopt,1,1,0,1)[2]*txrate(intopt,1,1,0,1)[5]*X[9+25]) + \
                    np.sum(txrate(intopt,0,1,0,1)[2]*txrate(intopt,0,1,0,1)[5]*X[10+25]) + \
                    np.sum(txrate(intopt,1,2,0,1)[2]*txrate(intopt,1,2,0,1)[5]*X[11+25]) + \
                    np.sum(txrate(intopt,0,2,0,1)[2]*txrate(intopt,0,2,0,1)[5]*X[12+25]) + \
                    np.sum(txrate(intopt,1,0,1,1)[2]*txrate(intopt,1,0,1,1)[5]*X[57+25]) + \
                    np.sum(txrate(intopt,0,0,1,1)[2]*txrate(intopt,0,0,1,1)[5]*X[58+25]) + \
                    np.sum(txrate(intopt,1,1,1,1)[2]*txrate(intopt,1,1,1,1)[5]*X[59+25]) + \
                    np.sum(txrate(intopt,0,1,1,1)[2]*txrate(intopt,0,1,1,1)[5]*X[60+25]) + \
                    np.sum(txrate(intopt,1,2,1,1)[2]*txrate(intopt,1,2,1,1)[5]*X[61+25]) + \
                    np.sum(txrate(intopt,0,2,1,1)[2]*txrate(intopt,0,2,1,1)[5]*X[62+25]) + \
                        np.sum(txrate(intopt,1,0,0,1)[0]*txrate(intopt,1,0,0,1)[4]*X[7+25]) + \
                        np.sum(txrate(intopt,0,0,0,1)[0]*txrate(intopt,0,0,0,1)[4]*X[8+25]) + \
                        np.sum(txrate(intopt,1,1,0,1)[0]*txrate(intopt,1,1,0,1)[4]*X[9+25]) + \
                        np.sum(txrate(intopt,0,1,0,1)[0]*txrate(intopt,0,1,0,1)[4]*X[10+25]) + \
                        np.sum(txrate(intopt,1,2,0,1)[0]*txrate(intopt,1,2,0,1)[4]*X[11+25]) + \
                        np.sum(txrate(intopt,0,2,0,1)[0]*txrate(intopt,0,2,0,1)[4]*X[12+25]) + \
                        np.sum(txrate(intopt,1,0,1,1)[0]*txrate(intopt,1,0,1,1)[4]*X[57+25]) + \
                        np.sum(txrate(intopt,0,0,1,1)[0]*txrate(intopt,0,0,1,1)[4]*X[58+25]) + \
                        np.sum(txrate(intopt,1,1,1,1)[0]*txrate(intopt,1,1,1,1)[4]*X[59+25]) + \
                        np.sum(txrate(intopt,0,1,1,1)[0]*txrate(intopt,0,1,1,1)[4]*X[60+25]) + \
                        np.sum(txrate(intopt,1,2,1,1)[0]*txrate(intopt,1,2,1,1)[4]*X[61+25]) + \
                        np.sum(txrate(intopt,0,2,1,1)[0]*txrate(intopt,0,2,1,1)[4]*X[62+25]) - \
                    np.sum(chg_arr[19+25,57+25]*txrate(intopt,1,0,0,1)[5]) - \
                    np.sum(chg_arr[20+25,58+25]*txrate(intopt,0,0,0,1)[5]) - \
                    np.sum(chg_arr[21+25,59+25]*txrate(intopt,1,1,0,1)[5]) - \
                    np.sum(chg_arr[22+25,60+25]*txrate(intopt,0,1,0,1)[5]) - \
                    np.sum(chg_arr[23+25,61+25]*txrate(intopt,1,2,0,1)[5]) - \
                    np.sum(chg_arr[24+25,62+25]*txrate(intopt,0,2,0,1)[5]) - \
                    np.sum(chg_arr[69+25,57+25]*txrate(intopt,1,0,1,1)[5]) - \
                    np.sum(chg_arr[70+25,58+25]*txrate(intopt,0,0,1,1)[5]) - \
                    np.sum(chg_arr[71+25,59+25]*txrate(intopt,1,1,1,1)[5]) - \
                    np.sum(chg_arr[72+25,60+25]*txrate(intopt,0,1,1,1)[5]) - \
                    np.sum(chg_arr[73+25,61+25]*txrate(intopt,1,2,1,1)[5]) - \
                    np.sum(chg_arr[74+25,62+25]*txrate(intopt,0,2,1,1)[5]) + \
                        np.sum((txrate(intopt,1,0,1,1)[4]+txrate(intopt,1,2,1,1)[3])*(12./dur_fail)*(1-defvfail_s)*X[i2+19]) + \
                        np.sum((txrate(intopt,0,0,1,1)[4]+txrate(intopt,1,2,1,1)[3])*(12./dur_fail)*(1-defvfail_s)*X[i2+20]) + \
                        np.sum((txrate(intopt,1,1,1,1)[4]+txrate(intopt,1,2,1,1)[3])*(12./dur_fail)*(1-defvfail_s)*X[i2+21]) + \
                        np.sum((txrate(intopt,0,1,1,1)[4]+txrate(intopt,1,2,1,1)[3])*(12./dur_fail)*(1-defvfail_s)*X[i2+22]) + \
                        np.sum((txrate(intopt,1,2,1,1)[4]+txrate(intopt,1,2,1,1)[3])*(12./dur_fail)*(1-defvfail_s)*X[i2+23]) + \
                        np.sum((txrate(intopt,0,2,1,1)[4]+txrate(intopt,1,2,1,1)[3])*(12./dur_fail)*(1-defvfail_s)*X[i2+24])
        return incnew, incretx, incinhnew, incinhretx, incmdrnew, incmdrretx, inctbhiv, TBmort, TBprev, HIVprev, cost

###########################################################
# 9. LOW INCIDENCE AND EMERGING MDR SCENARIO CALCULATOR
###########################################################

# One challenge is that equilibrium scenarios result in about 80% of all incident TB being due to recent infection.
# Low incidence settings have more TB due to reactivation and thus cannot be modeled at equilibrium,
#   as the proportion of incident TB reflecting reactivation is a key driver of impact.
# This section allows the model to generate low-incidence scenarios that appropriately have a lower proportion of
#   TB due to recent infection, with a higher proportion due to reactivation.

# If the target incidence is <50 per 100,000, the program generates an equilibrium 50 years in the past
#    with an incidence of 50 per 100,000/yr.
# It then reduces the transmission parameter by a set proportion at that time, searching until it finds the
#    reduction that allows it to re-create the target incidence as specified by the user.
# This allows the program to generate a situation in which the majority of incident TB is due to
#    reactivation rather than recent infection.

    equil_round = 0.
    if g.target_inc2 < 50:
        g.beta_s *= 1. + (g.target_inc2 - 50.)/82.
        while equil_round<50:
            pre_time = 50.
            time_burnin = np.arange(0,pre_time+1.,1.)
            EQUI = odeint(diffeq,equipop,time_burnin)
            if abs(incprevmort(EQUI[50,:])[0] + incprevmort(EQUI[50,:])[1] -g.target_inc2)<0.05:
                break
            print "calculating low-incidence scenario..."
            hb_log("calculating low-incidence scneario...")
            print equil_round, g.beta_s, incprevmort(EQUI[50,:])[0] + incprevmort(EQUI[50,:])[1]
            hb_log("{} {} {}".format(equil_round, g.beta_s, incprevmort(EQUI[50,:])[0] + incprevmort(EQUI[50,:])[1]))
            g.beta_s += (g.target_inc2 - incprevmort(EQUI[50,:])[0] - incprevmort(EQUI[50,:])[1])/(0.25*g.target_inc2)
            equil_round +=1
    if g.target_inc2<50:
        equipop = EQUI[50,:]


# Below also allows for an emerging MDR scenario, in which the final MDR-TB prevalence is
#    higher at the end of 5 years than at the beginning.

    equil_round = 0.
    if g.target_mdr2 > g.target_mdr:
        g.fit_m *= 1. + ((g.target_mdr2 - g.target_mdr)/0.25)*(1./g.fit_m)
        while equil_round<50:
            pre_time2 = 5.
            time_burnin2 = np.arange(0,pre_time2+0.1,0.1)
            EQUI2 = odeint(diffeq,equipop,time_burnin2)
            if abs(incprevmort(EQUI2[50,:])[4]/incprevmort(EQUI2[50,:])[0]-g.target_mdr2)<0.0005:
                break
            print "calculating emerging MDR scenario..."
            hb_log("calculating emerging MDR scenario...")
            print equil_round, g.fit_m, incprevmort(EQUI2[50,:])[4]/incprevmort(EQUI2[50,:])[0]
            hb_log("{} {} {}".format(equil_round, g.fit_m, incprevmort(EQUI2[50,:])[4]/incprevmort(EQUI2[50,:])[0]))
            g.fit_m += (g.target_mdr2 - incprevmort(EQUI2[50,:])[4]/incprevmort(EQUI2[50,:])[0])/(5.*g.target_mdr2)
            equil_round +=1


###########################################################
# 10. SOLVE DIFFERENTIAL EQUATIONS AND REPORT RESULTS
###########################################################

    duration = 5.                        # 5 year time horizon
    timestep = 0.01                      # calculate equations every 0.01 yrs
    time_range = np.arange(0, duration+timestep, timestep) # vector for calculations

    OUT = odeint(diffeq,equipop,time_range)  # solve the ODE's
    incvect = incprevmort(OUT[500,:])    # calculate the incidence/prevalence/mortality/cost at the end of year 5

    scenario_name = ['Baseline (Smear)','Xpert for smear-positive','Xpert for HIV+',
                     'Xpert for previously treated','Xpert for sm-neg HIV+ or prev tx',
                     'Xpert for all HIV+ or prev tx','Xpert for smear-negative',
                     'Xpert for all',g.ud_strat_name]

    """
    print " "
    print "BASELINE"
    print "Incidence, new:", np.around(incvect[0],decimals=1), "per 100,000"
    print "Incidence, retx:", np.around(incvect[1],decimals=1), "per 100,000"
    print "Incidence, total:", np.around(incvect[0] + incvect[1],decimals=1), "per 100,000"
    print "Incidence, INH new:", np.around(incvect[2]/incvect[0]*100,decimals=1), "%"
    print "Incidence, INH retx:", np.around(incvect[3]/incvect[1]*100,decimals=1), "%"
    print "Incidence, MDR new:", np.around(incvect[4]/incvect[0]*100,decimals=2), "%"
    print "Incidence, MDR retx:", np.around(incvect[5]/incvect[1]*100,decimals=2), "%"
    print "Incidence, MDR total:", np.around(incvect[4]+incvect[5],decimals=2), "per 100,000"
    print "Incidence, TB/HIV:", np.around(incvect[6]/(incvect[0]+incvect[1])*100,decimals=1), "%"
    print "TB mortality:", np.around(incvect[7],decimals=1), "per 100,000"
    print "TB duration:", np.around(incvect[8]/(incvect[0]+incvect[1]),decimals=2), "years"
    print "HIV prevalence:", np.around(incvect[9]/1000,decimals=1), "%"
    """
    costvect = np.zeros((100))
    for x in range(100):
        costvect[x] = incprevmort(OUT[1+x,:])[10]
    """
    print "Cost in Year 1: $", np.around(np.sum(costvect[:])/100.)
    print "Cost in Year 5: $", np.around(incvect[10])
    print " "

    outdata['0'] = {'inc_n':np.around(incvect[0],decimals=1),'inc_r':np.around(incvect[1],decimals=1),
                     'inc_t':np.around(incvect[0] + incvect[1],decimals=1),'inc_inh_n':np.around(incvect[2]/incvect[0]*100,decimals=1),
                     'inc_inh_r':np.around(incvect[3]/incvect[1]*100,decimals=1),'inc_mdr_n':np.around(incvect[4]/incvect[0]*100,decimals=2),
                     'inc_mdr_r':np.around(incvect[5]/incvect[1]*100,decimals=2), 'inc_mdr_t':np.around(incvect[4]+incvect[5],decimals=2),
                     'tb_mort':np.around(incvect[7],decimals=1),'tb_dur':np.around(incvect[8]/(incvect[0]+incvect[1]),decimals=2),
                     'hiv_prev': np.around(incvect[9]/1000,decimals=1), 'inc_tb_hiv': np.around(incvect[6]/(incvect[0]+incvect[1])*100,decimals=1),
                     'cost1':np.around(np.sum(costvect[:])/100.),'cost5':np.around(incvect[10]),
                     'name':scenario_name[0]}

    outdata['progress'] += 1

    json_write(tmp_filename, outdata)
    """

    if g.int_select==9:
        for abc in range(9):
            g.intopt = abc
            OUT = odeint(diffeq,equipop,time_range)  # solve the ODE's from the same equilibrium population
            #print "INTERVENTION", intopt
            hb_log('Running Intervention {}'.format(abc))
            incvect2 = incprevmort(OUT[500,:])
            """
            print "Incidence, new:", np.around(incvect2[0],decimals=1), "per 100,000"
            print "Incidence, retx:", np.around(incvect2[1],decimals=1), "per 100,000"
            print "Incidence, total:", np.around(incvect2[0] + incvect2[1],decimals=1), "per 100,000"
            print np.around((1-(incvect2[0] + incvect2[1])/(incvect[0] + incvect[1]))*100,decimals=1),"% reduction"
            print "Incidence, INH new:", np.around(incvect2[2]/incvect2[0]*100,decimals=1), "%"
            print "Incidence, INH retx:", np.around(incvect2[3]/incvect2[1]*100,decimals=1), "%"
            print "Incidence, MDR new:", np.around(incvect2[4]/incvect2[0]*100,decimals=2), "%"
            print "Incidence, MDR retx:", np.around(incvect2[5]/incvect2[1]*100,decimals=2), "%"
            print "Incidence, MDR total:", np.around(incvect2[4]+incvect2[5],decimals=2), "per 100,000"
            print np.around((1-(incvect2[4] + incvect2[5])/(incvect[4] + incvect[5]))*100,decimals=1),"% reduction"
            print "Incidence, TB/HIV:", np.around(incvect2[6]/(incvect2[0]+incvect2[1])*100,decimals=1), "%"
            print "TB mortality:", np.around(incvect2[7],decimals=1), "per 100,000"
            print np.around((1-(incvect2[7])/(incvect[7]))*100,decimals=1),"% reduction"
            print "TB duration:", np.around(incvect2[8]/(incvect2[0]+incvect2[1]),decimals=2), "years"
            print "HIV prevalence:", np.around(incvect2[9]/1000,decimals=1), "%"
            """
            # now also calculate costs in year 1 (from timestep 0 to timestep 99)
            costvect2 = np.zeros((100))
            for x in range(100):
                costvect2[x] = incprevmort(OUT[1+x,:])[10]
            """
            print "Cost in Year 1: $", np.around(np.sum(costvect2[:])/100.)
            print np.around(((np.sum(costvect2[:])/100)-(np.sum(costvect[:])/100.))/(np.sum(costvect[:])/100.)*100,decimals=1),"% increase"
            print "Cost at End of Year 5: $", np.around(incvect2[10])
            if incvect2[10]>incvect[10]:
                print np.around(((incvect2[10])/(incvect[10])-1)*100,decimals=1),"% increase"
            if incvect2[10]<=incvect[10]:
                print np.around((1-(incvect2[10])/(incvect[10]))*100,decimals=1),"% decrease"
            print " "
            """
            outdata[str(abc)] = {'inc_n':np.around(incvect2[0],decimals=1),'inc_r':np.around(incvect2[1],decimals=1),
                                    'inc_t':np.around(incvect2[0] + incvect2[1],decimals=1),'inc_inh_n':np.around(incvect2[2]/incvect2[0]*100,decimals=1),
                                    'inc_inh_r':np.around(incvect2[3]/incvect2[1]*100,decimals=1),'inc_mdr_n':np.around(incvect2[4]/incvect2[0]*100,decimals=2),
                                    'inc_mdr_r':np.around(incvect2[5]/incvect2[1]*100,decimals=2), 'inc_mdr_t':np.around(incvect2[4]+incvect2[5],decimals=2),
                                    'tb_mort':np.around(incvect2[7],decimals=1),'tb_dur':np.around(incvect2[8]/(incvect2[0]+incvect2[1]),decimals=2),
                                    'hiv_prev': np.around(incvect2[9]/1000,decimals=1), 'inc_tb_hiv': np.around(incvect2[6]/(incvect2[0]+incvect2[1])*100,decimals=1),
                                    'cost1':np.around(np.sum(costvect2[:])/100.),'cost5':np.around(incvect2[10]),
                                    'name':scenario_name[abc]}
            outdata['progress'] += 1
            json_write(tmp_filename, outdata)

    hb_log("------- Execution Ends -------")

if __name__=="__main__":

    hb_log("Run from command line")

    if len(sys.argv) != 2:
        hb_log( "Command line arguments invalid" )
        exit(0)

    tmp_filename = sys.argv[1]

    hb_log( "Loading data from: {}".format(tmp_filename) )
    with open(tmp_filename, 'r') as fd:
        jdata = json.load(fd)
    
    run ( jdata )
