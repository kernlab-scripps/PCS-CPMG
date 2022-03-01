#==============================================================================
# XPLOR -py Script for Finding the Structure of ADK using Expectation Max.
#==============================================================================
# The following script is a modified version of the standard refine.py script
# included in the XPLOR-NIH NMR tutorials.
# The goal is to use a set of PCS to refine ADK structure using rigid body
# restraints.
# Unlike standard PCS refinement, this script was built for PCS values
# extracted from CPMG experiments where there are total of four possible
# values per residue. Only one of these PCS values are correct.
# The script uses an psuedo-bayes expectation maximization algorith to
# determine the correct PCS values for each residues during the refinement.


#==============================================================================
# User Inputs:
#==============================================================================
# PCS File containing Four PCS values per residue w/ Uncertainity.
# Initial Structure used as the starting point for refinement.
ADK="4AKE"
Metals = ['Co']
PCS_Files = ['{}_Modeled_FourPCS_{}.txt'.format(ADK, M) for M in Metals]
Xax_Array = [16.8] # 1e-32 m^3 
Xrh_Array = [0.65] # 1e-32 m^3 
Starting_Structure = 'Prior_Best.pdb'

#==============================================================================
# Global Inputs
#==============================================================================
Metal_ID = 302
Global_Unc = {'Co': 0.01}
Prior_Prob = 0.25
Starting_ForceC = 0.0005
TensorValues = {'Co': [0, 0, 0]}


#==============================================================================
# Start Standard XPLOR Run
#==============================================================================
# Load standard XPLOR modules and prepare starting structure for refinement
xplor.requireVersion("2.34")                 # Check Version
(opts,args) = xplor.parseArguments(["slow"]) # check for command-line typos
quick=False                                  # Run Full Refinement
import protocol                              # General module for XPLOR

outFilename = "SCRIPT_STRUCTURE.sa"          # Output name
numberOfStructures=1                         # 1 Structure Per run

# Issue a totally random seed number
# Standard init random seed uses the time, which can be very similar
# if multiple runs are queued simultaneously.
import numpy as np, numpy.random
protocol.initRandomSeed(int(np.random.random()*10**7))

# Set up Protein Parameters:
command = xplor.command
protocol.initParams("protein")

# Next Steps involve loading structure and fixing geometries
protocol.loadPDB(Starting_Structure, deleteUnknownAtoms=True)
protocol.fixupCovalentGeom(maxIters=100,useVDW=1)

# This will calculate the total change between end and start structures.
from posDiffPotTools import create_PosDiffPot
refRMSD = create_PosDiffPot("refRMSD","name CA or name C or name N",
                            pdbFile=Starting_Structure,
                            cmpSel="not name H*")

#==============================================================================
# Preparing PCS Data to be used for refinement:
#==============================================================================
# During the refinement, we will be determining the probability of each
# PCS for every residue. Ideally, we'd like to determine that one of these
# is correct and other three (per residue) are incorrect.
# Once that is determined, we'd like to have XPLOR only fit the correct
# PCS and to ignore the other three.
# The easiest way to do this is have the force constant (K) go to zero
# for the incorrect PCS values.

# Make Dictionary to Store PCS values and Probabilities:
PCS_Dic, Prob_Dic = {}, {}
for M in Metals:
    PCS_Dic[M] = {}
    Prob_Dic[M] = {}

# Load PCS and Uncertainity from PCS File:
import math
for M, M_File in zip(Metals, PCS_Files):
    with open(M_File, 'rb') as fh: 
        for line in fh.readlines(): 
            s = line.split('\t')
            PCS_Dic[M].setdefault(s[0], []).append(float(s[1]))
    
for M in Metals:
    fname = 'Prior_File_{}.txt'.format(M)
    with open(fname, 'rb') as fh:
        for line in fh.readlines():
            s = line.split('\t')
            Prob_Dic[M][s[0]] = [float(i) for i in s[1:]]

# Create Objects to Store Results during Run:
Q_Sig_Array = []
Running_Probs, Running_PCS = {}, {}
for M in Metals: 
    Running_Probs[M] = {}
    Running_PCS[M] = {}
    for Res in PCS_Dic[M]:
        Running_PCS[M][Res] = []
        Running_Probs[M][Res] = []


#==============================================================================
# Create PCS Library
#==============================================================================
# Directory with each PCS (4*n total) are stored as individual Files:
import os
if not os.path.exists('PCS_tbl_Files'):
    os.mkdir('PCS_tbl_Files')

# Write individual PCS and Uncertainity Files:
for M in Metals:
    Global_Error = Global_Unc[M]/math.sqrt(0.25)
    for Res in PCS_Dic[M]: 
        for ent, PCS in enumerate(PCS_Dic[M][Res]):
            outstring = ('assign ( resid 500 and name 00 )\n'
                         '       ( resid 500 and name Z )\n'
                         '       ( resid 500 and name X )\n'
                         '       ( resid 500 and name Y )\n'
                         '       ( resid {0} and name ZN )\n'    
                         '       ( resid {1} and name HN )   {2}  {3}\n'
                         .format(Metal_ID, Res, PCS, Global_Error))
            with open('PCS_tbl_Files/Res{0}_{1}_{2}.tbl'.format(Res, ent, M),
                      'w') as out:
                out.write(outstring)

        
#==============================================================================
# Maximum Likelihood Elements
#==============================================================================
def Prob_with_Class_Av(Prob_Dic):
    # Determine Class Averages: 
    Class_Averages = {}
    for M in Prob_Dic:
        Prob_Array = []
        for Res in Prob_Dic[M]:
            Prob_Array.append(Prob_Dic[M][Res])
        Class_Averages[M] = np.mean(np.array(Prob_Array), axis=0)
    # Adjust Probability Dic with Class Averages and Normalize:
    for M in Prob_Dic:
        for Res in Prob_Dic[M]:
            Prob_Dic[M][Res] = Prob_Dic[M][Res]*Class_Averages[M]
            Prob_Dic[M][Res] = Prob_Dic[M][Res]/np.sum(Prob_Dic[M][Res])
    return Prob_Dic


def Calc_Global_Sigma(Square_Dic, Prob_Dic):
    # Determine Global Sigma Value:
    Sigma_Dic = {}
    for M in Square_Dic: 
        Sigma = []
        for Res in Square_Dic[M]:
            Sigma.append(np.sum(np.array(Square_Dic[M][Res])*Prob_Dic[M][Res]))
        Sigma_Dic[M] = math.sqrt(np.sum(Sigma)/(len(Sigma)-1.0))
    return Sigma_Dic


def Calc_Diamag_Signs(Square_Dic, Sigma): 
    Res_Likelihoods = {}
    Prior = 0.25
    for M in Metals:
        s1 = 1.0/(Sigma[M]*np.sqrt(2*np.pi))
        s2 = -1.0/(2.0*np.square(Sigma[M]))
        for Res in Square_Dic[M]:
            Square_Array = np.array(Square_Dic[M][Res])
            Gauss = s1*np.exp(s2*Square_Array)*Prior  
            Res_Likelihoods.setdefault(Res, []).append(Gauss)
    Sign_Dic = {}
    for Res in Res_Likelihoods:
        if len(Res_Likelihoods[Res]) == 1:
            Sign_Dic[Res] = np.array([0.5, 0.5, 0.5, 0.5])
        else: 
            Class_LL = np.log(np.array(Res_Likelihoods[Res]) + 1e-50).T
            Pos_Sign_Like = np.exp(np.sum(Class_LL[0]) + np.sum(Class_LL[1]))
            Neg_Sign_Like = np.exp(np.sum(Class_LL[2]) + np.sum(Class_LL[3]))
            Norm_Pos = Pos_Sign_Like/(Pos_Sign_Like + Neg_Sign_Like)
            Norm_Neg = Neg_Sign_Like/(Pos_Sign_Like + Neg_Sign_Like)
            Sign_Dic[Res] = np.array([Norm_Pos, Norm_Pos, Norm_Neg, Norm_Neg])
        if math.isnan(Sign_Dic[Res][0]) == True:
            print(Res, Res_Likelihoods[Res], Class_LL, Pos_Sign_Like, Norm_Pos,
                  Norm_Neg, Sign_Dic[Res])
    return Sign_Dic
            


def Adjust_Prob(Square_Dic, Sigma):
    # Make New Prob Dictionary
    New_Prob_Dic = {}
    for M in Square_Dic:
        New_Prob_Dic[M] = {}
        s1 = 1.0/(Sigma[M]*np.sqrt(2*np.pi))
        s2 = -1.0/(2.0*np.square(Sigma[M]))
        for Res in Square_Dic[M]:
            Gauss = s1*np.exp(s2*np.array(Square_Dic[M][Res]))
            New_Prob_Dic[M][Res] = Gauss / np.sum(Gauss)
    return New_Prob_Dic


def Calc_Total_Like(Square_Dic, Sigma, Prob_Dic):
    # Calculate Total Probability Weighted Log Likelihood - Q
    Q_All = 0
    # Calculate Log Likelihood on Individaul Residue and add to Q
    for M in Square_Dic: 
        Q_Vals = []
        s1 = -0.5*np.log(2.0*np.pi)
        s2 = -1.0*np.log(Sigma[M])
        s3 = -0.5/np.square(Sigma[M])
        for Res in Square_Dic[M]:
            Squared_Array = np.array(Square_Dic[M][Res])
            Q_Vals.append(np.sum((s1+s2+s3*Squared_Array) * Prob_Dic[M][Res]))
        print(M, np.sum(Q_Vals))
        Q_All += np.sum(Q_Vals)
    return Q_All


def Dirichlet(Prob_Dic):
    # To make the run less deterministic, the probabilities and be resampled
    # at each run. This is performed by sampling from the dirichlet distribution.
    # Where the Dirichlet pdf is the conjugate prior of a multinomial.
    # First calculate the four Class Averages and their Variances
    # Make an array of all Probilities (normalized likelihoods)
    Prob_Array = []
    for M in Metals:
        for Res in Prob_Dic[M]:
            Prob_Array.append(Prob_Dic[M][Res])
    # Calculate the average probability for each of the four positions(columns)
    Class_Vars = np.square(np.std(np.array(Prob_Array), axis = 0))
    # Calculate the standard variance of the mean.
    V = np.sum(Class_Vars / len(Prob_Array))
    # Resample Dictionary:
    for M in Metals: 
        for Res in Prob_Dic[M]:
            a1 = np.square(np.array(Prob_Dic[M][Res]))/V
            Prob_Dic[M][Res] = np.random.dirichlet(a1)
    return Prob_Dic


def Report_PCS_Calc(potList, Running_PCS):
    # Dictionary to Store Calculated Values
    PCS_Calc_Dic = {}
    for M in Running_PCS: 
       PCS_Calc_Dic[M] = {}
    # Extract Calculated Values for each Residue
    for M in Running_PCS:
        for Res in potList['PCSs_{}'.format(M)].energyReports():
            PCS_Calc = list(filter(None, potList['PCSs_{}'.format(M)][Res[0]].showViolations().
                                                                 split(' ')))[-4]
            Res_Number = Res[0].split('_')[2]
            M = Res[0].split('_')[0]
            PCS_Calc_Dic[M].setdefault(Res_Number, []).append(PCS_Calc)
    # Store current PCS calculated values.
    for M in Running_PCS:
        for Res in PCS_Calc_Dic[M]:
            Running_PCS[M][Res].append(PCS_Calc_Dic[M][Res][0])


def Output_Script(Sigma, Ending_Q, Prob_Dic, Q_Sig_Array, Running_Probs, 
                  Running_PCS, potList):
    Output = [Ending_Q]
    for M in Metals:
        Output.append(Sigma[M])
    Q_Sig_Array.append(Output)
    # Store PCS Values
    Report_PCS_Calc(potList, Running_PCS)
    # Store Current Probabilities
    for M in Metals: 
        for Res in Prob_Dic[M]:
            Running_Probs[M][Res].append(Prob_Dic[M][Res])

#==============================================================================
# Recalculate PCS Scales via EM Method
#==============================================================================
def ReCalculate_Scales(potList, PCS_Dic, Prob_Dic, Iteration,
                       Q_Sig_Array, Running_Probs, Running_PCS):
    # Orient Tensor and Determine PCS Values
    calcTensorOrientation(Co_Tensor['Co'])    
    Square_Dic = {}
    for M in PCS_Dic: 
        Square_Dic[M] = {}
        for i in potList['PCSs_{}'.format(M)].energyReports():
            s = i[0].split('_')
            M, Res = s[0], s[2]
            Squared_Difference = i[1]/potList['PCSs_{}'.format(M)][i[0]].scale()
            Square_Dic[M].setdefault(Res, []).append(Squared_Difference)
    if Square_Dic['Co'] == {}:
        return potList, Prob_Dic, Iteration, TensorValues
    
    # Conditional Expectation Maximization Until Convergence: 
    Likelihoods = [-1000, -999]
    while abs(Likelihoods[-1]-Likelihoods[-2]) > 0.001:
        Prob_Dic = Prob_with_Class_Av(Prob_Dic)
        Sigma = Calc_Global_Sigma(Square_Dic, Prob_Dic)
        Diamag_Dic = Calc_Diamag_Signs(Square_Dic, Sigma)
        Prob_Dic = Adjust_Prob(Square_Dic, Sigma)
        Likelihoods.append(Calc_Total_Like(Square_Dic, Sigma, Prob_Dic))
        if len(Likelihoods) == 100:
            with open('NoCovergence.txt', 'a+') as out:
                out.write('Did not converge - Iteration {}'.format(Iteration))
            break

    # Resampel with Dirchlet
    Prob_Dic = Dirichlet(Prob_Dic)

    # Write Output files
    Output_Script(Sigma, Likelihoods[-1], Prob_Dic, Q_Sig_Array, 
                  Running_Probs, Running_PCS, potList)
    # At the start of the simulating annealing, the force constants are kept
    # lower to allow the protein to sample a large conformational space
    # without much energetic restriction.
    # As cooling occurs, the force constant is ramped up to enforce the
    # PCS restraint.
    #Scale_Range = np.exp(np.linspace(0,np.log(2.0),860)) * np.linspace(0.04,1,860)
    Scale_Range = np.logspace(np.log10(0.001), np.log10(5.0), 1720)
    Scale_Range = np.concatenate((Scale_Range, np.array([5.0]*500)))
    Scale_Value = Scale_Range[Iteration]
    Iteration += 1
    # Update the potenial list with new Probabilities:
    for M in Square_Dic: 
        for Res in Prob_Dic[M]:
            for ent, Prob in enumerate(Prob_Dic[M][Res]):
                if Prob < 1e-100:
                    Scale = Scale_Value/((Global_Unc[M]/np.sqrt(1e-100))**2)
                else:
                    Scale = Scale_Value/((Global_Unc[M]/np.sqrt(Prob))**2)
                Restraint = '{}_NH_{}_{}'.format(M, Res, ent)
                potList['PCSs_{}'.format(M)][Restraint].setScale(Scale)
    return potList, Prob_Dic, Iteration, TensorValues

#==============================================================================
# Set up PCS Potential and Tensor
#==============================================================================
# Make Potentials List
from potList import PotList
potList = PotList()

# Setup ramped and high temperature specific parameters
from simulationTools import MultRamp, StaticRamp, InitialParams
rampedParams = []
highTempParams = []

# Create Paramagnetic Alignment Tensor
from varTensorTools import create_VarTensor
Co_Tensor = {}
for (M, Xax, Xrh) in zip(Metals, Xax_Array, Xrh_Array):
    if M == 'Co':
        oTensor = create_VarTensor(M)
        oTensor.setDaMax(6000)
        oTensor.setDa(Xax*10000/(12.0*3.14))
        oTensor.setRh(abs(Xrh/Xax))
        Co_Tensor[M] = oTensor    
 
    
# Add Individual PCS to Potential List - Called RDC because the RDC function
# is being used. They share the same mathemetical formulas.  
from rdcPotTools import create_RDCPot
PCS_Co_Pot = PotList('PCSs_Co')
for M in Metals:
    for Res in PCS_Dic[M]:
        for ent, PCS in enumerate(PCS_Dic[M][Res]):
            Restraint = '{0}_NH_{1}_{2}'.format(M, Res, ent)
            PCS_File = 'PCS_tbl_Files/Res{0}_{1}_{2}.tbl'.format(Res, ent, M)
            
            pcs_restr = create_RDCPot(Restraint, PCS_File, Co_Tensor[M])
            pcs_restr.setScale(Starting_ForceC/
                               (Global_Unc[M]/math.sqrt(0.25))**2)
            pcs_restr.setUseDistance(True)
            pcs_restr.setShowAllRestraints(1)
            pcs_restr.setThreshold(0.02)
            PCS_Co_Pot.append(pcs_restr)
potList.append(PCS_Co_Pot)

# At each step along the cooling process, probabilities will be recalculated.
Iteration = 0
Statement = ("potList, Prob_Dic, Iteration, TensorValues = ReCalculate_Scales(potList, "
             "PCS_Dic, Prob_Dic, Iteration, Q_Sig_Array, "
             "Running_Probs, Running_PCS)")
rampedParams.append(StaticRamp(Statement))

# Recalculate Euler Angles at each stage of cooling: 
from varTensorTools import calcTensorOrientation

# Fix Axial and Rhombic Components
for M in Co_Tensor.values(): 
    M.setFreedom("fixDa, fixRh") #fix tensor Rh, Da, vary orientation


#==============================================================================
# Set up the Remaining XPLOR Potentials
#==============================================================================
# hbdb - knowledge-based backbone hydrogen bond term
from xplorPot import XplorPot
protocol.initHBDB()
potList.append( XplorPot('HBDB') )

#New torsion angle database potential
from torsionDBPotTools import create_TorsionDBPot
torsionDB = create_TorsionDBPot('torsionDB')
potList.append( torsionDB )
rampedParams.append( MultRamp(.0002,2,"torsionDB.setScale(VALUE)") )

# setup parameters for atom-atom repulsive term. (van der Waals-like term)
from repelPotTools import create_RepelPot,initRepel
repel = create_RepelPot('repel')
potList.append(repel)
rampedParams.append( StaticRamp("initRepel(repel,use14=False)") )
rampedParams.append( MultRamp(.0004,4,  "repel.setScale( VALUE)") )
# nonbonded interaction only between CA atoms
highTempParams.append( StaticRamp("""initRepel(repel,
                                               use14=True,
                                               scale=0.0004,
                                               repel=1.2,
                                               moveTol=45,
                                               interactingAtoms='name CA'
                                               )""") )

# Selected 1-4 interactions.
import torsionDBPotTools
repel14 = torsionDBPotTools.create_Terminal14Pot('repel14')
potList.append(repel14)
highTempParams.append(StaticRamp("repel14.setScale(0)"))
rampedParams.append(MultRamp(0.0004, 4, "repel14.setScale(VALUE)"))
#
#
potList.append( XplorPot("BOND") )
potList.append( XplorPot("ANGL") )
potList['ANGL'].setThreshold( 5 )
rampedParams.append( MultRamp(0.04,1,"potList['ANGL'].setScale(VALUE)") )
potList.append( XplorPot("IMPR") )
potList['IMPR'].setThreshold( 5 )
rampedParams.append( MultRamp(0.01,1,"potList['IMPR'].setScale(VALUE)") )

# Give atoms uniform weights, except for the anisotropy axis
protocol.massSetup()


#==============================================================================
# IVM setup
#==============================================================================
from ivm import IVM
# At High Temperatures, use rigid body constraints
dyn_HT = IVM()
dyn_HT.group( AtomSel("resid 1:28 or (resid 62:113) or (resid 181:217)"))
dyn_HT.group( AtomSel("resid 32:55"))
dyn_HT.group( AtomSel("resid 128:157 or resid 302"))
dyn_HT.group( AtomSel("resid 116:125"))
dyn_HT.group( AtomSel("resid 165:178"))
# Use torsion movements:
protocol.torsionTopology(dyn_HT)

# During Cooling, use rigid body constraints:
dyn = IVM()
dyn.group( AtomSel("resid 1:28 or (resid 62:113) or (resid 181:217)"))
dyn.group( AtomSel("resid 32:55"))
dyn.group( AtomSel("resid 128:157 or resid 302"))
dyn.group( AtomSel("resid 116:125"))
dyn.group( AtomSel("resid 165:178"))
# Use torsion movements:
protocol.torsionTopology(dyn)

# Refine rigid bodys in cartesian space:
minc = IVM()
protocol.initMinimize(minc)
minc.group( AtomSel("resid 1:28 or (resid 62:113) or (resid 181:217)"))
minc.group( AtomSel("resid 32:55"))
minc.group( AtomSel("resid 128:157 or resid 302"))
minc.group( AtomSel("resid 116:125"))
minc.group( AtomSel("resid 165:178"))
# Use cartesian movements:
protocol.cartesianTopology(minc)

# Refine secondary structure bundles in torsion space:
dyn_cool = IVM()
protocol.initMinimize(dyn_cool)
dyn_cool.group( AtomSel("resid 1:28 or (resid 81:113) or (resid 181:198)"))
dyn_cool.group( AtomSel("resid 201:217"))
dyn_cool.group( AtomSel("resid 32:40"))
dyn_cool.group( AtomSel("resid 44:55"))
dyn_cool.group( AtomSel("resid 62:79"))
dyn_cool.group( AtomSel("resid 128:157 or resid 302"))
dyn_cool.group( AtomSel("resid 116:125"))
dyn_cool.group( AtomSel("resid 165:178"))
# Use torsion movements:
protocol.torsionTopology(dyn_cool)

# Refine secondary structure bundles in carteissan space:
minc_cool = IVM()
protocol.initMinimize(minc_cool)
minc_cool.group( AtomSel("resid 1:28 or (resid 81:113) or (resid 181:198)"))
minc_cool.group( AtomSel("resid 201:217"))
minc_cool.group( AtomSel("resid 32:40"))
minc_cool.group( AtomSel("resid 44:55"))
minc_cool.group( AtomSel("resid 62:79"))
minc_cool.group( AtomSel("resid 128:157 or resid 302"))
minc_cool.group( AtomSel("resid 116:125"))
minc_cool.group( AtomSel("resid 165:178"))
# Use cartessian movements:
protocol.cartesianTopology(minc_cool)

#==============================================================================
# Setup Simulated Annealing Temperatures
#==============================================================================
from simulationTools import AnnealIVM
init_t  = 2000
cool = AnnealIVM(initTemp =init_t,
                 finalTemp=25,
                 tempStep=(init_t-25.0)/1720.0,
                 ivm=dyn,
                 rampedParams = rampedParams)
cool2 = AnnealIVM(initTemp=25,
                 finalTemp=20,
                 tempStep=0.1,
                 ivm=dyn_cool,
                 rampedParams = rampedParams)

#==============================================================================
# Accept Structure if Ending Potential is good enough
#==============================================================================
def accept(potList):
    #return True if current structure meets acceptance criteria
    if potList['PCSs_Co'].rms()>1.2: #this might be tightened some
        return False
    return True


#==============================================================================
# Simulated Annealing Function
#==============================================================================
def calcOneStructure(loopInfo):
    InitialParams( rampedParams )
    InitialParams( highTempParams )
    protocol.initDynamics(dyn_HT,
                          potList=potList, # potential terms to use
                          bathTemp=init_t,
                          initVelocities=1,
                          finalTime=5,
                          numSteps=5000,   # whichever comes first
                          printInterval=100)
    dyn_HT.setETolerance( init_t/10 )  #used to det. stepsize. default:t/1000
    dyn_HT.run()

    # initialize parameters for cooling loop
    InitialParams( rampedParams )
    InitialParams( highTempParams )
    protocol.initDynamics(dyn,
                          potList=potList,
                          numSteps=5000,       #at each temp: 100 steps or
                          finalTime=0.5,
                          printInterval=5000)
    dyn.setETolerance( init_t/10)
    cool.run()


    # torsion angle minimization
    protocol.initMinimize(dyn,
                          printInterval=50)
    dyn.run()


    # all- atom minimization
    protocol.initMinimize(minc,
                          potList=potList,
                          dEPred=10)
    minc.run()

    #Secondary Structure Bundle Refinement:
    protocol.initDynamics(dyn_cool,
                          potList=potList, # potential terms to use
                          finalTime=0.5,
                          numSteps=5000,   # whichever comes first
                          printInterval=5000)
    dyn_cool.setETolerance(init_t/10)
    cool2.run()

    # final torsion angle minimization
    protocol.initMinimize(dyn_cool,
                          printInterval=50)
    dyn_cool.run()

    # final all- atom minimization
    protocol.initMinimize(minc_cool,
                          potList=potList,
                          dEPred=10)
    minc_cool.run()

    pass

#==============================================================================
# Run XPLOR structural calculation
#==============================================================================
from simulationTools import StructureLoop, FinalParams
StructureLoop(numStructures=numberOfStructures,
              doWriteStructures=True,
              calcMissingStructs=True,
              pdbTemplate=outFilename,
              structLoopAction=calcOneStructure,
              genViolationStats=True,
              averagePotList=potList,
              averageSortPots=[potList['BOND'],potList['ANGL'],potList['IMPR'],
                               PCS_Co_Pot],
              averageCrossTerms=refRMSD,
              averageTopFraction=0.5, #report only on best 50% of structs
              averageAccept=accept,   #only use structures which pass accept()
              averageContext=FinalParams(rampedParams),
              averageFilename="SCRIPT_ave.pdb",    #generate regularized ave structure
              averageFitSel="name CA",
              averageCompSel="not resname ANI and not name H*"     ).run()

#==============================================================================
# Write Outputs:
#==============================================================================
with open('QS_Output.txt', 'w') as out:
    out.write('\n'.join(['\t'.join([str(j) for j in i]) for i in Q_Sig_Array]))

# Write Probabilities:
for M in Metals:
    with open('Probabilities_{}.txt'.format(M), 'w') as out:
        for Res in Running_Probs[M]:
            out.write('Residue {}\n'.format(Res))
            out.write('\n'.join(['\t'.join([str(j) for j in i]) for
                                 i in Running_Probs[M][Res]]))
            out.write('\n')

for M in Metals: 
    with open('Calculate_PCS_{}.txt'.format(M), 'w') as out:
        out.write('\t'.join(list(Running_PCS[M].keys())) + '\n')
        PCS_Array = np.array([Running_PCS[M][Res] for Res in Running_PCS[M]]).T
        out.write('\n'.join(['\t'.join([str(j) for j in i]) for i in PCS_Array]))

for M in Metals: 
    with open('PCS_Probs_Final_{}.txt'.format(M), 'w') as out:
        out.write('\n'.join(['\t'.join([str(i) for i in [Res, Running_PCS[M][Res][-2],
                  Running_Probs[M][Res][-2][0], Running_Probs[M][Res][-2][1],
                  Running_Probs[M][Res][-2][2], Running_Probs[M][Res][-2][3]]])
                  for Res in sorted(Running_Probs[M])]))
