#!/bin/python3
################################################################################

#
import getopt
import math
import numpy
import copy
import time
#
import os
import sys
#include = os.environ['PYLIB']
#sys.path.append(include)

from bpstandard import oFileData as oFileData
from bpstandard import oStrings as oStrings
from bpstandard import general as general
from structures import structures as structures
from pwIn import pwIn as pwIn
from prepPwIn import prepPwIn as prepPwIn
from pwOut import pwOut as pwOut
from runPw import runPw as runPw
from randnum import RandomLCG as RandomLCG
from randnum import RandDist as RandDist
from gnuplot import gnuplot as gnuplot



###########################
# pwHeat
###########################
class pwHeat:
  iCount = 0  # Instance counter

  def __init__(self):
    self.cmdLog = []
    self.resultsArr = []
    self.startTime = time.time()


  def controlFile(self):
    print ("Reading control file")
    print ("==================================================")

    self.args = str(sys.argv)   # Store input argument as the templateFile
    self.cFileName = sys.argv[1]
    self.cFile = oFileData()
    self.cFile.loadFile(self.cFileName)
    self.tmpDir = "tmp"
    self.runType = "full"
    self.copies = 2
    self.bmstrain = 0.05
    self.orstrain = 0.02
    self.testrain = 0.02
    self.bmsteps = 10
    self.orsteps = 10
    self.testeps = 10
    self.outdir = ""
    self.ppdir = ""
    self.pwtemplate = None
    self.degauss = 0.1
    self.cell_xx = 1.0
    self.cell_xy = 0.0
    self.cell_xz = 0.0
    self.cell_yx = 0.0
    self.cell_yy = 1.0
    self.cell_yz = 0.0
    self.cell_zx = 0.0
    self.cell_zy = 0.0
    self.cell_zz = 1.0

    self.calculations = 10
    self.maxVariation = 0.0

    self.atomList = []
    self.massList = []
    self.ppList = []

    for i in range(0,self.cFile.lineCount):
      # Directories
      #######################
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$TMPDIR", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.tmpDir = fileRowArr[1]
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$CALCULATIONS", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.calculations = int(fileRowArr[1])
      keywordResult = self.checkKeyword(self.cFile.fileData[i], "$MAXVARIATION", True)
      if(keywordResult is not None):
        fileRowArr = keywordResult.split(" ")
        self.maxVariation = float(fileRowArr[1])

    print ("==================================================")
    print()
    print(self.atomList)

    ## Mk dir
    general.mkDir(self.tmpDir)
    ## Type of run
    self.runCode = 0 # Default
    if(self.runType.upper()=="FULL" or self.runType[0:1]=="1"):       # Full run from fresh
      self.runCode = 1
      self.runType = "FULL"
    if(self.runType.upper()=="CONTINUE" or self.runType[0:1]=="2"):   # Completes where there's no output file
      self.runCode = 2
      self.runType = "CONTINUE"
    if(self.runType.upper()=="NONE" or self.runType[0:1]=="3"):       # Don't run DFT
      self.runCode = 0


  @staticmethod
  def checkKeyword(lineIn, keyword, verbose=False):
    result = None
    if(lineIn != ""):
      lineInArr = lineIn.split("#")
      lineIn = lineInArr[0]
      if(lineIn != ""):
        lineInUC = lineIn.upper()
        keywordUC = keyword.upper()
        keywordLen = len(keywordUC)
        if(lineInUC[0:keywordLen]==keywordUC):
          result = oStrings.removeDouble(lineIn," ")
          if(verbose):
            print(result)
    return result


  def run(self):
    print ("Run")
    self.makeStructure()
    self.runCalcs()

  def makeStructure(self):       # Use vc-relax to fine optimum settings
    print ("1. Make atom structure")

    self.prepPwIn = prepPwIn(self.cFileName)
    self.prepPwIn.make()
    self.pwIn = self.prepPwIn.getPW()
    self.pwIn.outputFile(self.tmpDir+"/template.in")
    self.pwIn.extractData()
    self.nat = self.pwIn.getNAT()


  def runCalcs(self):       # Use vc-relax to fine optimum settings
    print ("2. Run")

    self.maxVariationList = []
    self.energyList = []
    self.forceList = []

    rdA = RandDist()
    rdA.randomSeed()
    rdA.flat()

    print ("")
    print ("")
    print ("Maximum Variation")
    variationFile = ""

    for i in range(0,self.calculations):
      j = i + 1
      dist = rdA.rng() * self.maxVariation
      self.maxVariationList.append(dist)
      print (j, dist)
      variationFile = variationFile + str(j) + "  " + str(dist) + "\n"

      pwHeatIn = pwIn()
      pwHeatIn.loadFile(self.tmpDir+"/template.in")
      pwHeatIn.addComment("Max variation: "+str(self.maxVariation)+"  "+str(dist))
      pwHeatIn.heat(dist)
      if(i<10):
        fileName = "heat_00"+str(j)
      elif(i<100):
        fileName = "heat_0"+str(j)
      else:
        fileName = "heat_"+str(j)
      # remove any old files
      try:
        os.remove(self.tmpDir+"/"+fileName+".in")
      except OSError:
        pass
      try:
        os.remove(self.tmpDir+"/"+fileName+".out")
      except OSError:
        pass
      # Make file
      pwHeatIn.outputFile(self.tmpDir+"/"+fileName+".in")
      # run
      runPw.run(fileName,self.tmpDir,1)
      # read
      pwHeatOut = pwOut(self.tmpDir+"/"+fileName+".out")
      self.energyList.append(pwHeatOut.getEnergyPerAtom())
      self.forceList.append(pwHeatOut.getForce())

    # Write variation file
    wFile = open(self.tmpDir+"/variation.dat", 'w')
    wFile.write(variationFile)

    # Make plots
    gp = gnuplot()
    gp.reset()
    gp.setDir(self.tmpDir+"/plots")
    gp.title("Energy - Heat")
    gp.axisLabel("x1","Maximum Variation (Bohr)")
    gp.axisLabel("y1","Energy (Ry)")
    gp.outputPlot("energy-heating")
    gp.addPlot(self.maxVariationList,self.energyList,"x1y1","EoS",1,1)
    gp.setPlotType(gp.getLastPlotID(),"points")
    gp.makePlot()

    # Force
    gp.reset()
    gp.setDir(self.tmpDir+"/plots")
    gp.title("Force - Heat")
    gp.axisLabel("x1","Maximum Variation (Bohr)")
    gp.axisLabel("y1","Force (Ry/Bohr)")
    gp.outputPlot("force-heating")
    gp.addPlot(self.maxVariationList,self.forceList,"x1y1","EoS",1,1)
    gp.setPlotType(gp.getLastPlotID(),"points")
    gp.makePlot()

########################################################################
## Make BP calc and run it

newHeat = pwHeat()
newHeat.controlFile()
newHeat.run()









################################################################################
