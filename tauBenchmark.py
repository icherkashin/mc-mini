#! /usr/bin/python
import sys
from subprocess import call
from math import sqrt, log
sys.path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
from visit import *

if (len(sys.argv) > 1):
  lower = int(sys.argv[1])
  upper = int(sys.argv[2])
else:
  lower = 1
  upper = 7

call(["rm", "output/tauBenchmark/tauBenchmark*"])
for refinement in range(lower, upper + 1):
  call(["./build/mc-mini", "paramFiles/tauBenchmark/tauBenchmark" + str(refinement)])

AddArgument("-nowin")

Launch()

OpenDatabase("output/tauBenchmark/tauBenchmark*.xdmf database")

DefinePythonExpression("TauBenchmarkUVelocity", \
                       ["UVelocity"], \
                       file="output/tauBenchmark/analysisScripts/tauBenchmarkUVelocity.py")
DefinePythonExpression("TauBenchmarkVVelocity", \
                       ["VVelocity"], \
                       file="output/tauBenchmark/analysisScripts/tauBenchmarkVVelocity.py")
DefineScalarExpression("UVelocityError", \
                       "abs(UVelocity - TauBenchmarkUVelocity)")
DefineScalarExpression("VVelocityError", \
                       "abs(VVelocity - TauBenchmarkVVelocity)")

DefineVectorExpression("VelocityError", \
                       "{UVelocityError, VVelocityError}")
DefineScalarExpression("ErrorMagnitude", \
                       "magnitude(VelocityError)")

DefineScalarExpression("ErrorSquare",
                       "ErrorMagnitude * ErrorMagnitude")

AddPlot("Pseudocolor", "ErrorMagnitude")

file = open('tauBenchmarkConvergence', 'w')
file.write("refinement\tMax Norm\t\tOne Norm\t\tTwo Norm\n")
DrawPlots()

oldMax = 0
oldOne = 0
oldTwo = 0

for state in range(TimeSliderGetNStates()):
  Query("NumZones")
  nZones = GetQueryOutputValue()

  refinement = lower + state

  Query("Max")
  maxNorm = GetQueryOutputValue()
  if (oldMax != 0):
    maxRate = log(oldMax / maxNorm) / log(2)
  else:
    maxRate = 0
  oldMax = maxNorm

  Query("Variable Sum")
  oneNorm = GetQueryOutputValue() / nZones
  if (oldOne != 0):
    oneRate = log(oldOne / oneNorm) / log(2)
  else:
    oneRate = 0
  oldOne = oneNorm

  AddPlot("Pseudocolor", "ErrorSquare")
  DrawPlots()
  Query("Variable Sum")
  twoNorm = sqrt (GetQueryOutputValue() / (nZones))
  if (oldTwo != 0):
    twoRate = log(oldTwo / twoNorm) / log(2)
  else:
    twoRate = 0
  oldTwo = twoNorm

  DeleteActivePlots()
  DrawPlots()
  file.write('{}\t\t{:.8f}\t{:.2f}\t{:.8f}\t{:.2f}\t{:.8f}\t{:.2f}\n'.format(refinement, maxNorm, maxRate, oneNorm, oneRate, twoNorm, twoRate))
  TimeSliderNextState()

file.close()

call(["cat", "tauBenchmarkConvergence"])
