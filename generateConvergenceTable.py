#! /usr/bin/python
# Generate a convergence table for a collection of output files by using
# Richardson Extrapolation to estimate the errors.
#
# usage: ./generateConvergenceTable.py <input directory> <file pattern> <lower> <upper> <output file>

from math import sqrt, log
from subprocess import call
from sys import argv, path
path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
from visit import *

if (len(argv) == 1):
  print("Usage: " + argv[0] + " <input directory> <file pattern> <output file> <lower> <upper>")
  quit(-1)

in_directory = argv[1]
file_pattern = argv[2]
out_file = argv[3]
lower = int(argv[4])
upper = int(argv[5])

call(["./richardsonExtrapolation.py", in_directory, file_pattern, in_directory + "/richardsonExtrapolation/", str(lower), str(upper)])

Launch()

OpenDatabase(in_directory + "/richardsonExtrapolation/" + file_pattern + "RE*.xdmf database")

DefineScalarExpression("Velocity", \
                       "sqrt(UVelocity * UVelocity + VVelocity * VVelocity)")
DefineScalarExpression("VelocityMagnitude", \
                       "UVelocity * UVelocity + VVelocity * VVelocity")

DefineScalarExpression("TemperatureMagnitude", \
                       "abs(Temperature)")

DefineScalarExpression("TemperatureSquare", \
                       "Temperature * Temperature");

AddPlot("Pseudocolor", "TemperatureMagnitude")

file = open(out_file, 'w')
file.write("refinement\tMaxNorm\t\tOneNorm\t\tTwoNorm\n")
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
  if (oldMax != 0) and (maxNorm != 0):
    maxRate = log(oldMax / maxNorm) / log(2)
  else:
    maxRate = 0
  oldMax = maxNorm

  Query("Variable Sum")
  oneNorm = GetQueryOutputValue() / nZones
  if (oldOne != 0) and (oneNorm != 0):
    oneRate = log(oldOne / oneNorm) / log(2)
  else:
    oneRate = 0
  oldOne = oneNorm

  AddPlot("Pseudocolor", "TemperatureSquare")
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

call(["cat", out_file])
