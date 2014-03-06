#! /usr/bin/python
import sys
from subprocess import call
sys.path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
from visit import *

if (len(sys.argv) > 1):
  lower = sys.argv[1]
  upper = sys.argv[2]
else:
  lower = 2
  upper = 6

for refinement in range(int(lower), int(upper) + 1):
  call(["./build/mc-mini", "paramFiles/solCXBenchmark/solCXBenchmark" + str(refinement)])

Launch()

dbFile = "output/solCXBenchmark/solCXBenchmark*.xdmf database"

OpenDatabase(dbFile)
md = GetMetaData(dbFile)

DefineVectorExpression("Velocity", "{UVelocity, VVelocity}")
DefineVectorExpression("Interpolation", "Velocity - pos_cmfe(<[1]id:Velocity>, mesh, 0)")
AddPlot("Vector", "Interpolation")

DrawPlots()

nStates = md.numStates
for state in range(nStates):
  print state
  TimeSliderNextState()
  Query("Max")

