#! /usr/bin/python
# Python script to perform Richardson Extrapolation on two mc-mini output files

import h5py as h5
import sys
from subprocess import call
sys.path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
import visit as v

if (len(sys.argv) == 1):
  print("Usage: " + sys.argv[0] + " <filename1> <filename2> <outfile>")
  quit(-1)
else:
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]
  outfilename = sys.argv[3]

file1 = h5.File(filename1)
file2 = h5.File(filename2)
outfile = h5.File(outfilename)

datasets = file1.items()

for (dataset_name, dataset_attrs) in datasets:
  if (dataset_name == "Viscosity"):
    break
  (M, N) = dataset_attrs.shape

  data1 = file1.get(dataset_name).value
  data2 = file2.get(dataset_name).value

  dataout = []
  for i in range(M):
    dataout.insert(i, [])
    for j in range(N):
      average = (data2[2 * i    ][2 * j] + data2[2 * i    ][2 * j + 1] + \
                 data2[2 * i + 1][2 * j] + data2[2 * i + 1][2 * j + 1]) / 4
      value = data1[i][j] - average
      dataout[i].insert(j, value)

  outfile.create_dataset(dataset_name, data=dataout)

