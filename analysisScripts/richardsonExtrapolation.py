#! /usr/bin/python
# Python script to perform Richardson Extrapolation on two mc-mini output files

import h5py as H5
import xml.etree.ElementTree as ET
from subprocess import Popen
from sys import argv

if (len(argv) == 1):
  print("Usage: " + argv[0] + " <input directory> <filename pattern> <output directory> <lower> <upper>")
  quit(-1)

in_directory = argv[1]
file_pattern = argv[2]
out_directory = argv[3]
lower = int(argv[4])
upper = int(argv[5])

Popen("rm " + out_directory + "/*", shell=True)

for index in range(lower, upper):
  file_1_name = in_directory + "/" + file_pattern + str(index) + ".h5"
  file_2_name = in_directory + "/" + file_pattern + str(index + 1) + ".h5"
  file_out_name = out_directory + "/" + file_pattern + "RE" + str(index) + ".h5"
  file_xdmf_name = in_directory + "/" + file_pattern + str(index) + ".xdmf"
  file_xdmf_out_name = out_directory + "/" + file_pattern + "RE" + str(index) + ".xdmf"

  file_1 = H5.File(file_1_name)
  file_2 = H5.File(file_2_name)
  file_out = H5.File(file_out_name)

  xdmf_data = ET.parse(file_xdmf_name)
  xdmf_out_file = open(file_xdmf_out_name, 'w')

  xdmf_out_root = xdmf_data.getroot().copy()
  xdmf_out_data = xdmf_out_root.find('./Domain/Grid')
  xdmf_out_data.remove(xdmf_out_data.find('./Attribute[@Name="Viscosity"]'))

  datasets = file_1.items()

  for (dataset_name, dataset_attrs) in datasets:
    if (dataset_name == "Viscosity"):
      break
    (M, N) = dataset_attrs.shape

    data_1 = file_1.get(dataset_name).value
    data_2 = file_2.get(dataset_name).value

    data_out = []
    for i in range(M):
      data_out.insert(i, [])
      for j in range(N):
        average = (data_2[2 * i    ][2 * j] + data_2[2 * i    ][2 * j + 1] + \
                   data_2[2 * i + 1][2 * j] + data_2[2 * i + 1][2 * j + 1]) / 4
        value = data_1[i][j] - average
        data_out[i].insert(j, value)

    xdmf_out_data.find('./Attribute[@Name="' + dataset_name + '"]/DataItem').text = file_pattern + "RE" + str(index) + ".h5:/" + dataset_name + "\n"

    file_out.create_dataset(dataset_name, data=data_out)
  xdmf_out_file.write('<?xml version="1.0"?>')
  xdmf_out_file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
  xdmf_out_file.write(ET.tostring(xdmf_out_root))
  xdmf_out_file.write('\n')
  xdmf_out_file.close()
