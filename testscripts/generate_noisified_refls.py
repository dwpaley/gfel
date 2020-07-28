import sys
from cctbx import crystal, miller

from random import gauss
import pickle

nref = int(sys.argv[1])
noise = float(sys.argv[2])
out_name = sys.argv[3]
OMIT_I_LIST = [0]
NO_NOISE_I_LIST = [1, 2, 5, 10, 15, 20]


  
uc = (5.88, 7.3, 29.1, 90, 95.8, 90)
sg = 'C2/c'
cs = crystal.symmetry(unit_cell=uc, space_group=sg)
ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=False, d_min=1.4)
d_spacings = ms.sort().d_spacings().data()
n_per_d = int(nref / d_spacings.size())

result = []
for i in range(d_spacings.size()):
  for j in range(n_per_d):
    d = d_spacings[i]
    if i not in NO_NOISE_I_LIST: d += gauss(0, noise * d_spacings[i])
    if i not in OMIT_I_LIST: result.append(d)

with open(out_name, 'wb') as f: pickle.dump(result, f)
print("{} reflections generated from {} unique. Result saved to {}".format(
    len(result), d_spacings.size(), out_name))

