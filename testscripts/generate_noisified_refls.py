import sys
from cctbx import crystal, miller

from random import gauss
import pickle

nref = int(sys.argv[1])
noise = float(sys.argv[2])
out_root = sys.argv[3]
OMIT_I_LIST = [0]
GOOD_I_LIST = [1, 2, 5, 10, 15, 20]


all_name = out_root + "_all.pkl"
good_name = out_root + "_good.pkl"
  
uc = (5.88, 7.3, 29.1, 90, 95.8, 90)
sg = 'C2/c'
cs = crystal.symmetry(unit_cell=uc, space_group=sg)
ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=False, d_min=1.4)
d_spacings = ms.sort().d_spacings().data()
n_per_d = int(nref / d_spacings.size())



all_refls = []
for i in range(d_spacings.size()):
  for j in range(n_per_d):
    d = d_spacings[i]
    if i not in OMIT_I_LIST: all_refls.append(d)

good_refls = [d_spacings[i] for i in GOOD_I_LIST]

with open(all_name, 'wb') as f: pickle.dump(all_refls, f)
with open(good_name, 'wb') as f: pickle.dump(good_refls, f)
print("{} reflections generated from {} unique".format(
    len(all_refls), d_spacings.size()))

