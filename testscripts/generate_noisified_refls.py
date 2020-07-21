from cctbx import crystal, miller

from random import gauss
import pickle

nref = int(sys.argv[1])
noise = float(sys.argv[2])
out_name = sys.argv[3]


  
uc = (5.88, 7.3, 29.1, 90, 95.8, 90)
sg = 'C2/c'
cs = crystal.symmetry(unit_cell=uc, space_group=sg)
ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=False, d_min=1.4)
d_spacings = ms.sort().d_spacings().data()
n_per_d = int(nref / d_spacings.size())

result = []
for d_spacing in d_spacings:
  for i in range(n_per_d):
    d = d_spacing + gauss(0, noise * d_spacing)
    result.append(d)

with open(out_name, 'wb') as f: pickle.dump(result, f)
print("{} reflections generated from {} unique. Result saved to {}".format(
    len(result), d_spacings.size(), out_name))

