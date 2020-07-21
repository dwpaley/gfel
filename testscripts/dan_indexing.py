
import GSASIIindex as gi
from cctbx import crystal, miller

from random import gauss


def gpeak_from_miller_set(ms, i_peak, wavl):
  ''' Return the i-th peak in a miller set in gsas-ii format'''
  print("gpeak")
  tta = ms.two_theta(wavl, deg=True).sort()
  twotheta = tta.data()[i_peak]
  #twotheta += gauss(0, twotheta*.005)
  d_spacing = tta.d_spacings().data()[i_peak]
  if i_peak not in [0,1]: d_spacing += gauss(0, d_spacing*.005)
  #peaks are [2th, I, use_flag, indexed_flag, h, k, l, d(obs), d(calc)]
  return [twotheta, 1000, True, False, 0, 0, 0, d_spacing, 0]

  
uc = (11,9,8,90,102,90)
sg = 'P21'
uc = (5.88, 7.3, 29.1, 90, 95.8, 90)
uc = (5.88, 7.3, 27.1, 90, 95.8, 90)
sg = 'C2/c'
# lattice codes at the bottom of the file
bravais = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,True,0,0]
wavl = 1.02
cs = crystal.symmetry(unit_cell=uc, space_group=sg)
ms = miller.build_set(crystal_symmetry=cs, anomalous_flag=False, d_min=.8)

peaks = [gpeak_from_miller_set(ms, i, wavl) for i in range(50)]

#Delete peaks that are too close together
tt_overlap_tol = .1
i_to_delete = []
for i in range(1, len(peaks)):
  if abs(peaks[i-1][0] - peaks[i][0]) < tt_overlap_tol:
    i_to_delete.append(i)
peaks = [peaks[i] for i in range(len(peaks)) if i not in i_to_delete]
peaks = peaks[0:30]
for p in peaks: print(p)


# The controls are magic
controls = [0, 0.0, 4, 200, 0, 'P1', 1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 1.0, 'P 1', []]


success, dmin, cells = gi.DoIndexPeaks(peaks, controls, bravais, None)
cells.sort(key=lambda x: x[2]/x[0])
for c in cells[0:10]: print(c)

import ipdb;ipdb.set_trace()




lattices = """
            * 0 F cubic
            * 1 I cubic
            * 2 P cubic
            * 3 R hexagonal (trigonal not rhombohedral)
            * 4 P hexagonal
            * 5 I tetragonal
            * 6 P tetragonal
            * 7 F orthorhombic
            * 8 I orthorhombic
            * 9 A orthorhombic
            * 10 B orthorhombic
            * 11 C orthorhombic
            * 12 P orthorhombic
            * 13 I monoclinic
            * 14 C monoclinic
            * 15 P monoclinic
            * 16 P triclinic
            """
