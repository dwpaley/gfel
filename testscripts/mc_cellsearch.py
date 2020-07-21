import pickle, random
import GSASIIindex as gi


REFLS = '20k_0p1_2.pkl'
N_UNIQ = 300
N_TOTAL = 10000
OUT_FILE = 'cells.pkl'
OVERLAP_TOL_FRAC = .01
N_SEARCH_PEAKS = 30
WAVL = 1.02

class Candidate_cell(object):
  def __init__(gcell_list):
    from cctbx import uctbx
    self.uc = uctbx.unit_cell(gcell_list[3:9])
    self.hits = [gcell_list]

  def matches_cell(uc2):
    return True if self.uc.similarity_transformations(uc2).size() > 0 else False

  def store_hit(gcell_list):
    self.hits.append(gcell_list)


def gpeak_from_miller_set(ms, i_peak, wavl):
  ''' Return the i-th peak in a miller set in gsas-ii format'''
  tta = ms.two_theta(wavl, deg=True).sort()
  twotheta = tta.data()[i_peak]
  d_spacing = tta.d_spacings().data()[i_peak]
  #peaks are [2th, I, use_flag, indexed_flag, h, k, l, d(obs), d(calc)]
  return [twotheta, 1000, True, False, 0, 0, 0, d_spacing, 0]

def gpeak_from_d_spacing(d, wavl):
  from cctbx.uctbx import d_as_d_star_sq, d_star_sq_as_two_theta
  twoth = d_star_sq_as_two_theta(d_as_d_star_sq(d), wavl, deg=True)
  return [twoth, 1000, True, False, 0, 0, 0, d, 0]


with open(REFLS, 'rb') as f: refls = pickle.load(f)


while True:
  trial_set = random.choices(refls, k=N_UNIQ)

  # Filter out overlapping reflections. Randomly filter one of the pair
  # to avoid bias
  to_skip = []
  for i in range(1, len(trial_set)):
    if (abs(trial_set[i]-trial_set[i-1]) / trial_set[i]) < OVERLAP_TOL_FRAC:
      to_skip.append(i-1 if random.random()>.5 else i)
  trial_set = [trial_set[i] for i in range(len(trial_set)) if i not in to_skip]
  trial_set = trial_set[:N_SEARCH_PEAKS]
  trial_peaks = [gpeak_from_d_spacing(d, WAVL) for d in trial_set]

  # lattice codes at the bottom of the file
  bravais = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,True,0,0]
  # The controls are magic
  controls = [0, 0.0, 4, 200, 0, 'P1', 1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 1.0, 'P 1', []]

  success, dmin, cells = gi.DoIndexPeaks(trial_peaks, controls, bravais, None)

  from cctbx import uctbx
  with open(OUT_FILE, 'rb+') as f:
    candidates = pickle.load(f)
    for cell in cells:
      done = False
      uc = uctbx.unit_cell(cell[3:9])
      for cand in candidates:
        if cand.matches_cell(uc):
          cand.store_hit(cell)
          done = True
          break
      if not done:
        candidates.append(Candidate_cell(cell))
    pickle.dump(candidates, f)

      




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
