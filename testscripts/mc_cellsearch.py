import pickle, random
import sys
import GSASIIindex as gi
from libtbx import easy_mp
from cctbx import uctbx


REFLS = '0p5_1.pkl'
N_UNIQ = 300
N_TOTAL = 10000
CELL_FILE = 'cells.pkl'
OVERLAP_TOL_FRAC = .01
N_SEARCH_PEAKS = 30
WAVL = 1.02

class Candidate_cell(object):
  def __init__(self, gcell):
    from cctbx import uctbx
    self.uc = uctbx.unit_cell(gcell[3:9])
    self.hit_count = 1
    self.hits = [gcell]
    self.best_score = self.cumul_score = self.hit_score(gcell)

  def matches_cell(self, uc2):
    return True if self.uc.similarity_transformations(uc2).size() > 0 else False

  def store_hit(self, gcell):
    self.hits.append(gcell)
    self.hit_count += 1
    score = self.hit_score(gcell)
    self.cumul_score += score
    if score > self.best_score: self.best_score = score

  def average_score(self):
    return self.cumul_score / self.hit_count

  @classmethod
  def hit_score(cls, gcell):
    m20, x20, nc = gcell[0:3]
    return m20/nc/(x20+1)

class Candidate_cell_manager(object):
  def __init__(self):
    self.cells = []
    self.min_score = 0

  def maintain(self, force=False):
    if len(self.cells) > 30 or force:
      self.cells.sort(key=lambda x: x.cumul_score, reverse=True)
      self.cells = self.cells[:20]
      self.min_score = min([c.average_score() for c in self.cells])

  def store_cell(self,gcell):
    self.maintain()
    uc = uctbx.unit_cell(gcell[3:9])
    score = Candidate_cell.hit_score(gcell)
    if score > self.min_score:
      found_match = False
      for cell in self.cells:
        if cell.matches_cell(uc):
          cell.store_hit(gcell)
          found_match = True
          break
      if not found_match:
        self.cells.append(Candidate_cell(gcell))

    




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




def call_gsas(min_score=None):
  trial_set = random.choices(refls, k=N_UNIQ)
  trial_set.sort(reverse=True)

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

  if min_score:
    cells = [c for c in cells if Candidate_cell.hit_score(c) > min_score]
  return cells


#comm.gather(all_cells, dest=0)
#comm.reduce(n, dest=0)
#if rank==0:
#  pass
if __name__=='__main__':
  with open(REFLS, 'rb') as f: refls = pickle.load(f)

  current_cells = easy_mp.parallel_map(
      call_gsas, 
      [None for _ in range(32)], 
      processes=32)

  cell_man = Candidate_cell_manager()

  current_cells_flat = []
  for l in current_cells: current_cells_flat.extend(l)
  for gcell in current_cells_flat:
    cell_man.store_cell(gcell)
  cell_man.maintain(force=True)
  min_score = cell_man.min_score


  current_cells = easy_mp.parallel_map(
      call_gsas,
      [min_score for _ in range(1000)],
      processes=32)


  current_cells_flat = []
  for l in current_cells: current_cells_flat.extend(l)
  for gcell in current_cells_flat:
    cell_man.store_cell(gcell)
  cell_man.maintain(force=True)

  with open(sys.argv[1], 'wb') as f: pickle.dump(cell_man, f)

      




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
