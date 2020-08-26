import pickle, random
import sys
import GSASIIindex as gi
from libtbx import easy_mp
from cctbx import uctbx

import warnings
warnings.simplefilter("error")
warnings.warn("mc_cellsearch is not in active use. Currently search_good_peaks" +
    " is the only option. The Candidate_cell and ..._manager classes are in" +
    " cell_manager.py.", DeprecationWarning)


REFLS = '../d_table.txt'
KNOWN_GOOD = '../known_good.txt'
N_UNIQ = 300
N_TOTAL = 10000
CELL_FILE = 'cells.pkl'
OVERLAP_TOL_FRAC = .01
N_SEARCH_PEAKS = 20
WAVL = 1.03232
NPROC = 60


class Candidate_cell(object):
  def __init__(self, gcell, latt_symbol=None):
    '''
    latt_symbol is a 2-letter string, e.g. "mP"
    '''
    from cctbx import uctbx, sgtbx
    self.uc = uctbx.unit_cell(gcell[3:9])
    if latt_symbol: 
      centering = latt_symbol[1]
      sg_string = centering + '1'
      self.sg = sgtbx.space_group(sg_string)
    else:
      self.sg = None
    self.hit_count = 1
    self.hits = [gcell]
    self.best_score = self.cumul_score = self.hit_score(gcell)

  def matches_cell(self, uc2):
    return (self.uc.similarity_transformations(uc2).size() > 0 and self.sg==uc2.sg)

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

  def powder_score(self, powder_pattern, d_min):
    '''Take a list of (d, counts) tuples and return a figure of merit (lower-better)
    '''
    assert self.sg is not None
    from cctbx import miller
    mig = miller.index_generator(self.uc, self.sg.type(), 0, 0.8*d_min)
    d_spacings = []
    for h in mig: d_spacings.append(self.uc.d(h))

    error_cumul = 0
    for x, y in data:
      best_match = min(d_spacings, key=lambda d: abs(x-d))
      error = abs(x-best_match)
      error_cumul += error*y

    return error_cumul
    

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

def prepare_trial_peaks(refls, known_good_all):
  # Throw out 1/3 of known_good
  known_good = random.sample(known_good_all, len(known_good_all)*2//3)

  trial_set = random.choices(refls, k=N_UNIQ)
  trial_set.sort(reverse=True)

  # Find groups of overlapping d-spacings. Keep only 1 of each group.
  i = 0
  i_to_keep = []
  while i<len(trial_set):
    i_matches = []
    datum = trial_set[i]
    while i<len(trial_set) and ((datum - trial_set[i]) / datum) < OVERLAP_TOL_FRAC:
      i_matches.append(i)
      i += 1
    i_to_keep.append(random.choice(i_matches))
  trial_set = [trial_set[i] for i in range(len(trial_set)) if i in i_to_keep]
      
  # Now filter out overlaps between random and known good.
  i_to_skip = []
  for i in range(len(trial_set)):
    for d in known_good:
      if (abs(trial_set[i]-d) / trial_set[i]) < OVERLAP_TOL_FRAC:
        i_to_skip.append(i)
        break
  trial_set = [trial_set[i] for i in range(len(trial_set)) if i not in i_to_skip]


  

#  # Filter out overlaps between random refls. Randomly filter one of the pair
#  # to avoid bias
#  to_skip = []
#  for i in range(1, len(trial_set)):
#    if (abs(trial_set[i]-trial_set[i-1]) / trial_set[i]) < OVERLAP_TOL_FRAC:
#      to_skip.append(i-1 if random.random()>.5 else i)
#  # Now filter out overlaps between random and known good.
#  for i in range(len(trial_set)):
#    for d in known_good:
#      if (abs(trial_set[i]-d) / trial_set[i]) < OVERLAP_TOL_FRAC:
#        to_skip.append(i)
#        break
#  trial_set = [trial_set[i] for i in range(len(trial_set)) if i not in to_skip]

  trial_set = random.sample(trial_set, N_SEARCH_PEAKS-len(known_good))
  trial_set.extend(known_good)
  trial_set.sort(reverse=True)
  trial_peaks = [gpeak_from_d_spacing(d, WAVL) for d in trial_set]
  return trial_peaks


def prepare_good_peaks(known_good_all):
  assert len(known_good_all) >= N_SEARCH_PEAKS
  trial_set = sorted(random.sample(known_good_all, N_SEARCH_PEAKS), reverse=True)
  trial_peaks = [gpeak_from_d_spacing(d, WAVL) for d in trial_set]
  return trial_peaks


def call_gsas(args):
  ''' 
  refls and known_good are lists of d-spacings. The idea is that refls
  will be the full list of noisy d-spacings measured from individual frames,
  while known_good will be a short list of sharp peaks accurately measured
  from the radial average. We will use the full list of known_good plus a
  random sample of refls.
  '''

  
  refls = args[0]
  weights = args[1]
  known_good_all = args[2]
  min_score = args[3]

  #trial_peaks = prepare_trial_peaks(refls, known_good_all)
  trial_peaks = prepare_good_peaks(known_good_all)


  # lattice codes at the bottom of the file
  bravais = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,True,0,0]
  # The controls are magic
  controls = [0, 0.0, 4, 800, 0, 'P1', 1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 1.0, 'P 1', []]

  success, dmin, cells = gi.DoIndexPeaks(trial_peaks, controls, bravais, None)

  if min_score:
    cells = [c for c in cells if Candidate_cell.hit_score(c) > min_score]
  return cells


#comm.gather(all_cells, dest=0)
#comm.reduce(n, dest=0)
#if rank==0:
#  pass
def run():
  with open(REFLS) as f:
    refls, weights = [], []
    for line in f.readlines():
      d, w = [float(x) for x in line.strip().split(',')]
      refls.append(d)
      weights.append(w)
  with open(KNOWN_GOOD) as f:
    known_good = [float(line.strip()) for line in f.readlines()]

  call_gsas((refls, weights, known_good, None))
  quit()
#  import profile
#  profile.runctx('call_gsas((refls, weights, known_good, None))', globals(), locals(), filename='call_gsas.prof')
#  quit()




  current_cells = easy_mp.parallel_map(
      call_gsas, 
      [(refls, weights, known_good, None) for _ in range(NPROC)], 
      processes=NPROC)

  cell_man = Candidate_cell_manager()

  current_cells_flat = []
  for l in current_cells: current_cells_flat.extend(l)
  for gcell in current_cells_flat:
    cell_man.store_cell(gcell)
  cell_man.maintain(force=True)
  min_score = cell_man.min_score


  current_cells = easy_mp.parallel_map(
      call_gsas,
      [(refls, weights, known_good, min_score) for _ in range(NPROC*3)],
      processes=NPROC)


  current_cells_flat = []
  for l in current_cells: current_cells_flat.extend(l)
  for gcell in current_cells_flat:
    cell_man.store_cell(gcell)
  cell_man.maintain(force=True)

  with open(sys.argv[1], 'wb') as f: pickle.dump(cell_man, f)


if __name__=='__main__':
  from ipdb import launch_ipdb_on_exception
  with launch_ipdb_on_exception():
    run()


def lattices():
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
