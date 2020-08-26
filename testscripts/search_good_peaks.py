from iotbx.phil import parse
from dials.util import show_mail_on_error
from dials.util.options import OptionParser
from cctbx import uctbx, miller, crystal
import sys
import GSASIIindex as gi
from libtbx import easy_mp
from cell_manager import Candidate_cell

n_triclinic = 24

help_message = "whatever"

phil_scope = parse(
    """
  search {
    lattices = cF cI cP hR hP tI tP oF oI oC oP mC mP aP
      .type = choice(multi=True)
      .help = "Bravais lattices to search"
    n_peaks = 20
      .type = int
      .help = "Number of d-spacings for unit cell search"
    wavl = 1.03
      .type = float
      .help = "GSASII wants 2th values in addition to d-spacings, so we need"
              "a wavelength. It doesn't seem to be used for anything."
    timeout = None
      .type = int
      .help = "Timeout the GSASII lattice search calls after this many seconds"
  }

  multiprocessing {
    nproc = 1
      .type = int
  }

  validate {
    method = *gsasii powder
      .type = choice
    d_min = 2
      .type = float
  }

  unit_cell = None
    .type = unit_cell
  space_group = None
    .type = space_group
  input {
    peak_list = None
      .type = str
      .help = "a list of d-spacings, 1 per line"
    powder_pattern = None
      .type = str
      .help = "A powder pattern in .xy format for validation of candidate cells"
  }
    """
)



  
def gpeak_from_d_spacing(d, wavl):
  from cctbx.uctbx import d_as_d_star_sq, d_star_sq_as_two_theta
  twoth = d_star_sq_as_two_theta(d_as_d_star_sq(d), wavl, deg=True)
  return [twoth, 1000, True, False, 0, 0, 0, d, 0]

def prepare_gpeaks(d_spacings, wavl, n_peaks=None):
  if n_peaks is not None:
    assert len(d_spacings) >= n_peaks
    trial_set = sorted(random.sample(d_spacings, n_peaks), reverse=True)
  else: 
    trial_set = d_spacings
  trial_gpeaks = [gpeak_from_d_spacing(d, wavl) for d in trial_set]
  return trial_gpeaks

def call_gsas(args):
  '''
  args is a tuple (d_spacings, bravais, powder_pattern, d_min):
  d_spacings: list of floats, the peaks for the cell search
  bravais: string, a lattice symbol like mP
  powder_pattern: a list of (x,y) tuples, the powder pattern for scoring
      candidates (x-axis must be d-spacing)
  d_min: float, d_min for scoring candidates against powder pattern
  '''

  symmorphic_sgs = ['F23', 'I23', 'P23', 'R3', 'P3', 'I4', 'P4', 'F222', 'I222',
      'A222', 'B222', 'C222', 'P222', 'I2', 'C2', 'P2', 'P1']

  d_spacings = args[0]
  bravais = args[1]
  powder_pattern = args[2]
  d_min = args[3]
  wavl = args[4]
  timeout = args[5]

  lattices = ['cF', 'cI', 'cP', 'hR', 'hP', 'tI', 'tP', 'oF', 'oI', 'oA', 'oB',
      'oC', 'oP', 'mI', 'mC', 'mP', 'aP']
  i_bravais = lattices.index(bravais)
  bravais_list = [i==i_bravais for i in range(17)]

  #TODO: adaptively set starting volume controls[3] based on number of peaks
  controls = [0, 0.0, 4, 200, 0, 'P1', 1.0, 1.0, 1.0, 90.0, 90.0, 90.0, 1.0,
      'P 1', []]

  trial_gpeaks = prepare_gpeaks(d_spacings, wavl)

  try:
    success, dmin, gcells = gi.DoIndexPeaks(trial_gpeaks, controls, bravais_list, None, timeout=timeout)
  except FloatingPointError: #this raises "invalid value encountered in double_scalars" sometimes
    print("############################################################\n"*10,
        "crash in search for {}".format(bravais))
    return []

  candidates = []
  for gcell in gcells:
    m20 = gcell[0]
    ibrav = gcell[2]
    npeaks = gcell[12]
    sg = symmorphic_sgs[ibrav]
    uc = gcell[3:9]
    cs = crystal.symmetry(unit_cell=uc, space_group_symbol=sg)
    candidate = Candidate_cell(cs, npeaks, m20)
    candidate.save_powder_score(powder_pattern, d_min)
    candidates.append(candidate)
  return candidates

def i_first_matching(cand1, cand_list):
  for i_cand, cand2 in enumerate(cand_list):
    if cand1.matches_cell(cand2):
      return i_cand
  raise RuntimeError





  


class Script(object):
  def __init__(self):
    usage = "whatever"
    self.parser = OptionParser(
         usage=usage,
         phil=phil_scope,
         epilog=help_message,
         check_format=False,
         read_reflections=True,
         read_experiments=True,
         )





        
  def run(self):
    params, options = self.parser.parse_args()


    # Load d-spacings and powder pattern from files
    with open(params.input.peak_list) as f:
      d_spacings = [float(l.strip()) for l in f.readlines()]
    with open(params.input.powder_pattern) as f:
      powder_pattern = []
      for l in f.readlines():
        x, y = l.split()
        powder_pattern.append((float(x), float(y)))
    d_min = params.validate.d_min
    wavl = params.search.wavl
    timeout = params.search.timeout

    
    lattices_todo = (
        ['cF'] * 2 +
        ['cI'] * 2 +
        ['cP'] * 2 +
        ['hR'] * 4 +
        ['hP'] * 4 +
        ['tI'] * 4 +
        ['tP'] * 4 +
        ['oF'] * 6 +
        ['oI'] * 6 +
        ['oC'] * 6 +
        ['oP'] * 6 +
        ['mC'] * 16 +
        ['mP'] * 16 +
        ['aP'] * n_triclinic
        )
    #lattices_todo = ['mC', 'mC', 'aP', 'aP']
    lattices_todo.reverse() # we want to start the longer jobs right away

    candidates = easy_mp.parallel_map(
        call_gsas,
        [(d_spacings, bravais, powder_pattern, d_min, wavl, timeout) 
            for bravais in lattices_todo],
        processes=params.multiprocessing.nproc)

    candidates_triclinic_flat, candidates_other_flat = [], []
    for c in candidates[:n_triclinic]:
      candidates_triclinic_flat.extend(c)
    for c in candidates[n_triclinic:]:
      candidates_other_flat.extend(c)

    from functools import partial

    print("Monoclinic and higher results:")
    i_in_other_function = partial(
        i_first_matching, cand_list=candidates_other_flat)
    i_in_other_list = easy_mp.parallel_map(
        i_in_other_function,
        candidates_other_flat,
        processes=params.multiprocessing.nproc)
    unique_cells_other = set(i_in_other_list)
    results_other = []
    for i_unique in unique_cells_other:
      matches = [
          (cand.powder_score/cand.m20, cand)
          for i_cand, cand in enumerate(candidates_other_flat)
          if i_in_other_list[i_cand] == i_unique
          ]
      best = min(matches, key=lambda m:m[0])
      results_other.append(best)
    results_other.sort(key=lambda r:r[0]) #, reverse=True)
    for c in results_other[:10]:
      print("{:.2f}\t{}".format(c[0], c[1]))

    print("Triclinic results:")
    i_in_triclinic_function = partial(
        i_first_matching, cand_list=candidates_triclinic_flat)
    i_in_triclinic_list = easy_mp.parallel_map(
        i_in_triclinic_function,
        candidates_triclinic_flat,
        processes=params.multiprocessing.nproc)
    unique_cells_triclinic = set(i_in_triclinic_list)
    results_triclinic = []
    for i_unique in unique_cells_triclinic:
      matches = [
          (cand.powder_score/cand.m20, cand)
          for i_cand, cand in enumerate(candidates_triclinic_flat)
          if i_in_triclinic_list[i_cand] == i_unique
          ]
      best = min(matches, key=lambda m:m[0])
      results_triclinic.append(best)
    results_triclinic.sort(key=lambda r:r[0])
    for c in results_triclinic[:10]:
      print("{:.2f}\t{}".format(c[0], c[1]))


    import IPython; IPython.embed()


    

if __name__=="__main__":
  script = Script()
  script.run()
