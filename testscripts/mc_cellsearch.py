import pickle, random


REFLS = 'test1.pkl'
N_UNIQ = 300
N_TOTAL = 10000
OUT_FILE = 'cells.pkl'
OVERLAP_TOL_FRAC = .01
N_SEARCH_PEAKS

def gpeak_from_miller_set(ms, i_peak, wavl):
  ''' Return the i-th peak in a miller set in gsas-ii format'''
  tta = ms.two_theta(wavl, deg=True).sort()
  twotheta = tta.data()[i_peak]
  #twotheta += gauss(0, twotheta*.005)
  d_spacing = tta.d_spacings().data()[i_peak]
  d_spacing += gauss(0, d_spacing*.002)
  #peaks are [2th, I, use_flag, indexed_flag, h, k, l, d(obs), d(calc)]
  return [twotheta, 1000, True, False, 0, 0, 0, d_spacing, 0]

with open(REFLS, 'rb') as f: refls = pickle.load(f)


while True:
  trial_set = random.choices(refls, N_UNIQ)

  to_skip = []
  for i in range(1, len(trial_set)):
    if (abs(trial_set[i]-trial_set[i-1]) / trial_set[i]) < OVERLAP_TOL_FRAC:
      to_skip.append(i-1 if random.random()>.5 else i)
  trial_set = [trial_set[i] for i in range(len(trial_set)) if i not in to_skip]
  trial_set = trial_set[:N_SEARCH_PEAKS]




