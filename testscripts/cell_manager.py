from cctbx import miller, crystal


class Candidate_cell(object):
  def __init__(self, cs, npeaks=None, m20=None):
    '''
    Constructed from an instance of cctbx.crystal.symmetry. Optionally, npeaks
    is the GSASIIindex quantity Nc, which is the number of peaks generated by
    this cell and lattice within a given resolution limit.
    '''
    self.cs = cs
    self.niggli_uc = cs.niggli_cell().unit_cell()
    self.npeaks = npeaks
    self.m20 = m20
#    self.hit_count = 1
#    self.hits = [gcell]
#    self.best_score = self.cumul_score = self.hit_score(gcell)

  def __str__(self):
    # cs.best_cell is really more like "better cell" so we call it a few
    # times to ensure we get the actual best cell
    uc = self.cs.best_cell().best_cell().best_cell().best_cell().unit_cell()
    return "{}\t{}".format(str(uc), str(self.sg.info()))

  def standardize(self):
    self.cs = self.cs.best_cell()

  def matches_cell(self, cell2):
    nc1 = self.niggli_uc
    nc2 = cell2.niggli_uc
    return nc1.similarity_transformations(nc2).size() > 0

  def calc_powder_score(self, powder_pattern, d_min):
    '''Take a list of (d, counts) tuples and return a figure of merit (lower-better)
    '''
    assert self.sg is not None
    mig = miller.index_generator(self.uc, self.sg.type(), 0, 0.8*d_min)
    d_spacings = []
    for h in mig: d_spacings.append(self.uc.d(h))

    error_cumul = 0
    for x, y in powder_pattern:
      if x < d_min: break
      best_match = min(d_spacings, key=lambda d: abs(x-d))
      error = abs(x-best_match)/x
      error_cumul += error*y

    return error_cumul

  def save_powder_score(self, powder_pattern, d_min):
    self.powder_score = self.calc_powder_score(powder_pattern, d_min)


  @property
  def uc(self):
    return self.cs.unit_cell()

  @property
  def sg(self):
    return self.cs.space_group()
    

  # Down here is Monte Carlo-related stuff we might not use anymore
  def store_hit(self, gcell):
    self.hits.append(gcell)
    self.hit_count += 1
    score = self.hit_score(gcell)
    self.cumul_score += score
    if score > self.best_score: self.best_score = score

  def average_score(self):
    return self.cumul_score / self.hit_count

  @staticmethod
  def hit_score(gcell):
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
