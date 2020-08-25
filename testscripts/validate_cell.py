from iotbx.phil import parse
from dials.util import show_mail_on_error
from dials.util.options import OptionParser
from cctbx import uctbx, miller
import sys

help_message = "haha good luck"
phil_scope = parse(
    """
  unit_cell = None
    .type = unit_cell
  space_group = None
    .type = space_group
  d_min = 2
    .type = float
  input {
    powder_pattern = None
      .type = str
  }
    """
)


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
    params, options = self.parser.parse_args(show_diff_phil=False)

    with open(params.input.powder_pattern) as f:
      data = []
      for line in f.readlines():
        x, y = [float(n) for n in line.split()]
        data.append((x, y))

    uc = params.unit_cell
    sg = params.space_group
    d_spacings = []
    mig = miller.index_generator(uc, sg.type(), 0, 0.8*params.d_min)
    for h in mig:
      d_spacings.append(uc.d(h))

    error_cumul = 0

    for x, y in data:
      best_match = min(d_spacings, key=lambda d: abs(x-d))
      error = abs(x-best_match)
      error_cumul += error*y

    print("Cumul. error: ", error_cumul)

if __name__ == "__main__":
  with show_mail_on_error():
    script = Script()
    script.run()

