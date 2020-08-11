import logging
from iotbx.phil import parse

from scitbx.array_family import flex

from dials.util import log
from dials.util import show_mail_on_error
from dials.util.options import OptionParser
from dials.util.version import dials_version

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

logger = logging.getLogger("dials.command_line.powder_from_spots")

help_message = """

Nobody can help you
"""

phil_scope = parse(
    """
  file_path = None
    .type = str
    .multiple = True
    .help = Files to read
  n_bins = 0
    .type = int
    .help = Number of bins in the radial average. Auto determined if set to 0
  d_max = 20
    .type = float
  d_min = 1.4
    .type = float
  verbose = True
    .type = bool
    .help = Extra logging information
  output_bins = True
    .type = bool
    .help = Whether to print values for each bin
  output_file = None
    .type = str
    .help = Output file for logging results
  plot_x_max = None
    .type = int
    .help = Max value for x axis
  plot_y_max = None
    .type = int
    .help = Max value for xyaxis
  low_max_two_theta_limit = None
    .type = float
    .help = Low two theta cutoff
  normalize = False
    .type = bool
    .help = Whether to normalize the Y values to 1
  show_plots = True
    .type = bool
    .help = Whether to show the radial average plot
  mask = None
    .type = str
    .help = DIALS style pixel mask. Average will skip these pixels
  median_filter_size = None
    .type = int
    .help = If not none, applies a scipy.ndimage median_filter to the average
  x_axis = *two_theta q resolution
    .type = choice
    .help = Units for x axis
  image_number = None
    .type = int
    .help = When supplying a composite file, which image to show. Otherwise \
            shows all images in the file.
  panel = None
    .type = int
    .help = Only use data from the specified panel
  max_images = None
    .type = int
    .help = When supplying a composite file, only show up to max_images images
  reference_geometry = None
    .type = path
    .help = Apply this geometry before creating average
  unit_cell = None
    .type = unit_cell
    .help = Show positions of miller indices from this unit_cell and space \
            group.
  space_group = None
    .type = space_group
    .help = Show positions of miller indices from this unit_cell and space \
            group.
output {
  log = dials.powder_from_spots.log
    .type = str
}
"""
)


class Script(object):
  def __init__(self):
    usage = "usage: make a powder pattern"
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

    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())


    assert len(params.input.reflections) == 1, "Please supply 1 reflections file"
    assert len(params.input.experiments) == 1, "Please supply 1 experiments file"

    # setup limits and bins
    assert params.n_bins, "Please supply n_bins for the pattern"
    n_bins = params.n_bins
    d_max, d_min = params.d_max, params.d_min
    d_inv_low, d_inv_high = 1/d_max, 1/d_min

    sums = flex.double(params.n_bins)

    refls = params.input.reflections[0].data
    expts = params.input.experiments[0].data

    import random
    print(len(expts))
    for i, expt in enumerate(expts):
        if random.random() < 0.01: print("experiment ", i)

        
        s0 = expt.beam.get_s0()
        sel = refls['id'] == i
        refls_sel = refls.select(sel)
        xyzobses = refls_sel['xyzobs.px.value']
        intensities = refls_sel['intensity.sum.value']
        panels = refls_sel['panel']
        shoeboxes = refls_sel['shoebox']

        for i_refl in range(len(refls_sel)):
            i_panel = panels[i_refl]
            panel = expt.detector[i_panel]
            sb = shoeboxes[i_refl]
            sbpixels = zip(sb.coords(), sb.values())

#            xy = xyzobses[i_refl][0:2]
#            intensity = intensities[i_refl]
#            res = panel.get_resolution_at_pixel(s0, xy)
#            res_inv = 1/res
#            i_bin = int(n_bins * (res_inv - d_inv_low) / (d_inv_high - d_inv_low))
#            if i_bin < 0 or i_bin >= n_bins: continue
#            sums[i_bin] += intensity
            for (x,y,_), intensity in sbpixels:
                res = panel.get_resolution_at_pixel(s0, (x,y))
                res_inv = 1/res
                i_bin = int(n_bins * (res_inv - d_inv_low) / (d_inv_high - d_inv_low))
                if i_bin < 0 or i_bin >= n_bins: continue
                sums[i_bin] += intensity

                
#        for i_refl in range(len(refls_sel)):
#            xy = xyzobses[i_refl][0:2]
#            i_panel = panels[i_refl]
#            intensity = intensities[i_refl]
#            panel = expt.detector[i_panel]
#            res = panel.get_resolution_at_pixel(s0, xy)
#            i_bin = int(n_bins * (res - d_max) / (d_min - d_max))
#            if i_bin < 0 or i_bin >= n_bins: continue
#            sums[i_bin] += intensity

    xvalues = np.linspace(d_inv_low, d_inv_high, n_bins)
    yvalues = np.array(sums)
    data = np.concatenate((xvalues, yvalues))
    np.save('out', data)
    fig, ax = plt.subplots()
    plt.plot(xvalues, yvalues)
    ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
        lambda x, _: "{:.3f}".format(1/x)))
#    plt.xlim(d_max, d_min)
    plt.show()


if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()

