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

import pickle

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
  peak_position = xyzobs
    .type = str
  peak_weighting = unit
    .type = str
  split_detectors = False
    .type = bool
output {
  log = dials.powder_from_spots.log
    .type = str
  d_table = d_table.pkl
    .type = str
  xy_file = None
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

    #sums = flex.double(params.n_bins)

    sums0 = flex.double(params.n_bins)
    sums1 = flex.double(params.n_bins)
    sums2 = flex.double(params.n_bins)
    sums3 = flex.double(params.n_bins)
    sums4 = flex.double(params.n_bins)
    sums5 = flex.double(params.n_bins)
    sums6 = flex.double(params.n_bins)
    sums7 = flex.double(params.n_bins)
    panelsums = {
        0: sums0,
        1: sums1,
        2: sums2,
        3: sums3,
        4: sums4,
        5: sums5,
        6: sums6,
        7: sums7,
        }
    d_table = []

    refls = params.input.reflections[0].data
    expts = params.input.experiments[0].data

    import random
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
            #if i_panel not in [1,2,5,6]: continue
            panel = expt.detector[i_panel]
            sb = shoeboxes[i_refl]
            sbpixels = zip(sb.coords(), sb.values())

            
            xy = xyzobses[i_refl][0:2]
            intensity = intensities[i_refl]
            res = panel.get_resolution_at_pixel(s0, xy)
            d_table.append((res, intensity))
            if params.peak_position=="xyzobs":
                res_inv = 1/res
                i_bin = int(n_bins * (res_inv - d_inv_low) / (d_inv_high - d_inv_low))
                if i_bin < 0 or i_bin >= n_bins: continue
                panelsums[i_panel][i_bin] += intensity if params.peak_weighting=="intensity" else 1
            if params.peak_position=="shoebox":
                for (x,y,_), value in sbpixels:
                    res = panel.get_resolution_at_pixel(s0, (x,y))
                    res_inv = 1/res
                    i_bin = int(n_bins * (res_inv - d_inv_low) / (d_inv_high - d_inv_low))
                    if i_bin < 0 or i_bin >= n_bins: continue
                    panelsums[i_panel][i_bin] += value if params.peak_weighting=="intensity" else 1

                

    xvalues = np.linspace(d_inv_low, d_inv_high, n_bins)
    fig, ax = plt.subplots()
    if params.split_detectors:
        offset = max(np.array(sums1))
        for i_sums, sums in enumerate([sums1, sums2, sums5, sums6]):
            yvalues = np.array(sums)
            plt.plot(xvalues, yvalues+0.5*i_sums*offset)
    else:
        yvalues = sum([v for v in panelsums.values()])
        plt.plot(xvalues, yvalues)
    ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
        lambda x, _: "{:.3f}".format(1/x)))

    if params.output.xy_file:
        with open(params.output.xy_file, 'w') as f:
            for x,y in zip(xvalues, yvalues):
                f.write("{:.6f}\t{}\n".format(1/x, y))
    plt.show()

    with open(params.output.d_table, 'wb') as f:
        pickle.dump(d_table, f)

if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()

