# https://tobiasraabe.github.io/blog/matplotlib-for-publications.html
# https://www.bastibl.net/publication-quality-plots/
# https://github.com/jbmouret/matplotlib_for_papers
# http://www.jesshamrick.com/2016/04/13/reproducible-plots/
def set_size(width, fraction=1):
  """ Set aesthetic figure dimensions to avoid scaling in latex.

  Parameters
  ----------
  width: float
          Width in pts
  fraction: float
          Fraction of the width which you wish the figure to occupy

  Returns
  -------
  fig_dim: tuple
          Dimensions of figure in inches
  """
  # Width of figure
  fig_width_pt = width * fraction

  # Convert from pt to inches
  inches_per_pt = 1 / 72.27

  # Golden ratio to set aesthetic figure height
  golden_ratio = (5 ** 0.5 - 1) / 2

  # Figure width in inches
  fig_width_in = fig_width_pt * inches_per_pt
  # Figure height in inches
  fig_height_in = fig_width_in * golden_ratio

  return fig_width_in, fig_height_in
