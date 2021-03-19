#!/usr/bin/env python3.9
import click

from src.time_series import time_series
from src.bifurcations import run_bifurcations
from src.phase_plot import phase_plot

time_serie = time_series()

@click.command()
@click.option("--plot", default="time-serie", help="Type of plot")
def main(plot):
  if plot == "time-series":
    return time_serie.plot()
  elif plot == "bifurcations":
    return run_bifurcations()
  elif plot == "phase-plot":
    return phase_plot()
  else:
    print("command not found.")

if __name__ == "__main__":
  main()
