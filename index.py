#!/usr/bin/env python3.9
import click

from src.time_series import time_series
#from src.bifurcations import bifurcations

time_serie = time_series()

@click.command()
@click.option("--plot", default="time-serie", help="Type of plot")
def main(plot):
  if plot == "time-serie":
    return time_serie.basic()
  if plot == "time-serie-noise":
    return time_serie.stochastic()
  if plot == "bifurcation":
    pass

if __name__ == "__main__":
  main()
