#!/usr/bin/env python3.9
import click

from src.time_series import time_series
from src.bifurcations import bifurcations

time_serie = time_series()

@click.command()
@click.option("--plot", default="time-serie", help="Type of plot")
def main(plot):
  if plot == "time-serie":
    return time_serie.plot()
  if plot == "bifurcation":
    return bifurcations()

if __name__ == "__main__":
  main()
