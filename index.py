#!/usr/bin/env python3.9
from PyInquirer import prompt, print_json
from src.time_series import time_series
from src.bifurcations import bifurcations

if __name__ == "__main__":
  start_questions = [
    {
      "type": "list",
      "name": "tipping-points",
      "message": "Points de cascade dans le système climatique - Bifurcation de Fold-Hopf",
      "choices": [
        "Bifurcation de Fold-Hopf - série temporelle (sans bruit)",
        "Bifurcation de Fold - diagramme de bifurcation",
        "Bifurcation de Hopf - diagramme de bifurcation"
      ]
    }
  ]

  start_answers = prompt(start_questions)
  start_answer = start_answers["tipping-points"]

  if (start_answer == "Bifurcation de Fold-Hopf - série temporelle (sans bruit)"):
    time_serie = time_series()
    time_serie.rk4()

  elif (start_answer == "Bifurcation de Fold - diagramme de bifurcation"):
    bifurcation = bifurcations()
    bifurcation.fold()

  elif (start_answer == "Bifurcation de Hopf - diagramme de bifurcation"):
    bifurcation = bifurcations()
    bifurcation.hopf()
