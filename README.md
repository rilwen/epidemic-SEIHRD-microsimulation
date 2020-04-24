# epidemic-SEIHRD-microsimulation
Python implementation of a simple microsimulation with compartmental epidemic dynamics (SEIHRD)  (C) Averisera Ltd 2020.

The model describes a simplistic microsimulationof a population undergoing the spread of a hypothetical contagion described by a compartmental (SEIHRD) model.
Each person in the population has a family of arbitrary size and interacts with arbitrary frequency with arbitrary number of randomly selected extrafamilial persons. The contagion spreads through contacts between individuals, who undergo transitions between different stages of the disease as described by the SEIHRD model: S(usceptible), E(exposed), I(nfected), R(ecovered) or (D)ead.
The model can be used to test the effects of varying the number and frequency of external contacts on the disease impact (e.g. number of deaths or hospitalised) as a simplistic simulation of the effects of self-isolation policy interventions.

Licensed under GPL v3.
