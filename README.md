ActivitySpace_JRSI2014
======================

This repository contains code related to the following paper:

Perkins TA, Garcia AJ, Paz Soldan VA, Stoddard ST, Reiner RC, Vazquez-Prokopec GM, Bisanzio D, Halsey ES, Kochel TJ, Morrison AC, Smith DL, Scott TW, Tatem AJ. 2014. **Theory and data for simulating fine-scale human movement in an urban environment**. *Journal of the Royal Society Interface* 11:20140642. doi:[10.1098/rsif.2014.0642](http://rsif.royalsocietypublishing.org/content/11/99/20140642.full)

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/).

================
  
### Overview

Running `script.R` simulates a list of activity spaces, with one individual's activity space corresponding to a single element of the list `sim_data`. Each list contains vectors with the `types` of locations in the activity space, the `time` spent at each location as an overall proportion of one's daytime hours, and the `codes` of the locations as they appear in the input file `locations.txt`.

The file `simulateActivitySpace.R` contains a function that simulates a sequence of activity spaces for individuals identified in an input vector `part_codes`. Each element of `part_codes` references a row in the input variable `Participants`, which contains information about all individuals of interest in the simulated population. Other inputs include the `Locations` variable, which should include all possible locations that an individual could visit, and the `dist` variable, which contains the distance between each home of each individual in `Participants` and each location in `Locations`. Other inputs to `simulate_aspace()` include `aspaceSize`, `timeHome`, `where`, `propLoc`, and `freqDurn`, which are lists of fitted models. The remaining inputs `aspaceSizeModel`, `propLocModel`, `whereModel`, and `fdModel` are strings that specify which model from the corresponding lists of models are to be used in the simulation.

The files `aspaceSize.RData`, `freqDurn.RData`, `propLoc.RData`, `timeHome.RData`, and `where.RData` contain fitted model objects that serve as inputs to `simulate_aspace()`. Throughout, for models where different groupings of location types are considered, models are named `LN` where N is a number 1-8 for the number of location type groupings. For example, `L1` means that the parameters for all location types are the same, and `L8` means that each location type has a different set of parameters.
