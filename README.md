# marsh-local-adapt

Code to reproduce model fitting and simulations from the paper, **Local thermal adaptation and local temperature regimes drive the performance of parasitic helminths under climate change: the case of Marshallagia marshalli from wild ungulates** by Oscar Aleuy, Steph Peacock, Péter Molnár, Kathreen Ruckstuhl, Susan Kutz (Global Change Biology 2023).

## Abstract

Across a species’ range, populations can differ genetically in their thermal histories and, consequently, have different phenotypic responses to similar thermal stimuli. Parasitic helminths often have broad geographic ranges and temperature-sensitive life stages, but little is known about whether and how local thermal adaptation can influence their response to climate change. We studied the thermal responses of the free-living stages of Marshallagia marshalli, a parasitic nematode of wild ungulates, along a latitudinal gradient. We first determined the distribution of M. marshalli in wild sheep species in North America. Then we cultured M. marshalli eggs from different locations at temperature regimes from 5 to 38°C to determine whether they showed evidence of local thermal adaptation. **We fit thermal performance curves based on the Metabolic Theory of Ecology to time series observations of individual development and mortality. We used parameter estimates in life-cycle based host-parasite models to understand how local thermal responses may influence parasite performance under general and location-specific climate-change temperature projections.** We found that M. marshalli has a wide latitudinal and host range, infecting wild sheep species from New Mexico (33°66' N) to Yukon (64°48' N). Increases in mortality rate and development time at higher ambient temperatures were most evident for isolates from northern locations. Despite differences in development and mortality rates among M. marshalli populations, when using site-specific climate change projections, there was a similar magnitude of impact on the relative performance of M. marshalli among populations. Climate change is predicted to decrease the expected lifetime reproductive output of M. marshalli in all populations while delaying its seasonal peak by approximately one month. Our research suggests that accurate projections of the impacts of climate change on broadly distributed species need to consider local adaptations of organisms together with local temperature profiles and climate projections.

## Contents

* `model-fitting` contains JAGS model fitting code to estimate development rate and survival of Marshallagia larvae sourced from four different locations. Each parameter has files containing the JAGS models, code to fit the models to data, and code to summarize and plot the output.

* `data` contains raw data that went into the analysis, most notably the temperature experiments (`Marsh_data.csv`) and climate change projections (summaries in the `climate-change` folder; raw data too large to upload here)

* `output` contains model estimates and simulation output; much of which was too large to upload here but can be produced using the R code provided or requested from the authors.

* `manuscript_figs.R` sources output files and produces the figures in the manuscript.

* `functions.R` and `R0_function.R` contain custom R functions written to simulate the parasite model and are called in `manuscript_figs.R` and other files.

* `R0_sim.R` simulates the parasite model under different scenarios (i.e., climate and parameter combinations)

## More information

The output files shared here are not complete as some are quite large. The code is not perfectly commented and there are many auxiliary files that were not uploaded...Please don't hesitate to reach out to the author Steph Peacock <stephanie.j.peacock at gmail.com> with questions or for copies of the model output.

This research was supported, in part, by an NSERC Postdoctoral Fellowship and Banting Fellowship to Steph Peacock.
