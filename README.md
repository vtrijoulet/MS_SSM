# MS_SSM
MultiSpecies State-Space Model

The model is a multispecies state-space age-structured stock assessment that can be used in estimation mode to fit to observed data or in simulation mode to simulate data. The model uses predator stomach contents. Predation can be turned off to get a multi-stock assessment model without trophic interactions between the fish species.

The model can be run via the "**MS_SSM.R**" script. This script controls for the different options that can be chosen in the model.

#### Please note:
The model depends on the installation of the TMB package.

#### References:

Trijoulet V., Fay G. and Miller T.J. (2019). Performance of a state-space multispecies model: what are the consequences of ignoring predation and process errors in stock assessments? Journal of Applied Ecology. https://doi.org/10.1111/1365-2664.13515

Trijoulet, V., Fay, G., Curti, K., Smith, B., & Miller, T. J. (2019). Performance of multispecies population models: insights on the influence of diet data. ICES Journal of Marine Science. https://doi.org/10.1093/icesjms/fsz053