# A Neural Network Approach to Latent Class Joint Models of Time-to-Event Outcomes and Non-linear Longitudinal Trajectories: An Application to Psychometric Scores in Patients with Alzheimer’s Disease
The code in this repo showcases the method proposed by Harris, Nakhaeirad and Skhosana (2026b). The main function is the `nnem.R` function.  This function is used to compute a joint latent class model for time-to-event and longitudinal data using an adapted EM algorithm under the flexible specification of nonlinear class-memberships and nonlinear longitudinal trajectories. This EM algorithm, uses a neural network to estimate the prior mixing proportions conditonal on a subjects covariates as well a network structure implemented from the `keras3` package. This approach allows users to flexibly model the mixing proportions without relying on the assumptions of a logistic regression model, which is used in the traditional approach (see Proust-Lima et. al (2017)). To showcase, the method using the data generated for simulation II, one can fit the model and plot the weighted means of the class-specific trajectories to visualize the average nonlinear trajectories per class.

Here is a display of the mixing proportions estimated as the algorithm runs for 20 iterations on simulation II for K=4 components:

![](https://github.com/tharris0924/nn-jlcm_p1/blob/main/Codes/mixing_proportions.gif)

And here are the fitted nonlinear trajectories:

![](https://github.com/tharris0924/nlme-nnjlcm/blob/main/Codes/fitted_long_traj.png)

This repo provides the code to be able to generate the data for the simulated data above in `nl_sim2.R`. The appropriate `keras3` dependancies can be installed using the requirements-r-keras.txt file with the following codes:

``
reticulate::virtualenv_create("r-keras", python = "3.10")
reticulate::virtualenv_install("r-keras", requirements = "requirements-r-keras.txt")
``


## References

Harris, T. J. E., Nakhaei Rad, N., & Skhosana, S. B. (2026). Deep learning embedded latent class joint modelling of time-to-event and longitudinal data. Statistical Methods in Medical Research, 1–39. https://doi.org/10.1177/09622802261465334 <br>

Proust-Lima, C., Philipps, V., & Liquet, B. (2017). Estimation of extended mixed models using latent classes and latent processes: The R package lcmm. Journal of Statistical Software, 78(2), 1–56. https://doi.org/10.18637/jss.v078.i02
