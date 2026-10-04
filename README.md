# spatiotemporal-Rt

Estimating the effective reproduction number is crucial for understanding and managing
infectious disease outbreaks. For vector-borne diseases like dengue, transmission depends on
environmental and spatial conditions: temperature affects the extrinsic incubation period in
mosquitoes, altering transmission timing, while spatial proximity can lead to clusters of
transmission. We integrated a temperature-dependent (TD) generation time (GT) distribution and
a spatial decay function weighting transmission likelihood by distance into the Wallinga & Teunis
estimation framework. Simulations compared scenarios where the true generation time for
spatially adjusted disease cases were TD or temperature independent (TI), versus their
corresponding model assumptions. Daily reproduction numbers (𝑅𝑡) estimated were evaluated via
percent error against true values. We found error to be predominantly driven by variance rather
than bias, indicating that stochastic uncertainty in infector assignment was the primary driver of
inaccuracies rather than systematic model miscalibration. Assuming a spatial weighting function
assigning higher probabilities of transmission to geographically close infector-infectee pairs
(Gaussian spatial decay) percent errors measuring deviation between estimated and true 𝑅𝑡
outperformed those of alternative spatial kernels. Mis-specifying temperature-dependence under a
Gaussian spatial kernel yielded higher errors when the ground truth was TD (16–38% vs 36–58%
deviation), and similar errors when it was TI (32–52% vs 35–56% deviation), indicating limited
sensitivity to temperature dependence when it was not present in the underlying transmission
process. Application to dengue case data in Singapore showed that spatially adjusted models
produced more variable 𝑅𝑡 trajectories than time-series smoothing approaches, with Gaussian
decay yielding more stable estimates than exponential decay and the TD model producing modest
refinements, suggesting that explicitly modeling spatial heterogeneity and TD transmission
dynamics increases responsiveness to even fluctuations in case counts. Our framework is robust
across climates, shown by simulations using a broader temperature range, and suggests a useful
application in retrospective analyses.

Files:

The main function file computes the individual reproduction number within a set of 10,000 cases within one year of time. Temperature dependence can be removed in either the ground-truth simulations or in the inference equations by removing the temperature term, T, from the corresponding likelihood functions. 

The spatial misspecified scenarios file utilizes previously-published generation time distributions for dengue transmission that are not temperature dependent, to test the effect of mis-specifying the generation time, while keeping the spatial component.

The non-spatial file contains the algorithm to reproduce the previous inference processes under the assumption that space is not a factor in viral transmission.
