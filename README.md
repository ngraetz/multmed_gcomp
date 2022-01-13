# G-computation with multiple sequential mediators
Functions for g-computation to calculate mediation estimands given multiple ordered mediators. Given data wide on repeated observations of exposure, any number of mediators, and outcome, this analysis uses g-computation (or VanderWeele & Lin's 2019 mediational g-formula) to estimate 2-way effect decomposition (NDE, NIE) or 3-way effect decomposition (CDE, PAI, PIE). 

Our implementation allows for post-treatment confounders to be included in the mediator/outcome models, but we do not specifically allow for post-treatment confounders affected by the treatment, or _treatment-induced mediator-outcome confounders_.These variables may only be included as _mediators_, which we believe makes most sense for applications to complex social exposures. In other words, it is difficult to imagine a scenario with a social exposure A where a path through a post-treatment variable (A->M1->Y) should not theoretically be treated as a mediating pathway of interest with its own set of mediated effects, though M1 may still confound the relation between a future M2 and Y (A->M1->M2->Y). Still, we are working on adding an option to more flexibly handle post-treatment variables. As mediation effects will currently be estimated through every post-treatment variable provided by the user, the CDE will be different compared to implementations that may treat some of these variables as treatment-induced mediator-outcome confounders rather than mediators.

An example is provided in __run_gcomp.R__ based on __data.RData__ and __dag.csv__ 

For application and conceptual discussion, see: https://journals.sagepub.com/doi/full/10.1177/00221465211066108

This repo is a work in progress and feedback is welcome; we are in the process of converting to a much cleaner R package for CRAN. 

# Comparison to other methods for multiple sequential mediators
This repo contains a script comparing our process to the estimation of joint mediated effects through a set of multiple sequential mediators (_CMAverse_: https://bs1125.github.io/CMAverse/articles/overview.html). Our code estimates separate mediated effects through M1 and M2 (i.e., a separate NIE via both M1 and M2), whereas _CMAverse_ (and other implementations) estimates a joint mediated effect through M1/M2. In doing so, our implementation is effectively putting the path directly through both mediators (A->M1->M2->Y) into the mediated effects via M2; these separate mediated effects then add up to the joint mediated effect. This might be more or less reasonable depending on the research question. Users should decide whether this implementation (and corresponding assumptions) is aligned to the theoretical quantities relevant to their question, compared to the joint mediated effect of the set of multiple sequential mediators (_CMAverse_ package) or each individual path-specific effect (_paths_ package); see Zhou & Yamamoto (Figure 1) for excellent discussion on this topic: https://scholar.harvard.edu/files/xzhou/files/zhou-yamamoto_paths.pdf.
