# Replication Files for Park and Kang (2020)

This Github repository contains replication files for [Park and Kang (2020)](https://arxiv.org/abs/2004.08950 "EffNet").


## Folders

### Main_Data

Main_Data folder contains Estimate folder, MainCode.R, MainFunction.R, and Summary.RData files. 

* MainCode.R file replicates the data analysis in Section 5 of the paper. 
	* The raw data is available at https://www.openicpsr.org/openicpsr/project/113783/version/V1/view.

* MainFunction.R file contains functions used in MainCode.R file.

* Estimate folder contains csv files titled as RESULT_Disc_B####_T#.csv. Here B#### represents the index of sample splitting procedure ranging from B0001 to B0100 and T# represents the cluster type ranging from T1 to T5; see Section 5 of the paper.

* Summary.RData file contains the IPW/nonparametric doubly robust estimates, their standard errors, and the corresponding Wald statistics across 2-dimensional grids of the counterfactual parameters.


### Supplementary

Supplementary folder contains Results15 folder, Function15.R, Function16.R, Section15.R, Section16.R, RESULT16.csv

* Section15.R replicates the simulation analysis in Section A.5 of the paper.

* Function15.R contains functions used in Section15.R file.

* Results15 folder contains csv files titled as RESULT_pA###_al###.csv. Here pA### represents the probability ![p_A^*](https://latex.codecogs.com/svg.image?p_A^* "p_A^*") taking ![p_A^*=0.3](https://latex.codecogs.com/svg.image?p_A^*=0.3 "p_A^*=0.3"),  ![p_A^*=0.5](https://latex.codecogs.com/svg.image?p_A^*=0.5 "p_A^*=0.5") , and  ![p_A^*=0.7](https://latex.codecogs.com/svg.image?p_A^*=0.7 "p_A^*=0.7"), and al### represents the counterfactual parameter ![\alpha](https://latex.codecogs.com/svg.image?\alpha "\alpha") ranging from ![\alpha=0.01](https://latex.codecogs.com/svg.image?\alpha=0.01 "\alpha=0.01") to ![\alpha=0.99](https://latex.codecogs.com/svg.image?\alpha=0.99 "\alpha=0.99"); see Section A.5 of the paper. 

* Section16.R replicates the simulation analysis in Section A.6 of the paper.

* Function16.R contains functions used in Section16.R file.

* RESULT16.csv contains the results of the simulation analysis in Section A.6 of the paper.
<!--- where variables are named as [A]\_[B]\_[C]. Here [A] takes B, S, and C and indicates bias, standard error, and whether 95% Wald-type confidence interval covers the true parameter, respectively. [B] takes True, MisOR, MisPS, MisBoth, Over, and Under which show the model specification in Table 1 of the paper. [C] takes DE and IE which show the target parameter ![DE](https://latex.codecogs.com/svg.image?\tau^{\rm&space;DE}(\alpha) "\tau^{\rm DE}(\alpha)") and ![IE](https://latex.codecogs.com/svg.image?\tau^{\rm&space;IE}(\alpha,\alpha') "\tau^{\rm IE}(\alpha,\alpha')"), respectively. -->

## References
Chan Park & Hyunseung Kang (2020) **Efficient Semiparametric Estimation of Network Treatment Effects Under Partial Interference**, _arXiv:2004.08950_ [[link](https://arxiv.org/abs/2004.08950 "EffNet")]