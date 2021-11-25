# Replication Files for Park and Kang (2021)

This Github repository contains replication files for [Park and Kang (2021)](https://arxiv.org/abs/2110.07740 "ClusterEff").


## Files

ClutserFtSource.R and MySL.R contain functions that are used in the simulation and the data analysis.

## Folders

### Simulation

Simulation folder replicates Table 5.1 and 5.2 in the paper and here are some details of files and folders.

* 1.GenData.R generates the simulation data sets which are saved in GenData folder. The data sets are titled as Data_N[A]\_[B]PS\_B[C].csv where [A] indicates the number of cluster ([A]=25,50,100,250,500), [B] indicates the propensity model ([B]=No,Strong), and [C] indicates the simulation replicates ([C]=1,...,210).

* 2.OurMethod.R implements our method and summarizes the result in Ours folder. The result files are titled as N[A]\_[B]PS\_G[C].csv. [A] indicates the number of cluster ([A]=25,50,100,250,500), [B] indicates the propensity model ([B]=No,Strong), and [C] indicates the magnitude of intra-cluster correlation ([C]=1 and [C]=2 indicates ![s1](https://latex.codecogs.com/svg.image?\sigma_U=0.5 "\sigma_U=0.5") and ![s2](https://latex.codecogs.com/svg.image?\sigma_U=1.5 "\sigma_U=1.5")).
	* Parallel computing is strongly recommended.

* 3.OtherMethod.R implements competing methods and summarizes the result in Others folder. The result files are titled as Others_N0500\_sV\_[A]\_B[B].csv where [A] indicates the propensity model ([A]=No,Strong) and [B] indicates the simulation replicates ([B]=1,...,210).
	* Parallel computing is strongly recommended.

* 4.Summary.R summarizes the results in Ours and Others folders.

### Simulation\_ICC

Simulation\_ICC replicates Figure A.1 in the supplementary material. 1.Illustration.R generates Illustration_Grid.csv summarizing the simulation results and draws Figure A.1.

### Simulation\_NonNormal

Simulation\_NonNormal replicates Table 5.3 in the paper and here are some details of files and folders.

* 1.GenData.R generates the simulation data sets which are saved in GenData folder. The data sets are titled as Data_N500\_[A]PS\_B[B].csv where [A] indicates the propensity model ([A]=No,Strong) and [B] indicates the simulation replicates ([C]=1,...,200).

* 2.OurMethod.R implements our method and summarizes the result in Ours folder. The result files are titled as N0500\_[A]PS\_G1.csv where [A] indicates the propensity model ([A]=No,Strong).
	* Parallel computing is strongly recommended.

* 3.OtherMethod.R implements competing methods and summarizes the result in Others folder. The result files are titled as Others_N0500\_sV\_[A]\_B[B].csv where [A] indicates the propensity model ([A]=No,Strong) and [B] indicates the simulation replicates ([B]=1,...,200).
	* Parallel computing is strongly recommended.

* 4.Summary.R summarizes the results in Ours and Others folders.

### Data\_ECLSK

Data_ECLSK replicates Section 6.1 and A.3 of the paper and here are some details of files and folders.

* Download data from https://nces.ed.gov/ecls/dataproducts.asp and make ECLSK\_Kto8\_child\_STATA.dta.

* By running 1.ECLSK.R, the data analysis in Section 6.1 of the paper is replicated. In particular,
	* 1.ECLSK.R cleans the raw data and makes Reading1\_ECLSK.csv.
	* 1.ECLSK.R implements our method and saves the files in Reading\_u folder titled Reading1\_[A]\_NF[B].csv. Here [A] indicates the nuisance functions and main/auxiliary samples in cross-fitting procedures ([A]=CPS,IPS,OR,SSlist) and [B] indicates the sample split replicates ([B]=1,...,100). The results in Reading\_u folder are summarized as Reading\_Ours.csv. 
	* 1.ECLSK.R implements competing methods and summarizes as Reading\_RU.csv and Reading\_GRF.csv.
	* Table 6.1, Table A.2, and Figure A.2 are generated.

### Data\_HIV

Data_HIV replicates Section 6.2 in the paper and here are some details of files and folders.

* Download data from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/CVOPZL and use bio\_for\_analysis.dta.

* By running 1.HIV.R, the data analysis in Section 6.2 of the paper is replicated. In particular,
	* 1.HIV.R cleans the raw data and makes Data\_Cleaned.csv.
	* 1.HIV.R implements our method and saves the files in Outcome folder titled Dupas\_OR2\_NF[A].csv and SSlist\_NF[A].csv. Here [A] indicates the sample split replicates ([A]=1,...,100). The results in Outcome folder are summarized as HIV\_Ours.csv. 
	* 1.HIV.R implements competing methods and summarizes as HIV\_RU.csv and HIV\_GRF.csv.
	* Table 6.2 is generated.

## References
Chan Park & Hyunseung Kang (2021) **More Efficient, Doubly Robust, Nonparametric Estimators of Treatment Effects in Multilevel Studies**, _arXiv:2110.07740_ [[link](https://arxiv.org/abs/2110.07740 "ClusterEff")]