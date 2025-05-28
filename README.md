# Best Subset Selection: Optimal Pursuit for Feature Selection and Elimination


This set of Matlab (Vision R2021b) functions contain the core code to reproduce the results of the above paper.

arXiv: https://arxiv.org/abs/2501.16815.
***

## Main Fuctions 
 The following three main functions are located in the main folder. When running the **DEMO folders**, ensure that these three main functions **are added to the path** of the demo subfolders.
* `op.m` ---> the main function of the Optimal Pursuit (OP) algorithm.

* `Cosaop.m` ---> the main function of the Compressive Sampling Optimal Pursuit (CoSaOP) algorithm.

* `op_bess.m` ---> the main function of the Optimal Pursuit-Best Subset Selection (OP-BESS) algorithm.

***


## Folder 1: DEMO_Synthetic
* `demo_synthetic_CS.m`   ---> generate a compressed sensing demo for synthetic sparse data using three enhanced algorithms.

***

## Folder 2: DEMO_Audio

* `demo_Audio_CS.m`    ---> generate a compressed sensing demo for AudioSet from [1].
    
* _0bN5mYLXb0.mav [1]    ---> An audio in WAV format for testing purposes are also provided here.
  
***
  
## Folder 3: DEMO_Sparse_Regression
* `demo_Boston_housing.m`   ---> generate a sparse regression demo for Boston Housing dataset.
* data_boston_X.mat, data_boston_y.mat        ---> The provided Boston housing dataset [2] augments the feature space by constructing polynomial features.
* data_superconductivty_X_y.mat               ---> The superconductivty dataset [3]. The last column of this dataset represents the response vector y, while the remaining columns is the design matrix X.



##




**For bug reports, please contact me at email: zhihanzhu@buaa.edu.cn.**


Authors: Zhihan Zhu, Yanhao Zhang, Yong Xia.

Beihang University,  May, 26, 2025.



##

# References

[1] Gemmeke, Jort F., et al. "Audio set: An ontology and human-labeled dataset for audio events." 2017 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). IEEE, 2017.

[2] Pedregosa, Fabian, et al. "Scikit-learn: Machine learning in Python." the Journal of Machine Learning Research 12 (2011): 2825-2830.

[3] Hamidieh, Kam. "Superconductivty Data." UCI Machine Learning Repository, 2018, https://doi.org/10.24432/C53P47.