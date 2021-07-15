---
output:
  pdf_document: default
  html_document: default
---
This file documents the reproducing guidelines for the work:

> Rui Tuo, Shiyuan He, Arash Pourhabib, Yu Ding and Jianhua Z. Huang (2021+). A Reproducing Kernel Hilbert Space Approach to Functional Calibration of Computer Models.

## Prerequisite


Reproducing the numerical results of this work requires the following R packages. 

- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)
- [snowfall](https://cran.r-project.org/web/packages/snowfall/index.html)
- [rlecuyer](https://cran.r-project.org/web/packages/rlecuyer/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)


Install the developed package for this work via

```bash
R CMD INSTALL calibration_1.0.tar.gz
```


## Steps to reproduce the results


The RMarkdown file is generated on cached data files. The cached data files can be reproduced by the following steps.

Execute the following command in the terminal. The code will run the cheap code of each simulation setting. The generated RData files can be found in the `data` folder.

```bash
Rscript ./code/simu_1.R
Rscript ./code/simu_2.R
Rscript ./code/simu_3.R
Rscript ./code/laGP/laGP_simu3.R
Rscript ./code/simu_4.R
Rscript ./code/laGP/laGP_simu4.R
```

Run the RMarkdown file `NumericalResults.Rmd` to produce the final numerical results. You will get a pdf file `NumericalResults.pdf`.

The script to apply the proposed method over the real data can be found at `./code/real_data.R`.
