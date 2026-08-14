# Partitioning-Based Least Squares Regression Methods

The package `lspartition` implements tuning parameter selection, point estimation,
and pointwise and uniform inference for partitioning-based least squares series
estimators, including B-splines, compactly supported wavelets, and piecewise
polynomials.

- `lsprobust`: point estimation of regression functions and their derivatives, with robust bias-corrected pointwise and uniform inference.
- `lspkselect`: data-driven selection of the IMSE-optimal number of partitioning knots.
- `lsprobust.plot`: graphical presentation of estimates, pointwise confidence intervals, and uniform confidence bands.
- `lsplincom`: estimation and robust bias-corrected inference for user-specified linear combinations of regression functions across groups.


## R Implementation

To install/update in R type:
```r
install.packages('lspartition')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/lspartition/lspartition.pdf), [CRAN repository](https://cran.r-project.org/package=lspartition).

- Examples/data: [lspartition illustration](R/lspartition_illustration.R), [bike-sharing data](R/bikeSharing.csv).


## References

For overviews and introductions, see the [nppackages website](https://nppackages.github.io).

### Software and Implementation

- Cattaneo, Farrell and Feng (2020): [lspartition: Partitioning-Based Least Squares Regression](https://nppackages.github.io/references/Cattaneo-Farrell-Feng_2020_R.pdf). _R Journal_ 12(1): 172-187.

### Technical and Methodological

- Cattaneo and Farrell (2013): [Optimal Convergence Rates, Bahadur Representation, and Asymptotic Normality of Partitioning Estimators](https://nppackages.github.io/references/Cattaneo-Farrell_2013_JoE.pdf). _Journal of Econometrics_ 174(2): 127-143. [Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Farrell_2013_JoE--Supplemental.pdf).
- Cattaneo, Farrell and Feng (2020): [Large Sample Properties of Partitioning-Based Series Estimators](https://nppackages.github.io/references/Cattaneo-Farrell-Feng_2020_AoS.pdf). _Annals of Statistics_ 48(3): 1718-1741. [Supplemental Appendix](https://nppackages.github.io/references/Cattaneo-Farrell-Feng_2020_AoS--Supplemental.pdf).

## Funding

This work was supported in part by the National Science Foundation through grant [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931).
