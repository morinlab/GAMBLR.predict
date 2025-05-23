![Build check](https://github.com/morinlab/GAMBLR.predict/actions/workflows/build_check.yaml/badge.svg)
![GitHub R package version](https://img.shields.io/github/r-package/v/morinlab/GAMBLR.predict)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/morinlab/GAMBLR.predict)
![GitHub last commit](https://img.shields.io/github/last-commit/morinlab/GAMBLR.predict)


<p align="center" width="100%">
    <img width="33%" src="GAMBLR.predict.png"> 
</p>

This is a helper repo of the package Genomic Analysis of Mature B-cell Lymphomas in R ([GAMBLR](https://github.com/morinlab/GAMBLR.open)) developed by the Morin Lab. This repo contains the set of functions and helpers to operate on matrices for classification of B-cell lymphomas like BL, FL, and DLBCL. Please refer to the main package for  more information.

## cFL/dFL random forest classifier

This pakage contains functionaliry to apply the random forest classifier to identify genetic subgroups of Follicular Lymphoma cFL/dFL described in [Blood (2023)](https://ashpublications.org/blood/article/142/6/561/495422/Genetic-subdivisions-of-follicular-lymphoma). The function that allows to do that is `classify_fl()`. It either accepts the metadata and mutations in maf format as input data frames to assemble the matrix and perform classification, or it can also accept the matrix prepared in advance. Please refer to the function [documentation](https://github.com/morinlab/GAMBLR.predict/blob/0f95c6b1801b9a3bee4b908eb6529c7cb4ca1551/R/classifiers.R#L1) for the examples and description of supported parameters.
