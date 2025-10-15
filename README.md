![Build check](https://github.com/morinlab/GAMBLR.predict/actions/workflows/build_check.yaml/badge.svg)
![GitHub R package version](https://img.shields.io/github/r-package/v/morinlab/GAMBLR.predict)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/morinlab/GAMBLR.predict)
![GitHub last commit](https://img.shields.io/github/last-commit/morinlab/GAMBLR.predict)


<p align="center" width="100%">
    <img width="33%" src="GAMBLR.predict.png"> 
</p>

This repository is part of the family of packages Genomic Analysis of Mature B-cell Lymphomas in R ([GAMBLR](https://github.com/morinlab/GAMBLR.open)) developed and maintained by the Morin Lab. This repository contains functionality for molecular/genetic classification of B-cell lymphomas (e.g. BL, FL, and DLBCL). Please refer to the main package for  more information.

## DLBCLone

This package contains the functionality to apply the DLBCLone classifier to identify genetic subgroups of Diffuse Large B-cell Lymphoma (DLBCL) described in the preprint [Klossok *et al* (2025)](https://www.medrxiv.org/content/10.1101/2025.09.18.25335809v1). The repository includes a Shiny interface with some pre-trained models. You can launch this by running `GAMBLR.predict:::DLBCLone_shiny()`. If you want to classify your own samples using an existing model, the relevant functions to get started are: `DLBCLone_load_optimized()` and `DLBCLone_predict()`.

## cFL/dFL random forest classifier

This pakage contains functionaliry to apply the random forest classifier to identify genetic subgroups of Follicular Lymphoma cFL/dFL described in [Blood (2023)](https://ashpublications.org/blood/article/142/6/561/495422/Genetic-subdivisions-of-follicular-lymphoma). The function that allows to do that is `classify_fl()`. It either accepts the metadata and mutations in maf format as input data frames to assemble the matrix and perform classification, or it can also accept the matrix prepared in advance. Please refer to the function [documentation](https://github.com/morinlab/GAMBLR.predict/blob/0f95c6b1801b9a3bee4b908eb6529c7cb4ca1551/R/classifiers.R#L1) for the examples and description of supported parameters.

## Docker

To build the docker image, run the following command in the root of the repository:

```bash
docker build -t gamblr.predict .
```

To run the docker image, use the following command:

```bash
docker run -it --rm -e PASSWORD='some_password' -p 8787:8787 gamblr.predict
```

Then open your browser and go to `http://localhost:8787`. You can log in with the username `rstudio` and the password you set in the command above.