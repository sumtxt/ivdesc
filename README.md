# ivdesc

Estimating the mean and variance of a covariate for the complier, never-taker and always-taker subpopulation in the context of instrumental variable estimation. This package implements the method described in [Marbach and Hangartner (2020)](https://doi.org/10.1017/pan.2019.48) and [Hangartner, Marbach, Henckel, Maathuis, Kelz and Keele 2021](https://arxiv.org/abs/2103.06328).


### Install R Package 

You can install the package directly from CRAN: 

```R
install.packages("ivdesc")
```

Or install the latest version from Github using:  

```R
remotes::install_github("sumtxt/ivdesc/R/ivdesc")
```

You may need to install the `remotes` package first. 


### Install STATA Package

To install the STATA package directly from SSC use: 

```STATA
ssc install ivdesc
```

Or install the latest version from Github using:  

```STATA
github install sumtxt/ivdesc
```

This will only work if you installed the [GITHUB module](https://github.com/haghish/github) via `net install github, from("https://haghish.github.io/github/")` first. 
