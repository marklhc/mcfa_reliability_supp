Testing implementation in `semTools::compRelSEM()`
================

-   [Load Packages](#load-packages)
-   [Configural Model in Lai (2021)](#configural-model-in-lai-2021)
    -   [Data Import](#data-import)
    -   [Comparison](#comparison)
-   [Saturated Within-Cluster Model](#saturated-within-cluster-model)
    -   [Data Import](#data-import-1)
    -   [Comparison](#comparison-1)
-   [Simultaneous Shared-and-Configural
    Model](#simultaneous-shared-and-configural-model)
    -   [Comparison](#comparison-2)

## Load Packages

``` r
library(tidyverse)
library(lavaan)
library(semTools)
```

## Configural Model in [Lai (2021)](https://doi.org/10.1037/met0000287)

### Data Import

``` r
timss_usa <- haven::read_sav("../../timss/data/asgUSAm4.sav") %>%
    select(
        IDCNTRY, IDSCHOOL,
        AS4MAMOR, AS4MAENJ, AS4MALIK, AS4MABOR
    ) %>%
    drop_na() %>%
    mutate(AS4MABORr = 5 - AS4MABOR)
```

### Comparison

``` r
# Note: need to use the same name for the factors at level-1 and level-2
mcfa11 <- 
   'level: 1
     f1 =~ NA * AS4MAMOR + l1 * AS4MAMOR + l2 * AS4MAENJ + l3 * AS4MALIK + l4 * AS4MABORr
     AS4MAMOR ~~ ev1w * AS4MAMOR
     AS4MAENJ ~~ ev2w * AS4MAENJ
     AS4MALIK ~~ ev3w * AS4MALIK
     AS4MABORr ~~ ev4w * AS4MABORr
     f1 ~~ 1 * f1
   level: 2
     f1 =~ NA * AS4MAMOR + l1 * AS4MAMOR + l2 * AS4MAENJ + l3 * AS4MALIK + l4 * AS4MABORr
     # fixed residual variances
     AS4MAMOR ~~ ev1b * AS4MAMOR
     AS4MAENJ ~~ 0 * AS4MAENJ
     AS4MALIK ~~ ev3b * AS4MALIK
     AS4MABORr ~~ ev4b * AS4MABORr
     f1 ~~ vf1b * f1
   # tilde omega values:
   tilomgb := (l1 + l2 + l3 + l4)^2 * vf1b /
              ((l1 + l2 + l3 + l4)^2 * vf1b + ev1b + 0 + ev3b + ev4b)
   tilomgw := (l1 + l2 + l3 + l4)^2 * 1 /
              ((l1 + l2 + l3 + l4)^2 * 1 + ev1w + ev2w + ev3w + ev4w)
   # score reliability:
   omg2l := (l1 + l2 + l3 + l4)^2 * (1 + vf1b) /
            ((l1 + l2 + l3 + l4)^2 * (1 + vf1b) + 
             ev1b + 0 + ev3b + ev4b + ev1w + ev2w + ev3w + ev4w)
   omgb := (l1 + l2 + l3 + l4)^2 * vf1b /
           ((l1 + l2 + l3 + l4)^2 * vf1b + ev1b + 0 + ev3b + ev4b + 
            (ev1w + ev2w + ev3w + ev4w + (l1 + l2 + l3 + l4)^2) / 25.1)'

mcfa11_fit <- cfa(mcfa11, data = timss_usa, cluster = "IDSCHOOL")
```

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "AS4MABORr" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 245

``` r
summary(mcfa11_fit)
```

    ## lavaan 0.6-11 ended normally after 40 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        20
    ##   Number of equality constraints                     4
    ##                                                       
    ##   Number of observations                          7475
    ##   Number of clusters [IDSCHOOL]                    257
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                                40.111
    ##   Degrees of freedom                                 8
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## 
    ## Level 1 [within]:
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   f1 =~                                               
    ##     AS4MAMOR  (l1)    0.787    0.011   68.933    0.000
    ##     AS4MAENJ  (l2)    0.822    0.009   87.990    0.000
    ##     AS4MALIK  (l3)    0.881    0.010   87.944    0.000
    ##     AS4MABORr (l4)    0.748    0.011   65.527    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .AS4MAMOR          0.000                           
    ##    .AS4MAENJ          0.000                           
    ##    .AS4MALIK          0.000                           
    ##    .AS4MABORr         0.000                           
    ##     f1                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .AS4MAMO (ev1w)    0.563    0.011   51.246    0.000
    ##    .AS4MAEN (ev2w)    0.240    0.007   36.414    0.000
    ##    .AS4MALI (ev3w)    0.273    0.008   35.825    0.000
    ##    .AS4MABO (ev4w)    0.622    0.012   52.861    0.000
    ##     f1                1.000                           
    ## 
    ## 
    ## Level 2 [IDSCHOOL]:
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   f1 =~                                               
    ##     AS4MAMOR  (l1)    0.787    0.011   68.933    0.000
    ##     AS4MAENJ  (l2)    0.822    0.009   87.990    0.000
    ##     AS4MALIK  (l3)    0.881    0.010   87.944    0.000
    ##     AS4MABORr (l4)    0.748    0.011   65.527    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .AS4MAMOR          2.236    0.022  103.117    0.000
    ##    .AS4MAENJ          1.836    0.019   98.944    0.000
    ##    .AS4MALIK          1.874    0.020   93.584    0.000
    ##    .AS4MABORr         1.937    0.019  101.709    0.000
    ##     f1                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .AS4MAMO (ev1b)    0.027    0.005    5.945    0.000
    ##    .AS4MAEN           0.000                           
    ##    .AS4MALI (ev3b)    0.001    0.002    0.951    0.342
    ##    .AS4MABO (ev4b)    0.005    0.003    1.949    0.051
    ##     f1      (vf1b)    0.081    0.011    7.427    0.000
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     tilomgb           0.962    0.008  126.236    0.000
    ##     tilomgw           0.861    0.003  321.628    0.000
    ##     omg2l             0.867    0.003  328.512    0.000
    ##     omgb              0.622    0.032   19.645    0.000

``` r
coef(mcfa11_fit, type = "user")[c("tilomgw", "omg2l", "omgb")]
```

    ##   tilomgw     omg2l      omgb 
    ## 0.8605494 0.8674610 0.6221596

``` r
comp_rel <- compRelSEM(mcfa11_fit,
    obs.var = FALSE,
    config = c("f1"), shared = "f1"
)
comp_rel
```

    ## $config
    ## $config$f1
    ##   omega_W  omega_2L 
    ## 0.8605494 0.8674610 
    ## 
    ## 
    ## $shared
    ## $shared$f1
    ##   omega_B       IRR 
    ## 0.6223438 0.6467215

Discrepancy in the between reliability likely due to rounding error in
the harmonic mean of cluster sizes, in which case `compRelSEM()` should
be more accurate.

## Saturated Within-Cluster Model

### Data Import

``` r
timss15_eng <- haven::read_sav("../../timss/data/MSGUSAM3.sav") %>%
  select(
    IDCNTRY, IDSCHOOL, IDCLASS,
    MSBM18A, MSBM18B, MSBM18C, MSBM18H, MSBM18I, MSBM18M
  )
```

### Comparison

``` r
mcfas1 <- 
  'level: 1
     MSBM18A ~~ c11w * MSBM18A + c12w * MSBM18B + c13w * MSBM18C + 
                c14w * MSBM18H + c15w * MSBM18I + c16w * MSBM18M
     MSBM18B ~~ c22w * MSBM18B + c23w * MSBM18C + c24w * MSBM18H + 
                c25w * MSBM18I + c26w * MSBM18M
     MSBM18C ~~ c33w * MSBM18C + c34w * MSBM18H + c35w * MSBM18I + 
                c36w * MSBM18M
     MSBM18H ~~ c44w * MSBM18H + c45w * MSBM18I + c46w * MSBM18M
     MSBM18I ~~ c55w * MSBM18I + c56w * MSBM18M
     MSBM18M ~~ c66w * MSBM18M
   level: 2
     f1b =~ NA * MSBM18A + l1 * MSBM18A + l2 * MSBM18B + l3 * MSBM18C + 
            l4 * MSBM18H + l5 * MSBM18I + l6 * MSBM18M
     # fixed residual variances
     MSBM18A ~~ ev1b * MSBM18A
     MSBM18B ~~ ev2b * MSBM18B
     MSBM18C ~~ ev3b * MSBM18C
     MSBM18H ~~ ev4b * MSBM18H
     MSBM18I ~~ ev5b * MSBM18I
     MSBM18M ~~ ev6b * MSBM18M
     f1b ~~ 1 * f1b
   # score reliability:
   omgb := (l1 + l2 + l3 + l4 + l5 + l6)^2 /
              ((l1 + l2 + l3 + l4 + l5 + l6)^2 + 
               ev1b + ev2b + ev3b + ev4b + ev5b + ev6b + 
           (c11w + c22w + c33w + c44w + c55w + c66w + 
            2 * (c12w + c13w + c14w + c15w + c16w + 
                 c23w + c24w + c25w + c26w + 
                 c34w + c35w + c36w + c45w + c46w + c56w)) / 9.15)'
mcfas1_fit <- cfa(mcfas1, data = timss15_eng, cluster = "IDCLASS")
```

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18A" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 201 401 501 3901
    ##     4301 4601 4801 6901 9101 9201 9501 18301 18401 21501 23201 23801

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18B" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 401 901 1001 1101
    ##     1301 1501 1801 2001 3001 3301 3901 4201 4301 4801 4901 5601 5701
    ##     6401 6901 7601 8801 9101 9201 9501 9601 9901 10301 12201 12301
    ##     12501 12601 14501 15001 16301 18001 18301 18401 18601 20001 20301
    ##     21501 22601 23701 24001

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18C" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 401 501 4801 7801
    ##     9101 9201 10901 11701 15301 16301 18301 18401 23201

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18H" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 201 501 2101 3901
    ##     4301 4801 6901 9101 9201 10301 12201 13401 15301 17601 18001 18301
    ##     18601 22701

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18I" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 201 301 401 501
    ##     1301 2101 3901 4301 4601 4801 9101 9201 9501 9901 10301 11701
    ##     12201 12601 12701 15001 15201 15301 16201 16301 17101 18001 18301
    ##     18401 18601 19801

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18M" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 9101 12501 12601
    ##     20701

``` r
summary(mcfas1_fit)
```

    ## lavaan 0.6-11 ended normally after 77 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        39
    ##                                                       
    ##                                                   Used       Total
    ##   Number of observations                          2891        2954
    ##   Number of clusters [IDCLASS]                     240            
    ##                                                                   
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                                 9.146
    ##   Degrees of freedom                                 9
    ##   P-value (Chi-square)                           0.424
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## 
    ## Level 1 [within]:
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##  .MSBM18A ~~                                          
    ##    .MSBM18B (c12w)    0.212    0.009   23.085    0.000
    ##    .MSBM18C (c13w)    0.357    0.013   27.521    0.000
    ##    .MSBM18H (c14w)    0.279    0.012   24.175    0.000
    ##    .MSBM18I (c15w)    0.345    0.012   28.333    0.000
    ##    .MSBM18M (c16w)    0.310    0.014   22.738    0.000
    ##  .MSBM18B ~~                                          
    ##    .MSBM18C (c23w)    0.237    0.010   24.153    0.000
    ##    .MSBM18H (c24w)    0.191    0.009   21.576    0.000
    ##    .MSBM18I (c25w)    0.207    0.009   22.991    0.000
    ##    .MSBM18M (c26w)    0.176    0.010   17.145    0.000
    ##  .MSBM18C ~~                                          
    ##    .MSBM18H (c34w)    0.284    0.012   23.486    0.000
    ##    .MSBM18I (c35w)    0.399    0.013   30.067    0.000
    ##    .MSBM18M (c36w)    0.309    0.014   21.626    0.000
    ##  .MSBM18H ~~                                          
    ##    .MSBM18I (c45w)    0.302    0.012   26.064    0.000
    ##    .MSBM18M (c46w)    0.273    0.013   20.792    0.000
    ##  .MSBM18I ~~                                          
    ##    .MSBM18M (c56w)    0.296    0.013   22.211    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A           0.000                           
    ##    .MSBM18B           0.000                           
    ##    .MSBM18C           0.000                           
    ##    .MSBM18H           0.000                           
    ##    .MSBM18I           0.000                           
    ##    .MSBM18M           0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A (c11w)    0.533    0.015   36.327    0.000
    ##    .MSBM18B (c22w)    0.338    0.009   36.674    0.000
    ##    .MSBM18C (c33w)    0.596    0.016   36.347    0.000
    ##    .MSBM18H (c44w)    0.517    0.014   36.394    0.000
    ##    .MSBM18I (c55w)    0.514    0.014   36.427    0.000
    ##    .MSBM18M (c66w)    0.748    0.021   36.446    0.000
    ## 
    ## 
    ## Level 2 [IDCLASS]:
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   f1b =~                                              
    ##     MSBM18A   (l1)    0.344    0.023   14.855    0.000
    ##     MSBM18B   (l2)    0.183    0.016   11.641    0.000
    ##     MSBM18C   (l3)    0.390    0.026   14.871    0.000
    ##     MSBM18H   (l4)    0.287    0.021   13.393    0.000
    ##     MSBM18I   (l5)    0.389    0.024   15.994    0.000
    ##     MSBM18M   (l6)    0.268    0.030    9.019    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A           1.602    0.027   59.919    0.000
    ##    .MSBM18B           1.300    0.017   78.097    0.000
    ##    .MSBM18C           1.705    0.030   56.320    0.000
    ##    .MSBM18H           1.526    0.024   64.272    0.000
    ##    .MSBM18I           1.530    0.029   52.428    0.000
    ##    .MSBM18M           1.959    0.029   66.908    0.000
    ##     f1b               0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A (ev1b)    0.002    0.002    0.967    0.334
    ##    .MSBM18B (ev2b)    0.002    0.002    1.345    0.179
    ##    .MSBM18C (ev3b)    0.009    0.003    2.881    0.004
    ##    .MSBM18H (ev4b)    0.004    0.003    1.475    0.140
    ##    .MSBM18I (ev5b)    0.002    0.002    0.983    0.326
    ##    .MSBM18M (ev6b)    0.061    0.010    6.170    0.000
    ##     f1b               1.000                           
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     omgb              0.719    0.026   27.346    0.000

``` r
coef(mcfas1_fit, type = "user")["omgb"]
```

    ##      omgb 
    ## 0.7192835

``` r
compRelSEM(mcfas1_fit,
    obs.var = FALSE,
    shared = "f1b"
)
```

    ## $shared
    ## $shared$f1b
    ##   omega_B       IRR 
    ## 0.9768044 1.0000000

The documentation of `compRelSEM()` recommends not using a saturated
within-cluster model. Using such a model, the computation assumes the
within-variance is zero, which overestimates the reliability.

## Simultaneous Shared-and-Configural Model

### Comparison

``` r
mcfa11shared <- 
  'level: 1
     f1 =~ NA * MSBM18A + l1 * MSBM18A + l2 * MSBM18B + l3 * MSBM18C + 
            l4 * MSBM18H + l5 * MSBM18I + l6 * MSBM18M
     MSBM18A ~~ ev1w * MSBM18A
     MSBM18B ~~ ev2w * MSBM18B
     MSBM18C ~~ ev3w * MSBM18C
     MSBM18H ~~ ev4w * MSBM18H
     MSBM18I ~~ ev5w * MSBM18I
     MSBM18M ~~ ev6w * MSBM18M
     f1 ~~ 1 * f1
   level: 2
     f1 =~ NA * MSBM18A + l1 * MSBM18A + l2 * MSBM18B + l3 * MSBM18C + 
            l4 * MSBM18H + l5 * MSBM18I + l6 * MSBM18M
     f1bs =~ MSBM18A + l2s * MSBM18B + l3s * MSBM18C + 
             l4s * MSBM18H + l5s * MSBM18I + l6s * MSBM18M
     # fixed residual variances
     MSBM18A ~~ ev1b * MSBM18A
     MSBM18B ~~ ev2b * MSBM18B
     MSBM18C ~~ ev3b * MSBM18C
     MSBM18H ~~ ev4b * MSBM18H
     MSBM18I ~~ ev5b * MSBM18I
     MSBM18M ~~ ev6b * MSBM18M
     f1 ~~ (.05 / (1 - .05)) * f1
     f1bs ~~ v1bs * f1bs
     f1 ~~ 0 * f1bs
   # score reliability:
   omgbs := (1 + l2s + l3s + l4s + l5s + l6s)^2 * v1bs /
            ((l1 + l2 + l3 + l4 + l5 + l6)^2 * .05 / (1 - .05) + 
             (1 + l2s + l3s + l4s + l5s + l6s)^2 * v1bs + 
             ev1b + ev2b + ev3b + ev4b + ev5b + ev6b + 
             (ev1w + ev2w + ev3w + ev4w + ev5w + ev6w + 
              (l1 + l2 + l3 + l4 + l5 + l6)^2) / 9.15)'
mcfa11sh_fit <- cfa(mcfa11shared, data = timss15_eng, cluster = "IDCLASS")
```

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18A" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 201 401 501 3901
    ##     4301 4601 4801 6901 9101 9201 9501 18301 18401 21501 23201 23801

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18B" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 401 901 1001 1101
    ##     1301 1501 1801 2001 3001 3301 3901 4201 4301 4801 4901 5601 5701
    ##     6401 6901 7601 8801 9101 9201 9501 9601 9901 10301 12201 12301
    ##     12501 12601 14501 15001 16301 18001 18301 18401 18601 20001 20301
    ##     21501 22601 23701 24001

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18C" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 401 501 4801 7801
    ##     9101 9201 10901 11701 15301 16301 18301 18401 23201

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18H" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 201 501 2101 3901
    ##     4301 4801 6901 9101 9201 10301 12201 13401 15301 17601 18001 18301
    ##     18601 22701

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18I" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 201 301 401 501
    ##     1301 2101 3901 4301 4601 4801 9101 9201 9501 9901 10301 11701
    ##     12201 12601 12701 15001 15201 15301 16201 16301 17101 18001 18301
    ##     18401 18601 19801

    ## Warning in lav_data_full(data = data, group = group, cluster = cluster, : lavaan WARNING:
    ##     Level-1 variable "MSBM18M" has no variance within some clusters.
    ##     The cluster ids with zero within variance are: 9101 12501 12601
    ##     20701

``` r
summary(mcfa11sh_fit)
```

    ## lavaan 0.6-11 ended normally after 72 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        36
    ##   Number of equality constraints                     6
    ##                                                       
    ##                                                   Used       Total
    ##   Number of observations                          2891        2954
    ##   Number of clusters [IDCLASS]                     240            
    ##                                                                   
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                               134.373
    ##   Degrees of freedom                                18
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## 
    ## Level 1 [within]:
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   f1 =~                                               
    ##     MSBM18A   (l1)    0.573    0.012   45.989    0.000
    ##     MSBM18B   (l2)    0.362    0.011   33.888    0.000
    ##     MSBM18C   (l3)    0.632    0.013   49.422    0.000
    ##     MSBM18H   (l4)    0.489    0.013   37.778    0.000
    ##     MSBM18I   (l5)    0.610    0.012   52.375    0.000
    ##     MSBM18M   (l6)    0.506    0.016   31.424    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A           0.000                           
    ##    .MSBM18B           0.000                           
    ##    .MSBM18C           0.000                           
    ##    .MSBM18H           0.000                           
    ##    .MSBM18I           0.000                           
    ##    .MSBM18M           0.000                           
    ##     f1                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A (ev1w)    0.206    0.007   29.109    0.000
    ##    .MSBM18B (ev2w)    0.208    0.006   33.786    0.000
    ##    .MSBM18C (ev3w)    0.196    0.007   27.006    0.000
    ##    .MSBM18H (ev4w)    0.279    0.009   32.590    0.000
    ##    .MSBM18I (ev5w)    0.141    0.006   24.407    0.000
    ##    .MSBM18M (ev6w)    0.491    0.014   34.108    0.000
    ##     f1                1.000                           
    ## 
    ## 
    ## Level 2 [IDCLASS]:
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   f1 =~                                               
    ##     MSBM18A   (l1)    0.573    0.012   45.989    0.000
    ##     MSBM18B   (l2)    0.362    0.011   33.888    0.000
    ##     MSBM18C   (l3)    0.632    0.013   49.422    0.000
    ##     MSBM18H   (l4)    0.489    0.013   37.778    0.000
    ##     MSBM18I   (l5)    0.610    0.012   52.375    0.000
    ##     MSBM18M   (l6)    0.506    0.016   31.424    0.000
    ##   f1bs =~                                             
    ##     MSBM18A           1.000                           
    ##     MSBM18B  (l2s)    0.515    0.041   12.487    0.000
    ##     MSBM18C  (l3s)    1.160    0.063   18.454    0.000
    ##     MSBM18H  (l4s)    0.829    0.054   15.370    0.000
    ##     MSBM18I  (l5s)    1.164    0.059   19.829    0.000
    ##     MSBM18M  (l6s)    0.768    0.089    8.602    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   f1 ~~                                               
    ##     f1bs              0.000                           
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A           1.602    0.027   60.268    0.000
    ##    .MSBM18B           1.300    0.017   78.461    0.000
    ##    .MSBM18C           1.705    0.030   56.182    0.000
    ##    .MSBM18H           1.526    0.024   64.722    0.000
    ##    .MSBM18I           1.530    0.029   52.337    0.000
    ##    .MSBM18M           1.959    0.029   66.810    0.000
    ##     f1                0.000                           
    ##     f1bs              0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .MSBM18A (ev1b)    0.003    0.002    1.170    0.242
    ##    .MSBM18B (ev2b)    0.002    0.002    1.223    0.221
    ##    .MSBM18C (ev3b)    0.009    0.003    2.789    0.005
    ##    .MSBM18H (ev4b)    0.005    0.003    1.647    0.100
    ##    .MSBM18I (ev5b)    0.002    0.002    0.705    0.481
    ##    .MSBM18M (ev6b)    0.062    0.010    6.163    0.000
    ##     f1                0.053                           
    ##     f1bs    (v1bs)    0.098    0.016    6.165    0.000
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     omgbs             0.607    0.037   16.329    0.000

``` r
# Reliability for shared construct
coef(mcfa11sh_fit, type = "user")["omgbs"]
```

    ##     omgbs 
    ## 0.6065982

``` r
compRelSEM(mcfa11sh_fit,
    obs.var = FALSE,
    config = "f1",
    shared = c("f1", "f1bs")
)
```

    ## $config
    ## $config$f1
    ##   omega_W  omega_2L 
    ## 0.8687155 0.7017764 
    ## 
    ## 
    ## $shared
    ##                f1      f1bs
    ## omega_B 0.1108770 0.8253698
    ## IRR     0.7348633 1.0000000

The discrepancy is due to level-1 variances not included for the
construct reliability of `f1bs`.

``` r
par2 <- lavInspect(mcfa11sh_fit, what = "est")$IDCLASS
lam2 <- par2$lambda
phi2 <- par2$psi
theta2 <- par2$theta
sum_lam2 <- colSums(lam2)
# Common variance for f1 at between-level
comvar_f1 <- sum_lam2[1]^2 * phi2[1, 1]
# Common variance for f1bs at between-level
comvar_f1bs <- sum_lam2[2]^2 * phi2[2, 2]
c(comvar_f1, comvar_f1bs) / (sum(comvar_f1, comvar_f1bs, theta2))
```

    ##        f1      f1bs 
    ## 0.1508811 0.8253698

Accounting for the level-1 variances:

``` r
par1 <- lavInspect(mcfa11sh_fit, what = "est")$within
lam1 <- par1$lambda
phi1 <- par1$psi
theta1 <- par1$theta
sum_lam1 <- colSums(lam1)
# Common variance for f1 at within-level
comvar_f1w <- sum_lam1[1]^2 * phi1[1, 1]
# Add in the level-1 variances
c(comvar_f1, comvar_f1bs) /
  (sum(comvar_f1, comvar_f1bs, theta2) + sum(comvar_f1w, theta1) / 9.15)
```

    ##        f1      f1bs 
    ## 0.1108887 0.6065982
