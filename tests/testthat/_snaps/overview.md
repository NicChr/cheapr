# overview

    Code
      overview(numeric(), hist = TRUE)
    Output
      obs: 0 
      cols: NA 
      
      ----- Numeric -----
        col n_missng p_complt n_unique mean p0 p25 p50 p75 p100 iqr sd hist
      1   x       NA       NA       NA   NA NA  NA  NA  NA   NA  NA NA     

---

    Code
      overview(numeric(), hist = FALSE)
    Output
      obs: 0 
      cols: NA 
      
      ----- Numeric -----
        col n_missng p_complt n_unique mean p0 p25 p50 p75 p100 iqr sd
      1   x       NA       NA       NA   NA NA  NA  NA  NA   NA  NA NA

---

    Code
      overview(airquality, hist = FALSE)
    Output
      obs: 153 
      cols: 6 
      
      ----- Numeric -----
            col n_missng p_complt n_unique   mean   p0    p25   p50    p75  p100
      1   Ozone       37     0.76       67  42.13    1     18  31.5  63.25   168
      2 Solar.R        7     0.95      117 185.93    7 115.75   205 258.75   334
      3    Wind        0        1       31   9.96  1.7    7.4   9.7   11.5  20.7
      4    Temp        0        1       40  77.88   56     72    79     85    97
      5   Month        0        1        5   6.99    5      6     7      8     9
      6     Day        0        1       31   15.8    1      8    16     23    31
           iqr    sd
      1  45.25 32.99
      2    143 90.06
      3    4.1  3.52
      4     13  9.47
      5      2  1.42
      6     15  8.86

---

    Code
      overview(iris, hist = FALSE)
    Output
      obs: 150 
      cols: 5 
      
      ----- Numeric -----
                 col n_missng p_complt n_unique mean  p0 p25  p50 p75 p100 iqr   sd
      1 Sepal.Length        0        1       35 5.84 4.3 5.1  5.8 6.4  7.9 1.3 0.83
      2  Sepal.Width        0        1       23 3.06   2 2.8    3 3.3  4.4 0.5 0.44
      3 Petal.Length        0        1       43 3.76   1 1.6 4.35 5.1  6.9 3.5 1.77
      4  Petal.Width        0        1       22  1.2 0.1 0.3  1.3 1.8  2.5 1.5 0.76
      
      ----- Categorical -----
            col n_missng p_complt n_unique n_levels    min       max
      1 Species        0        1        3        3 setosa virginica

---

    Code
      overview(iris2, hist = FALSE)
    Output
      obs: 100 
      cols: 7 
      
      ----- Logical -----
          col n_missng p_complt n_true n_false p_true
      1 large        0        1     24      76   0.24
      
      ----- Numeric -----
                 col n_missng p_complt n_unique mean  p0 p25  p50  p75 p100  iqr   sd
      1 Sepal.Length        0        1       28 5.47 4.3   5  5.4  5.9    7  0.9 0.64
      2  Sepal.Width        0        1       23  3.1   2 2.8 3.05  3.4  4.4  0.6 0.48
      3 Petal.Length        0        1       28 2.86   1 1.5 2.45 4.32  5.1 2.83 1.45
      4  Petal.Width        0        1       15 0.79 0.1 0.2  0.8  1.3  1.8  1.1 0.57
      
      ----- Categorical -----
             col n_missng p_complt n_unique n_levels    min        max
      1  Species        0        1        2        3 setosa versicolor
      2 Species2        0        1        2       NA setosa versicolor

---

    Code
      overview(warpbreaks, hist = FALSE)
    Output
      obs: 54 
      cols: 3 
      
      ----- Numeric -----
           col n_missng p_complt n_unique  mean p0   p25 p50 p75 p100   iqr   sd
      1 breaks        0        1       31 28.15 10 18.25  26  34   70 15.75 13.2
      
      ----- Categorical -----
            col n_missng p_complt n_unique n_levels min max
      1    wool        0        1        2        2   A   B
      2 tension        0        1        3        3   L   H

---

    Code
      overview(ToothGrowth, hist = FALSE)
    Output
      obs: 60 
      cols: 3 
      
      ----- Numeric -----
         col n_missng p_complt n_unique  mean  p0   p25   p50   p75 p100  iqr   sd
      1  len        0        1       43 18.81 4.2 13.07 19.25 25.27 33.9 12.2 7.65
      2 dose        0        1        3  1.17 0.5   0.5     1     2    2  1.5 0.63
      
      ----- Categorical -----
         col n_missng p_complt n_unique n_levels min max
      1 supp        0        1        2        2  OJ  VC

---

    Code
      overview(df)
    Output
      obs: 25 
      cols: 3 
      
      ----- Time-Series -----
               col n_missng p_complt n_unique  mean    p0   p25  p50  p75 p100  iqr
      1          y        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
      2          x        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
      3 z_Series 1        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
      4 z_Series 2        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
      5 z_Series 3        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
      6 z_Series 4        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
      7 z_Series 5        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59 1.26
          sd  hist grwth_rt
      1 1.15 ▂▃▅▆▇  -14.38%
      2 1.15 ▂▃▅▆▇  -14.38%
      3 1.15 ▂▃▅▆▇  -14.38%
      4 1.15 ▂▃▅▆▇  -14.38%
      5 1.15 ▂▃▅▆▇  -14.38%
      6 1.15 ▂▃▅▆▇  -14.38%
      7 1.15 ▂▃▅▆▇  -14.38%

---

    Code
      overview(ts(matrix(x, ncol = 5)))
    Output
      obs: 5 
      cols: 5 
      
      ----- Time-Series -----
             col n_missng p_complt n_unique  mean    p0   p25    p50  p75 p100  iqr
      1 Series 1        0        1        5   0.2 -1.24 -0.36   0.66 0.82 1.14 1.18
      2 Series 2        0        1        5  0.76 -0.39  0.77   0.81  1.3 1.33 0.53
      3 Series 3        0        1        5 -0.67 -2.53  -1.6  -0.11 0.24 0.64 1.84
      4 Series 4        0        1        5  0.21  -0.7 -0.44  0.076 0.53 1.59 0.97
      5 Series 5        0        1        5 -0.24 -2.36 -1.23  0.028 0.86 1.51 2.09
          sd  hist grwth_rt
      1 0.98 ▅▅▁▅▇   -8.04%
      2  0.7 ▅▁▁▇▇  -12.77%
      3 1.34 ▅▅▁▅▇   97.37%
      4 0.91 ▇▅▅▁▅      NA%
      5 1.57 ▅▅▁▅▇      NA%

---

    Code
      overview(EuStockMarkets)
    Output
      obs: 1860 
      cols: 4 
      
      ----- Time-Series -----
         col n_missng p_complt n_unique    mean      p0     p25     p50     p75
      1  DAX        0        1     1774 2530.66 1402.34  1744.1 2140.56 2722.37
      2  SMI        0        1     1725 3376.22  1587.4 2165.62 2796.35 3812.43
      3  CAC        0        1     1617 2227.83    1611 1875.15  1992.3 2274.35
      4 FTSE        0        1     1729 3565.64    2281 2843.15  3246.6 3993.57
           p100     iqr      sd  hist grwth_rt
      1 6186.09  978.26 1084.79 ▇▂▂▁▁   0.065%
      2    8412  1646.8 1663.03 ▇▃▁▁▁   0.082%
      3  4388.5   399.2  580.31 ▇▂▁▁▁   0.044%
      4    6179 1150.43  976.72 ▇▇▂▂▂   0.043%

