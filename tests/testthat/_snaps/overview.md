# overview

    Code
      overview(airquality, hist = FALSE)
    Output
      obs: 153 
      cols: 6 
      
      ----- Numeric -----
            col  class n_missng p_complt n_unique   mean   p0    p25   p50    p75
      1   Ozone integr       37     0.76       67  42.13    1     18  31.5  63.25
      2 Solar.R integr        7     0.95      117 185.93    7 115.75   205 258.75
      3    Wind numerc        0        1       31   9.96  1.7    7.4   9.7   11.5
      4    Temp integr        0        1       40  77.88   56     72    79     85
      5   Month integr        0        1        5   6.99    5      6     7      8
      6     Day integr        0        1       31   15.8    1      8    16     23
         p100    iqr    sd
      1   168  45.25 32.99
      2   334    143 90.06
      3  20.7    4.1  3.52
      4    97     13  9.47
      5     9      2  1.42
      6    31     15  8.86

---

    Code
      overview(iris, hist = FALSE)
    Output
      obs: 150 
      cols: 5 
      
      ----- Numeric -----
                 col  class n_missng p_complt n_unique mean  p0 p25  p50 p75 p100 iqr
      1 Sepal.Length numerc        0        1       35 5.84 4.3 5.1  5.8 6.4  7.9 1.3
      2  Sepal.Width numerc        0        1       23 3.06   2 2.8    3 3.3  4.4 0.5
      3 Petal.Length numerc        0        1       43 3.76   1 1.6 4.35 5.1  6.9 3.5
      4  Petal.Width numerc        0        1       22  1.2 0.1 0.3  1.3 1.8  2.5 1.5
          sd
      1 0.83
      2 0.44
      3 1.77
      4 0.76
      
      ----- Categorical -----
            col  class n_missng p_complt n_unique n_levels    min       max
      1 Species factor        0        1        3        3 setosa virginica

---

    Code
      overview(iris2, hist = FALSE)
    Output
      obs: 100 
      cols: 7 
      
      ----- Logical -----
          col  class n_missng p_complt n_true n_false p_true
      1 large logicl        0        1     24      76   0.24
      
      ----- Numeric -----
                 col  class n_missng p_complt n_unique mean  p0 p25  p50  p75 p100
      1 Sepal.Length numerc        0        1       28 5.47 4.3   5  5.4  5.9    7
      2  Sepal.Width numerc        0        1       23  3.1   2 2.8 3.05  3.4  4.4
      3 Petal.Length numerc        0        1       28 2.86   1 1.5 2.45 4.32  5.1
      4  Petal.Width numerc        0        1       15 0.79 0.1 0.2  0.8  1.3  1.8
         iqr   sd
      1  0.9 0.64
      2  0.6 0.48
      3 2.83 1.45
      4  1.1 0.57
      
      ----- Categorical -----
             col  class n_missng p_complt n_unique n_levels    min        max
      1  Species factor        0        1        2        3 setosa versicolor
      2 Species2 chrctr        0        1        2       NA setosa versicolor

---

    Code
      overview(warpbreaks, hist = FALSE)
    Output
      obs: 54 
      cols: 3 
      
      ----- Numeric -----
           col  class n_missng p_complt n_unique  mean p0   p25 p50 p75 p100   iqr
      1 breaks numerc        0        1       31 28.15 10 18.25  26  34   70 15.75
          sd
      1 13.2
      
      ----- Categorical -----
            col  class n_missng p_complt n_unique n_levels min max
      1    wool factor        0        1        2        2   A   B
      2 tension factor        0        1        3        3   L   H

---

    Code
      overview(ToothGrowth, hist = FALSE)
    Output
      obs: 60 
      cols: 3 
      
      ----- Numeric -----
         col  class n_missng p_complt n_unique  mean  p0   p25   p50   p75 p100  iqr
      1  len numerc        0        1       43 18.81 4.2 13.07 19.25 25.27 33.9 12.2
      2 dose numerc        0        1        3  1.17 0.5   0.5     1     2    2  1.5
          sd
      1 7.65
      2 0.63
      
      ----- Categorical -----
         col  class n_missng p_complt n_unique n_levels min max
      1 supp factor        0        1        2        2  OJ  VC

---

    Code
      overview(df)
    Output
      obs: 25 
      cols: 3 
      
      ----- Time-Series -----
               col class n_missng p_complt n_unique  mean    p0   p25  p50  p75 p100
      1          y    ts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
      2          x    ts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
      3 z_Series 1   mts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
      4 z_Series 2   mts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
      5 z_Series 3   mts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
      6 z_Series 4   mts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
      7 z_Series 5   mts        0        1       25 0.053 -2.53 -0.44 0.24 0.82 1.59
         iqr   sd  hist
      1 1.26 1.15 ▂▃▅▆▇
      2 1.26 1.15 ▂▃▅▆▇
      3 1.26 1.15 ▂▃▅▆▇
      4 1.26 1.15 ▂▃▅▆▇
      5 1.26 1.15 ▂▃▅▆▇
      6 1.26 1.15 ▂▃▅▆▇
      7 1.26 1.15 ▂▃▅▆▇

---

    Code
      overview(ts(matrix(x, ncol = 5)))
    Output
      obs: 5 
      cols: 5 
      
      ----- Time-Series -----
             col class n_missng p_complt n_unique  mean    p0   p25    p50  p75 p100
      1 Series 1   mts        0        1        5   0.2 -1.24 -0.36   0.66 0.82 1.14
      2 Series 2   mts        0        1        5  0.76 -0.39  0.77   0.81  1.3 1.33
      3 Series 3   mts        0        1        5 -0.67 -2.53  -1.6  -0.11 0.24 0.64
      4 Series 4   mts        0        1        5  0.21  -0.7 -0.44  0.076 0.53 1.59
      5 Series 5   mts        0        1        5 -0.24 -2.36 -1.23  0.028 0.86 1.51
         iqr   sd  hist
      1 1.18 0.98 ▅▅▁▅▇
      2 0.53  0.7 ▅▁▁▇▇
      3 1.84 1.34 ▅▅▁▅▇
      4 0.97 0.91 ▇▅▅▁▅
      5 2.09 1.57 ▅▅▁▅▇

---

    Code
      overview(EuStockMarkets)
    Output
      obs: 1860 
      cols: 4 
      
      ----- Time-Series -----
         col class n_missng p_complt n_unique    mean      p0     p25     p50     p75
      1  DAX   mts        0        1     1774 2530.66 1402.34  1744.1 2140.56 2722.37
      2  SMI   mts        0        1     1725 3376.22  1587.4 2165.62 2796.35 3812.43
      3  CAC   mts        0        1     1617 2227.83    1611 1875.15  1992.3 2274.35
      4 FTSE   mts        0        1     1729 3565.64    2281 2843.15  3246.6 3993.57
           p100     iqr      sd  hist
      1 6186.09  978.26 1084.79 ▇▂▂▁▁
      2    8412  1646.8 1663.03 ▇▃▁▁▁
      3  4388.5   399.2  580.31 ▇▂▁▁▁
      4    6179 1150.43  976.72 ▇▇▂▂▂

