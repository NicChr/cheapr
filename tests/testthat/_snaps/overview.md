# overview

    Code
      overview(airquality, hist = TRUE)
    Output
      rows: 153 cols: 6 
      
      ----- Numeric -----
            col   class n_missing p_complete n_unique   mean  p0    p25  p50    p75
      1   Ozone integer        37       0.76       67  42.13   1     18 31.5  63.25
      2 Solar.R integer         7       0.95      117 185.93   7 115.75  205 258.75
      3    Wind numeric         0          1       31   9.96 1.7    7.4  9.7   11.5
      4    Temp integer         0          1       40  77.88  56     72   79     85
      5   Month integer         0          1        5   6.99   5      6    7      8
      6     Day integer         0          1       31   15.8   1      8   16     23
        p100   iqr    sd  hist
      1  168 45.25 32.99 ▇▃▂▁▁
      2  334   143 90.06 ▅▃▅▇▅
      3 20.7   4.1  3.52 ▂▇▇▃▁
      4   97    13  9.47 ▂▃▇▇▃
      5    9     2  1.42 ▇▇▇▇▇
      6   31    15  8.86 ▇▇▇▇▆

---

    Code
      overview(iris, hist = TRUE)
    Output
      rows: 150 cols: 5 
      
      ----- Numeric -----
                 col   class n_missing p_complete n_unique mean  p0 p25  p50 p75 p100
      1 Sepal.Length numeric         0          1       35 5.84 4.3 5.1  5.8 6.4  7.9
      2  Sepal.Width numeric         0          1       23 3.06   2 2.8    3 3.3  4.4
      3 Petal.Length numeric         0          1       43 3.76   1 1.6 4.35 5.1  6.9
      4  Petal.Width numeric         0          1       22  1.2 0.1 0.3  1.3 1.8  2.5
        iqr   sd  hist
      1 1.3 0.83 ▆▇▇▅▂
      2 0.5 0.44 ▁▆▇▂▁
      3 3.5 1.77 ▇▁▆▇▂
      4 1.5 0.76 ▇▁▇▅▃
      
      ----- Categorical -----
            col  class n_missing p_complete n_unique n_levels    min       max
      5 Species factor         0          1        3        3 setosa virginica

---

    Code
      overview(iris2, hist = TRUE)
    Output
      rows: 100 cols: 7 
      
      ----- Logical -----
          col   class n_missing p_complete n_true n_false p_true
      6 large logical         0          1     24      76   0.24
      
      ----- Numeric -----
                 col   class n_missing p_complete n_unique mean  p0 p25  p50  p75
      1 Sepal.Length numeric         0          1       28 5.47 4.3   5  5.4  5.9
      2  Sepal.Width numeric         0          1       23  3.1   2 2.8 3.05  3.4
      3 Petal.Length numeric         0          1       28 2.86   1 1.5 2.45 4.32
      4  Petal.Width numeric         0          1       15 0.79 0.1 0.2  0.8  1.3
        p100  iqr   sd  hist
      1    7  0.9 0.64 ▅▇▇▃▂
      2  4.4  0.6 0.48 ▂▅▇▃▁
      3  5.1 2.83 1.45 ▇▁▁▃▅
      4  1.8  1.1 0.57 ▇▁▂▅▂
      
      ----- Categorical -----
             col     class n_missing p_complete n_unique n_levels    min        max
      5  Species    factor         0          1        2        3 setosa versicolor
      7 Species2 character         0          1        2       NA setosa versicolor

---

    Code
      overview(warpbreaks, hist = TRUE)
    Output
      rows: 54 cols: 3 
      
      ----- Numeric -----
           col   class n_missing p_complete n_unique  mean p0   p25 p50 p75 p100
      1 breaks numeric         0          1       31 28.15 10 18.25  26  34   70
          iqr   sd  hist
      1 15.75 13.2 ▇▆▃▁▁
      
      ----- Categorical -----
            col  class n_missing p_complete n_unique n_levels min max
      2    wool factor         0          1        2        2   A   B
      3 tension factor         0          1        3        3   L   H

---

    Code
      overview(ToothGrowth, hist = TRUE)
    Output
      rows: 60 cols: 3 
      
      ----- Numeric -----
         col   class n_missing p_complete n_unique  mean  p0   p25   p50   p75 p100
      1  len numeric         0          1       43 18.81 4.2 13.07 19.25 25.27 33.9
      3 dose numeric         0          1        3  1.17 0.5   0.5     1     2    2
         iqr   sd  hist
      1 12.2 7.65 ▅▅▅▇▂
      3  1.5 0.63 ▇▇▁▁▇
      
      ----- Categorical -----
         col  class n_missing p_complete n_unique n_levels min max
      2 supp factor         0          1        2        2  OJ  VC

---

    Code
      overview(EuStockMarkets, hist = TRUE)
    Output
         col class n_missing p_complete n_unique     mean      p0      p25      p50
      1  DAX array         0          1     1774 2530.657 1402.34 1744.102 2140.565
      2  SMI array         0          1     1725 3376.224 1587.40 2165.625 2796.350
      3  CAC array         0          1     1617 2227.828 1611.00 1875.150 1992.300
      4 FTSE array         0          1     1729 3565.643 2281.00 2843.150 3246.600
             p75    p100      iqr        sd  hist
      1 2722.367 6186.09  978.265 1084.7927 ▇▂▂▁▁
      2 3812.425 8412.00 1646.800 1663.0265 ▇▃▁▁▁
      3 2274.350 4388.50  399.200  580.3142 ▇▂▁▁▁
      4 3993.575 6179.00 1150.425  976.7155 ▇▇▂▂▂

