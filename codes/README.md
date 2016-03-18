存放项目中所用到的代码

## CG permutation处理程序

CG permutation根据测序数据产生随机的数据，用来作为数据分析的参照。

### mtbr2csv.R

该脚本用来产生randomCG.py所需要的输入csv文件，输入文件为mtbr的Rdata文件。

### randomCG.py

该程序用来产生permutation随机mtbr。输入文件为实际测序数据的mtbr文件(csv格式，使用mtbr2csv.R生成)以及参考序列文件。

### process_random.R

该R脚本用来批量处理和分析randomCG.py产生的随机mtbr文件，并生成结果。

## CG Hist2D处理与绘图程序

### cghist2d.lib.R

该脚本用来分析和准备CG density的绘图数据，包括Normalization，slicing等等。

### cghist2d.R

该脚本使用plotly用来绘制CG density的hist2d图。
 
## fig1cmd

该文件夹用来存放生成figure1的程序，figurePlot-hilert-distance.R生成hilbert图和distance的图，
shift-distance中的代码为生成figures/figure1/shift-distance中各物种的CG距离的代码。

## Fig3 cmd

### Corplot.R

该脚本是计算methylation 和 CpGdensity 的相关性，包括实验数据和随机数据。


### FeatureBin.R

该脚本是分bin 观察 methylation 和 CpGdensity 的关系，包括 gene，Line和LTR。

## Fig4 cmd

### Complement.R

该脚本是绘制 methylation 和 CpGdensity 符合情况。

## Promoter 和 Repeat cmd

### mm9Promoter.R

该脚本是为了绘制  promoter heatmap。

### Icmp.R

该脚本是统计 bonemarrow overlap or gap region 在 other tissues 是不是 overlap or gap， 同时绘制饼图。

### OverlapwithRepeat.R

该脚本是计算 bonemarrow overlap 在 Repeat 区域的分布情况，同时绘制饼图。 
