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
