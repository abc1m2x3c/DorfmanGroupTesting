# Dorfman testing for COVID19 data
This repository contains R codes for the article "Using Dorfman group testing during the COVID-19 outbreak in
the United States".
To use the program in this repository, please include the following command at the beginning of the program: 
```
source("https://raw.githubusercontent.com/abc1m2x3c/DorfmanGroupTesting/master/src.R")
```
In addition, this program requires 3 external R packages 'maxLik', 'rgl', and 'reshape'.

## A simple example 
Below is an example of R codes to manually input and analyze a Dorfman testing data set.

### Step 1. Clean the memory, load the required R packages, and source the code in this repository  
```
rm(list=ls(all=TRUE))
library(maxLik)
library(rgl)
library(reshape)
source("https://raw.githubusercontent.com/abc1m2x3c/DorfmanGroupTesting/master/src.R")
```
### Step 2. Input the data
```
my.data=interative.input()
```
This is an interactive procedure. Users can input as much data as needed. A sample display of this procedure is shown below.
```
##################################
######Please input your data######
##################################
Enter pool size: 5
Enter total number of pools: 30
Enter number of postive pools: 9
Enter number of positive individuals in each postive pool: 1,1,1,1,1,2,1,2,1
Is this the end of your data (Y/N)? N
##################################
######Please input your data######
##################################
Enter pool size: 4
Enter total number of pools: 1
Enter number of postive pools: 0
Is this the end of your data (Y/N)? Y
```

### Step 3. Data analysis
First, calculate the maximum likelihood estimates and the confidence intervals for each of the parameters p (prevalence), Se (sensitivity), and Sp (specificity).
```
DT.fit(my.data)
```
The sample output is shown below. The columns from left to right are: maximum likelihood estimation, standard error, lower bound of the 95% Wald confidence interval, upper bound of the 95% Wald confidence interval, lower bound of the 95% profile confidence interval, upper bound of the 95% profile confidence interval.
```
          est    std err  95% Wald L 95% Wald U 95% Profile L 95% Profile U
p  0.07172404 0.03293992 0.007162988  0.1362851    0.03483162     0.1238392
Se 0.99999864 0.10675137 0.790769793  1.2092275    0.76901045     1.0000000
Sp 0.99999726 0.04123922 0.919169876  1.0808246    0.90546264     1.0000000
```
Secondly, plot the Wald-type joint confidence region of (p, Se, Sp).
```
DT.plot(my.data,type='Wald')
```
The sample output is shown below. It is a 3D plot that can be rotated interactively using R command 'plot3D'.
![](https://github.com/abc1m2x3c/DorfmanGroupTesting/blob/master/Wald_example2.JPG)


Lastly, calculate the Wald confidence intervals of group testing charateristics for each pool size.
```
DT.chr(my.data)
```
The sample output is
```
$`pool size 5`
             L         U
E(T) 0.3448607 0.6766145
E(C) 0.9893016 1.0106967
PPV  0.7306390 1.2693427
NPV  0.9676666 1.0323330

$`pool size 4`
             L         U
E(T) 0.3680058 0.6469566
E(C) 0.9849983 1.0150003
PPV  0.7906560 1.2093298
NPV  0.9676666 1.0323330
```
