# Dorfman testing for COVID19 data
This repository contains R codes for the article "Using Dorfman group testing during the COVID-19 outbreak in
the United States".
To use the program in this repo, please include the following command at the beginning of the program: 
```
source("https://raw.githubusercontent.com/abc1m2x3c/DorfmanGroupTesting/master/src.R")
```
In addition, this program requires 3 external R packages 'maxLik', 'rgl', and 'reshape'.

## A simple example 
Below is an example of R codes to manually input and analyze a Dorfman testing Data set.

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
A sample display of the input procedure is
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
Calculate the maximum likelihood estimates and the confidence intervals for each of the parameters p (prevalence), Se (sensitivity), and Sp (specificity).
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
To plot the Wald-type joint confidence region of (p, Se, Sp),
```
DT.plot(my.data,type='Wald')
```
The sample output is
![Optional Text](../master/Example.png)
