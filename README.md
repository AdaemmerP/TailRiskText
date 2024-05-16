# TailRiskText
This repository contains the source code and data for the paper **Forecasting Macroeconomic Tail Risk in Real Time: Do Textual Data Add Value?** by [Adämmer, Prüser and Schüssler (2024)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4372186). 
Please cite the paper when using the code and data from this repository. 
We used newspaper articles from *The New York Times* and *The Washington Post*. 
The data were collected using the online tool [NexisUni](https://www.lexisnexis.com/en-int/products/nexis-uni). 
Due to licensing agreements, we do not upload the newspaper articles itself, but all computed time series derived from it. 
The vintage data we used are publicly available [here](https://research.stlouisfed.org/econ/mccracken/fred-databases/).

## Data
This folder contains all the time series data we used for quantile forecasting.
The 'FREDMD' folder contains the (vintage) macroeconomic data. The 'SentTopics' folder contains 
all the time series data based on the basis of textual data. The file 'sentiment_daily.RData' contains for each newspaper article the computed sentiment score. The data are separated by topic numbers. 
For example, the folder 'K80' contains all calculated data based on K=80 topics. The file 
'ctm_80_10_topicsonlymonthly.csv' corresponds to a CTM model with 80 topics and the top 10k words based on tf-idf.
The top words for each topic model are also included in the subfolders. 

## Code_Figures
This folder contains all R codes to replicate each Figure in the paper. The scripts are self-explanatory. 
At the beginning of each script you have to set the appropriate path, which is explained at the beginning. 
The folder also contains an "renv.lock", which contains information on the package versions used. 

## Code_Matlab
This folder contains all Matlab scripts and data to replicate the empirical results:

 - Use 'outofsample_main.m' to produce the predictions.
 - Use 'results_summary.m' to evaluate the predictions and save the results.
 - Use 'beta_table.m' to select the important variables for the predictions and save these.
 - The folder "data" contains the time series based on textual data.
 - The folder "Historical FRED-MD Vintages Final" contains the fred vintage data.
 - The folder "functions" contains the estimation functions.
 - The folder "Results" contains the calculated predictions and the actual values over time as Matlab-files.


## Results
This folder contains all empirical results based on the code in the folder "Code_Matlab".


