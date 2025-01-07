README

Documentation for data for the article `Belief Shocks and Implications of Expectations about Growth-at-Risk` by Maximilian Boeck and Michael Pfarrhofer. In this document, we describe the data sources involved.


Survey of Professional Forecasters by the Federal Reserve Bank of Philadelphia:
	- Link: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters
	- meanGrowth.xlsx (Link: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/mean-forecasts)
	- medianGrowth.xlsx (Link: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts)
	- SPFmicrodata.xlsx (Link: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/individual-forecasts)
	- Documentation is online available.

Real-Time Data Set by the Federal Reserve Bank of Philadelphia:
	- REAL GNP/GDP (ROUTPUT) routput_first_second_third.xlsx (Link: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/first-second-third)

Federal Reserve Economic Data (FRED) by the Federal Reserve Bank of St. Louis:
	- NFCI.xls: Chicago Fed National Financial Conditions Index (Link: https://fred.stlouisfed.org/series/NFCI)
	- OUTBS.xls: Business Sector: Real Value-Added Output for All Workers (Link: https://fred.stlouisfed.org/series/OUTBS)

facs_realtime:
	- extracted factors from the real-time data set by the Federal Reserve Bank of Philadelphia and series from the FRED database
	- all variables are listed in Table B1 in the Appendix
	- saved in facs_realtime.rda