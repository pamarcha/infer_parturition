infer\_parturition
================
Pascal MARCHAND

## Code for paper: “A standardised biologging approach to infer parturition: an application in large herbivores across the hider-follower continuum”

Marchand P., Garel M., Morellet N., Benoit L., Chaval Y., Itty C., Petit
E., Cargnelutti B., Hewison A.J.M., Loison A. (202X) Journal, X(Y),
XXX-YYY.

#### Paper can be accessed at: <https://doi.org/XXX>

### Using this code:

The R script
‘Marchandetal\_infer\_parturition\_R\_script\_biologging\_data\_mouflon.R’
is given for readers that wish to replicate our analyses or applicate
our approach to other species or other behavioural metrics. The data
file contains biologging data and information on the reproductive status
of Mediterranean mouflon analysed in the manuscript (training and test
data sets; see text for details). .

#### IMPORTANT NOTE:

Due to changes in the default method for generating random sample from a
discrete uniform distribution (sample function) in R 3.6.0, a R version
\>=3.6.0 is absolutely needed to run the following script without error
messages and to get the same results as those presented in the
manuscript.ython 2.7 is required and nonstandard python packages
necessary are: numpy, scipy, keras, tables, and hdf5storage

For details:

<https://community.rstudio.com/t/getting-different-results-with-set-seed/31624/5>

<https://github.com/wch/r-source/blob/7f6cc784523dfa69087958633f7deb309d9c8718/doc/NEWS.Rd#L150-L161>:

<https://stackoverflow.com/questions/48626086/same-seed-different-os-different-random-numbers-in-r>

#### File Details:

  - Marchandetal\_infer\_parturition\_biologging\_data\_mouflon.RData
    This file contains several R objects:
    
        * info_train --> id and estimated date of parturition (based on change point analysis, see step 1) of the 29 females observed both before and after parturition = INFORMATION FOR TRAINING DATA SET
        
        * info_test --> id and observed reproductive status of the 28 females observed with a lamb at heel after parturition only (hence = parturient), repeatedly observed without lamb (and considered non parturient), or not observed after parturition (reproductive status unknown) = INFORMATION FOR TEST DATA SET
        
        * data_train = behavioural metrics of individuals (id) included in the training data set, derived from GPS locations and associated biologgers at each point in time (date) during the study period. The name of the columns correspond to the behavioural metrics as described in the Methods section of the manuscript
        
        * data_test = behavioural metrics of individuals (id) included in the test data set, derived from GPS locations and associated biologgers at each point in time (date) during the study period. The name of the columns correspond to the behavioural metrics as described in the Methods section of the manuscript

in both data\_train and data\_test objects, individual data are centered
= for each individual and each metric, we substracted the average value
computed over the study period
