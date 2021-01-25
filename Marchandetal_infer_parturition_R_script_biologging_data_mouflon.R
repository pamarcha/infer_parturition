##########################################################
##     A standardised biologging approach to infer      ##
##    parturition in wildlife, with application to      ##
## large herbivores across the hider-follower continuum ##
##########################################################

                                        # Marchand et al., submitted to  Methods in Ecology and Evolution

# This R script is given for readers that wish to replicate our analyses or applicate our
# approach to other species or other behavioural metrics. The data file contains biologging
# data and information on the reproductive status of Mediterranean mouflon analysed in the
                                        # manuscript (training and test data sets; see text for details).

####################
## IMPORTANT NOTE ##
####################

## Due to changes in the default method for generating random sample from a discrete
## uniform distribution (sample function) in R 3.6.0,  R version >=3.6.0 is absolutely
## needed to run the following script without error messages and to get the same results
## as those presented in the manuscript.
## For details:
## https://community.rstudio.com/t/getting-different-results-with-set-seed/31624/5
## https://github.com/wch/r-source/blob/7f6cc784523dfa69087958633f7deb309d9c8718/doc/NEWS.Rd#L150-L161:
## https://stackoverflow.com/questions/48626086/same-seed-different-os-different-random-numbers-in-r


###############
## LIBRARIES ##
###############

library(rfPermute)

###############
## LOAD DATA ##
###############

load("Marchandetal_infer_parturition_biologging_data_mouflon.RData")

                                        # This file contains :

                                        # info_train --> id and estimated date
                                        # of parturition (based on change point analysis,
                                        # see step 1) of the 29 females observed both
                                        # before and after parturition = INFORMATION FOR TRAINING DATA SET

                                        # info_test --> id and observed reproductive
                                        # status of the 28 females observed with a lamb at
                                        # heel after parturition only (hence =
                                        # parturient), repeatedly observed without lamb
                                        # (and considered non parturient), or not observed
                                        # after parturition (reproductive status unknown)
                                        # = INFORMATION FOR TEST DATA SET

                                        # data_train = behavioural metrics of individuals
                                        # (id) included in the training data set, derived
                                        # from GPS locations and associated biologgers at
                                        # each point in time (date) during the study
                                        # period. The name of the columns correspond to
                                        # the behavioural metrics as described in the
                                        # Methods section of the manuscript

                                        # data_test = behavioural metrics of individuals
                                        # (id) included in the test data set, derived
                                        # from GPS locations and associated biologgers at
                                        # each point in time (date) during the study
                                        # period. The name of the columns correspond to
                                        # the behavioural metrics as described in the
                                        # Methods section of the manuscript

                                        # in both data_train and data_test objects,
                                        # individual data are scaled = for each
                                        # individual and each metric, we substracted
                                        # the average value computed over the study period


behaviouralmetrics <- c("dist","RT100","HR","slope","rock","SS","FB","AI")
                                        # the name of the columns/behavioural metrics used to identify parturition in mouflon

###########################################################
##      PRELIMINARY STEP NEEDED FOR DECISION RULES :     ##
## DEFINE ERROR RATES OF THE BEST RANDOM FOREST APPROACH ##
###########################################################

## COMPUTE TIME DIFFERENCE BETWEEN EACH DATA AND THE ESTIMATED PARTURITION DATE (TRAINING DATA SET)
for(id in unique(data_train$id)){
    data_train$difftime_parturition[data_train$id %in% id] <- difftime(data_train$date[data_train$id %in% id],as.POSIXlt(info_train$estimated_date_of_parturition[info_train$id %in% id]),units="hours")
}

## CREATE CORRESPONDING COLUMNS IN TRAINING DATA SET
for(i in c(seq(24, 240, 12))){
    data_train[,paste("parturition_", i, "h", sep="")] <- ifelse(abs(data_train$difftime_parturition) <=i, paste("parturition +/-", i, sep=""), "rest of the study period")
    data_train[,paste("parturition_", i, "h", sep="")] <- as.factor(data_train[,paste("parturition_", i, "h", sep="")])
}


## FOR TIME WINDOWS BETWEEN 24 AND 240H, COMPUTE ERROR RATES FROM RANDOM FORESTS TRAINED WITH MOVEMENT+HABITAT_USE+ACTIVITY METRICS
resboot <- list()

for(i in c(seq(24, 240, 12))){
    time_window <- paste("parturition_", i, "h", sep="")
    behaviouralmetricsform <- paste(behaviouralmetrics, collapse=" + ")

    formul <- as.formula(paste(time_window, behaviouralmetricsform, sep="~"))

    sub_train <- data_train[apply(data_train[,behaviouralmetrics], 1, function(x){all(!is.na(x))}),]

    tw_train<-sub_train[!sub_train[,time_window] %in% "rest of the study period",]
    rsp_train<-sub_train[sub_train[,time_window] %in% "rest of the study period",] # rsp = rest of the study period

    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind="Rejection"); RNGkind(sample.kind="Rejection")
    rsp_train <-  rsp_train[sample(1:nrow(rsp_train), nrow(tw_train), replace=F),]
    train <- rbind.data.frame(tw_train, rsp_train)

    rf<-randomForest(formul, data=train)

    res <- as.data.frame(rfPermute::confusionMatrix(rf, threshold=0.5))

    global=(100-res$pct.correct[3])/100
    sup_global=(100-res$LCI_0.95[3])/100
    inf_global=(100-res$UCI_0.95[3])/100
    fp=(100-res$pct.correct[2])/100
    sup_fp=(100-res$LCI_0.95[2])/100
    inf_fp=(100-res$UCI_0.95[2])/100
    fn=(100-res$pct.correct[1])/100
    sup_fn=(100-res$LCI_0.95[1])/100
    inf_fn=(100-res$UCI_0.95[1])/100

    resboot[[as.character(i)]] <- data.frame(per=i, global=global,inf_global=inf_global,sup_global=sup_global,fp=fp,inf_fp=inf_fp,sup_fp=sup_fp,fn=fn,inf_fn=inf_fn,sup_fn=sup_fn)

    ##print(paste("get error rates from RF and different time windows, step ", which(seq(24, 240, 12) %in% i), "/", length(seq(24, 240, 12)), " --> ", i, "h", sep=""))
}

simulrf <- do.call(rbind, resboot)
names(simulrf) <- c("tw",paste(c("avg", "inf", "sup"), rep(c("error_rate","fp","fn"), each=3), sep="_"))

##################################################################################################################
##                         PRELIMINARY STEP NEEDED FOR DECISION RULES :                                         ##
## DEFINE PROBABILITIES THAT EACH DATA WAS RECORDED IN THE FOCAL TIME WINDOW BY THE BEST RANDOM FOREST APPROACH ##
##################################################################################################################

sequence <- seq(48, 72, 12)

for(i in sequence){

    time_window <- paste("parturition_", i, "h", sep="")

    ## the formulae including all the behavioural metrics (movement+habitat_use+activity)
    behaviouralmetricsform <- paste(behaviouralmetrics, collapse=" + ")
    formul <- as.formula(paste(time_window, behaviouralmetricsform, sep="~"))

    sub_train <- data_train[apply(data_train[,behaviouralmetrics], 1, function(x){all(!is.na(x))}),]

    ## obtain balanced samples between data recorded during/outside the focal time window by resampling
    tw_train<-sub_train[!sub_train[,time_window] %in% "rest of the study period",] # tw = time window
    rsp_train<-sub_train[sub_train[,time_window] %in% "rest of the study period",] # rsp = rest of the study period

    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind="Rejection"); RNGkind(sample.kind="Rejection")
    rsp_train <-  rsp_train[sample(1:nrow(rsp_train), nrow(tw_train), replace=F),]
    train <- rbind.data.frame(tw_train, rsp_train)

    ## run random forest
    rf<-randomForest(formul, data=train)

    ## get probabilities predicted by RF for the training data set
    data_train[,paste("pred_prob", i, sep="")] <- predict(rf, newdata=data_train, type="prob")[,1]

    ##print(paste("training data set - predict probabilities to be in the focal time window, step ", which(sequence %in% i), "/", length(sequence), " --> ", i, "h", sep=""))

}

############################################################
## PREDICT PARTURITION PROBABILITIES IN THE TEST DATA SET ##
############################################################

for(i in sequence){

    time_window <- paste("parturition_", i, "h", sep="")

    ## the formulae including all the behavioural metrics (movement+habitat_use+activity)
    behaviouralmetricsform <- paste(behaviouralmetrics, collapse=" + ")
    formul <- as.formula(paste(time_window, behaviouralmetricsform, sep="~"))

    sub_train <- data_train[apply(data_train[,behaviouralmetrics], 1, function(x){all(!is.na(x))}),]

    tw_train<-sub_train[!sub_train[,time_window] %in% "rest of the study period",]
    rsp_train<-sub_train[sub_train[,time_window] %in% "rest of the study period",] # rsp = rest of the study period


    ## obtain balanced samples between data recorded during/outside the focal time window by resampling
    set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind="Rejection"); RNGkind(sample.kind="Rejection")
    rsp_train <-  rsp_train[sample(1:nrow(rsp_train), nrow(tw_train), replace=F),]
    train <- rbind.data.frame(tw_train, rsp_train)

    ## run random forest
    rf<-randomForest(formul, data=train)

    ## get probabilities predicted by RF for the test data set
    data_test[,paste("pred_prob", i, sep="")] <- predict(rf, newdata=data_test, type="prob")[,1]

    ##print(paste("test data set - predict probabilities to be in the focal time window, step ", which(sequence %in% i), "/", length(sequence), " --> ", i, "h", sep=""))

}

##########################
## APPLY DECISION RULES ##
##########################

library(zoo)

respred <- list()
for(id in unique(data_test$id)){
    sub<-data_test[data_test$id %in% id,]

    for(i in sequence){ ## the decisions rules are applied several times = for several time windows around parturition


        ## point in time where probabilities predicted by RF are above a defined threshold
        ## are considered as parturition

        ##the threshold is the 1% quantile of the
        ## probabilities predicted by RF on data from the corresponding time window in the
        ## training data set
        threshold <- quantile(data_train[abs(data_train$difftime_parturition) < i, paste("pred_prob",i,sep="")], 0.01, na.rm=T)


        ## which proportion of data are predicted as parturition in a sliding window corresponding to the time window
                slw <- i/2 ## length of the sliding window (/2 because we have 1 location/biologging data every 2 hours in mouflon)

                sub[,paste("mb",i, sep="")] <- as.vector(rollapply(zoo(sub[,paste("pred_prob",i,sep="")]), slw, function(x){mean(sum(as.numeric(x > threshold), na.rm=T), na.rm=T)}, fill=NA))/as.vector(rollapply(zoo(sub[,paste("pred_prob",i,sep="")]), slw, function(x){length(x[!is.na(x)])}, fill=NA))

        ## if there is a period when such proportion is higher than 1-error rate --> parturition identified !
        toto<-rle(as.vector(as.numeric(sub[,paste("mb",i,sep="")] > simulrf$sup_fn[simulrf$tw %in% i])))

        ## if parturition periods identified, keep only the most probable
        if(length(which(toto$values==1)) > 0){
            toto<-data.frame(end=cumsum(toto$lengths)[toto$values %in% "1"], duration=toto$lengths[toto$values %in% "1"])
            toto$parturition_start<-sub$date[toto$end - toto$duration + 1]
            toto$parturition_end<-sub$date[toto$end]
            toto$duration<-as.vector(difftime(toto$parturition_end, toto$parturition_start, unit="hours"))
            toto$id <- id
            toto$tw <- i
            toto<-toto[,c("id","tw","parturition_start","parturition_end","duration")]
            for(line in 1:nrow(toto)){
                toto$prop[line] <- mean(sub[sub$date >=toto$parturition_start[line] & sub$date <=toto$parturition_end[line],paste("mb",i,sep="")])
                toto$propmb[line] <- mean(sub[sub$date >=toto$parturition_start[line] & sub$date <=toto$parturition_end[line],paste("pred_prob",i,sep="")])
                toto$max[line] <- max(sub[sub$date >=toto$parturition_start[line] & sub$date <=toto$parturition_end[line],paste("pred_prob",i,sep="")])
                toto$mb[line] <- as.character(weighted.mean(sub$date[sub$date >=toto$parturition_start[line] & sub$date <=toto$parturition_end[line]], w=sub[sub$date >=toto$parturition_start[line] & sub$date <=toto$parturition_end[line],paste("pred_prob",i,sep="")],  na.rm=T))
            }
            toto <- toto[toto$max %in% max(toto$max),][1,]
            tmp <- toto
            respred[[paste(id, i, sep="-")]] <- tmp
        }
    }
}
detach("package:zoo")

respred <- do.call(rbind, respred)
respred <- droplevels(respred)
respred$mb <- strptime(respred$mb, "%Y-%m-%d %H:%M:%S")
respred <- respred[!duplicated(respred$id),] ## keep the shortest parturition period if several identified


## Compare results with reproductive status determined from repeated observations
res_test <- merge(info_test, respred, by="id", all.x=T)
res_test$predicted_parturition <- ifelse(!is.na(res_test$mb), "parturient", "non parturient")
table(predicted_RF=res_test$predicted_parturition, observed =res_test$observed_reproductive_status)

##                observed
## predicted_RF     non parturient parturient unknown
##   non parturient              2          2       0
##   parturient                  1         22       1

                                        #see Table 1 in the core text of the manuscript
