#   MAT8480 Data Mining: Project Phase II
#   Lift Maximizers: Aaron Siegel, Sebastian Tilson, Dan Preston, and Zack Scalyer.
#   File: ProjectScoring.R
#   Created 11/25/2019
#
#   Purpose - To apply data preparation and predictive modeling to the scoring data.
#
#   input
#   %%%%%%%%
#   scoring_data.csv: data for scoring with out target.
#   RF_quickGMM.RData
#     - RF.fits: list of three (was a named list to compare best fits in complete tuning scrip)
#         [[1]] c(fscore.valid, fscore.test, fscore.oob): not used.
#         [[2]] RF.Final.VIM: not used.
#         [[3]] RF.Final: RF final fit
#   DM_GMM_fit_subj.RData
#     - fit_subj: Mclust Gaussian Mixture Modeling survey data inputs.
#   data_clean.RData
#     - data.clean: not used.
#     - cleanMeasure: list of character vectors identifying measurements corrections and
#         rejected inputs by column names to be used in scoring.
#   data_prep.RData
#     - data.prep: not used.
#     - inp.n: not used.
#     - indx: not used.
#     - split: not used.
#     - split.valid: not used.
#     - split.test: inot used.
#     - prepInput: list of keys to recode inputs to be used in scoring.
#   data_imp.RData
#     - data.imp: not used.
#     - inp.n: not used.
#     - indx: not used.
#     - split: not used.
#     - split.valid: not used.
#     - split.test: not used.
#     - imp.parms.num: imputation parameter median for num inputs
#     - imp.parms.cat: imputation parameter mode for cat inputs
#
#   output
#   %%%%%%%%
#   LiftMax_Score_RFfinal.csv: data.frame with EvalID and predict.class.
#
#   Notes
#   %%%%%%%%
#   1. All libraries used within this script are:
#         tidyverse, randomForest, caret, mlr, mclust, lubridate,
#         & tree.bins.
#
#   Updates
#   %%%%%%%%
#   1. Phase I implamented rule for setting zero values in ordenal inputs (survey data)
#        as NA; this rule is no longer inplament and comment out below.
#   2. 12/08/2019 & 12/14/2019 updated with clusters from GMM fit.
#
#%%%%%%%%%%%%#
#% Clean up %#
#%%%%%%%%%%%%#
rm(list = ls()) # clear the workspace
cat("\014")     # clear the console
graphics.off()  # clear all plots
#%%%%%%%%%%%%#
#% Start up %#
#%%%%%%%%%%%%#
library(tidyverse)
library(randomForest)
library(caret)
library(mlr)    # masked _by_ ‘.GlobalEnv’:impute,
                # masked from ‘package:caret’: train
library(mclust) # masked from ‘package:purrr’:map

#%%%%%%%%%%#
#% Set up %#
#%%%%%%%%%%#
homeDir <- "C:/Users/zack0/Documents/School/Villanova/FA_2019/MAT_8480/Project/LiftMax_DB/Final/Submittion"
setwd(homeDir) # Set the working directory
getwd()        # working directory
## helper functions ----
inx <- function (data, inp.n) { # data: current dataframe; inp.n: position for non-inputs
  # numeric input indicator
  num <- sapply(data, is.numeric)
  num[inp.n] <- FALSE

  # nominal input indicator
  cat <- sapply(data, is.factor)
  cat[inp.n] <- FALSE

  # missing value indicator
  na <- sapply(data,function(x) any(is.na(x)))
  na[inp.n] <- FALSE

  data.frame(num, cat, na)
}

z.impute <- function(x,y) { # impute x with y
  x[is.na(x)]<-y
  x
} # was renamed to aviod conflict with mlr::impute

funscore <- function(datScore, cleanMeasure, prepInput,
                     imp.parms.num, imp.parms.cat, na.inputs,
                     GMM.fit, RF.Final){
  raw<-datScore
  #% measurement corrections
  data <- datScore %>%
    # data cleaning
    mutate_at(.vars = cleanMeasure$measure.char, as.character) %>%         # measurement correction
    mutate_at(.vars = cleanMeasure$measure.cat[-length(cleanMeasure$measure.cat)],
              as.factor) %>%                                               # measurement correction
    mutate_at(.vars = cleanMeasure$measure.date,                           # measurement correction
              ~as.Date(.,format = "%m/%d/%Y")) %>%
    mutate_at(.vars = cleanMeasure$measure.time,                           # measurement correction
              ~lubridate::hms(., quiet = TRUE)) %>%
    # mutate_at(.vars = cleanMeasure$measure.ord[-c(10,11)], # data correction (ref Update: 1)
    #           ~na_if(.,0)) %>%                             #   where NA is assumed Zero
    mutate_at(.vars = vars(contains("Connections")),       # data correction
              abs) %>%                                     #   for {"Outbound_Connections", "Return_Connections"}
    mutate_at(.vars = vars(contains("Connect_Time")),      # data correction (step 1)
              ~recode(., `-1` = 0L)) %>%                   #   for "Outbound..." & "Return..." "_Connect_Time_Mins..." "_1" & "_2"
    mutate_at(.vars = vars(contains("Connect_Time")),      # data correction (step 2)
              abs) %>%                                     #   where step 1 (change -1 to zero) then step 2 (take absoulte value)
    #% cleaning factor levels
    mutate(Grp_Size_Cat = recode(Grp_Size_Cat, !!!prepInput$key.Grp_Size_Cat,
                                 .default = NA_character_), # all values not otherwise matched to NA
           Tour_Region = recode(Tour_Region, !!!prepInput$key.Tour_Region,
                                .default = "All Other"), # all values not otherwise matched to "All Other"
           # TourCode = recode(TourCode, !!!prepInput$key.TourCode,
           #                   .default = "All Other"), # all values not otherwise matched to "All Other"
           Book_Months = recode_factor(Book_Months, !!!prepInput$key.Book_Months,
                                       .default = NA_character_), # all values not otherwise matched to NA
           Age = recode_factor(Age, !!!prepInput$key.Age,
                               .default = NA_character_), # all values not otherwise matched to NA
           DB_Enter_Months = recode_factor(DB_Enter_Months, !!!prepInput$key.DB_Enter_Months,
                                           .default = NA_character_), # all values not otherwise matched to NA
           Email = recode(Email, !!!prepInput$key.Email,
                          .default = NA_character_), # all values not otherwise matched to NA
           Past_Trips = recode(Past_Trips, !!!prepInput$key.Past_Trips,
                               .default = "2 Or More"), # all values not otherwise matched to "2 Or More"
           TourPriceCat = recode_factor(TourPriceCat, !!!prepInput$key.TourPriceCat,
                                        .default = NA_character_), # all values not otherwise matched to NA
           TravelAgain = recode_factor(TravelAgain, !!!prepInput$key.TravelAgain,
                                       .default = NA_character_)) %>%  # all values not otherwise matched to NA
    mutate(Capacity_Cat = recode(Capacity, !!!prepInput$key.Capacity_Cat, # Define new input, Capacity_Cat
                                 .default = "All Other"),       #   all values not otherwise matched to "All Other"
           Capacity_Cat = as.factor(Capacity_Cat),              #   Capacity_Cat as factor
           Contact_Event = if_else(Eval_Contact_Days > 0,1, 0), # Define new input, Contact_Event
           Contact_Event = as.factor(Contact_Event),            #   Contact_Event as factor
           Total_Outbound_Connect_Time = Outbound_Connect_Time_Mins_1+Outbound_Connect_Time_Mins_2, # Define new input, Total_Outbound_Connect_Time
           Total_Return_Connect_Time = Return_Connect_Time_Mins_1+Return_Connect_Time_Mins_2, # Define new input, Total_Return_Connect_Time
           Outbound_Connect_cat = as.factor((!is.na(Outbound_Connect_Gateway1)) + (!is.na(Outbound_Connect_Gateway2))), # Define new input see Note: 4
           Return_Connect_cat = as.factor((!is.na(Return_Connect_Gateway1)) + (!is.na(Return_Connect_Gateway2))), # Define new input see Note: 4
           TourDate_WeekYear = strftime(TourDate, format="%V-%Y"), # Aggregate TourDate at the weekly level, returns character vector
           TourDate_WeekYear = as.factor(TourDate_WeekYear)) %>% #  TourDate_WeekYear as factor
    mutate_at(.vars = cleanMeasure$measure.time,                   # Define new input,
              list(AM = ~lubridate::am(.))) %>%       #   time mesurments as logical flag for AM or not
    mutate_at(.vars = sapply(cleanMeasure$measure.time, function(x) paste0(x, "_AM"), USE.NAMES = F),
              as.factor) %>%                          #   *flight*_Time_AM to factor
    mutate_at(.vars = sapply(cleanMeasure$measure.time, function(x) paste0(x, "_AM"), USE.NAMES = F),
              ~recode(.,`FALSE` = "0", `TRUE` = "1")) #   with binary labels

  #% recategorize inputs from lkup.list
  # Bins defined on training data may not define all possible levels
  data %>%
    filter(!State %in% unique(prepInput$lkup.list[[1]]$State)) %>%
    drop_na(State) %>% # do not want to impute na values here
    pull(State) %>%
    unique() %>%
    as.character() -> list.nonlevls_State
  if (!(length(list.nonlevls_State) == 0)){
    prepInput$lkup.list[[1]] <- prepInput$lkup.list[[1]] %>%
      add_row(State = list.nonlevls_State,
              Categories = rep("bin.2",length(list.nonlevls_State)))
  }

  data %>%
    filter(!SourceType %in% unique(prepInput$lkup.list[[2]]$SourceType)) %>%
    drop_na(SourceType) %>% # do not want to impute na values here
    pull(SourceType) %>%
    unique() %>%
    as.character() -> list.nonlevls_SourceType
  if (!(length(list.nonlevls_SourceType) == 0)){
    prepInput$lkup.list[[2]] <- prepInput$lkup.list[[2]] %>%
      add_row(SourceType = list.nonlevls_SourceType,
              Categories = rep("bin.3",length(list.nonlevls_SourceType)))
  }

  data %>%
    filter(!TourCode %in% unique(prepInput$lkup.list[[3]]$TourCode)) %>%
    drop_na(TourCode) %>% # do not want to impute na values here
    pull(TourCode) %>%
    unique() %>%
    as.character() -> list.nonlevls_TourCode
  if (!(length(list.nonlevls_TourCode) == 0)){
    prepInput$lkup.list[[3]] <- prepInput$lkup.list[[3]] %>%
      add_row(TourCode = list.nonlevls_TourCode,
              Categories = rep("bin.1",length(list.nonlevls_TourCode)))
  }

  data %>%
    filter(!TourDate_WeekYear %in% unique(prepInput$lkup.list[[4]]$TourDate_WeekYear)) %>%
    drop_na(TourDate_WeekYear) %>% # do not want to impute na values here
    pull(TourDate_WeekYear) %>%
    unique() %>%
    as.character() -> list.nonlevls_TDWY
  if (!(length(list.nonlevls_TDWY) == 0)){
    prepInput$lkup.list[[4]] <- prepInput$lkup.list[[4]] %>%
      add_row(TourDate_WeekYear = list.nonlevls_TDWY,
              Categories = rep("bin.1",length(list.nonlevls_TDWY)))
  }

  data %>%
    filter(!Outbound_Domestic_Gateway %in% unique(prepInput$lkup.list[[5]]$Outbound_Domestic_Gateway)) %>%
    drop_na(Outbound_Domestic_Gateway) %>% # do not want to impute na values here
    pull(Outbound_Domestic_Gateway) %>%
    unique() %>%
    as.character() -> list.nonlevls_ODG
  if (!(length(list.nonlevls_ODG) == 0)){
    prepInput$lkup.list[[5]] <- prepInput$lkup.list[[5]] %>%
      add_row(Outbound_Domestic_Gateway = list.nonlevls_ODG,
              Categories = rep("bin.1",length(list.nonlevls_ODG)))
  }

  data %>%
    filter(!Outbound_Intr_Gateway %in% unique(prepInput$lkup.list[[6]]$Outbound_Intr_Gateway)) %>%
    drop_na(Outbound_Intr_Gateway) %>% # do not want to impute na values here
    pull(Outbound_Intr_Gateway) %>%
    unique() %>%
    as.character() -> list.nonlevls_OIG
  if (!(length(list.nonlevls_OIG) == 0)){
    prepInput$lkup.list[[6]] <- prepInput$lkup.list[[6]] %>%
      add_row(Outbound_Intr_Gateway = list.nonlevls_OIG,
              Categories = rep("bin.1",length(list.nonlevls_OIG)))
  }

  data %>%
    filter(!Return_Domestic_Gateway %in% unique(prepInput$lkup.list[[7]]$Return_Domestic_Gateway)) %>%
    drop_na(Return_Domestic_Gateway) %>% # do not want to impute na values here
    pull(Return_Domestic_Gateway) %>%
    unique() %>%
    as.character() -> list.nonlevls_RDG
  if (!(length(list.nonlevls_RDG) == 0)){
    prepInput$lkup.list[[7]] <- prepInput$lkup.list[[7]] %>%
      add_row(Return_Domestic_Gateway = list.nonlevls_RDG,
              Categories = rep("bin.1",length(list.nonlevls_RDG)))
  }

  data %>%
    filter(!Return_Intr_Gateway %in% unique(prepInput$lkup.list[[8]]$Return_Intr_Gateway)) %>%
    drop_na(Return_Intr_Gateway) %>% # do not want to impute na values here
    pull(Return_Intr_Gateway) %>%
    unique() %>%
    as.character() -> list.nonlevls_RIG
  if (!(length(list.nonlevls_RIG) == 0)){
    prepInput$lkup.list[[8]] <- prepInput$lkup.list[[8]] %>%
      add_row(Return_Intr_Gateway = list.nonlevls_RIG,
              Categories = rep("bin.2",length(list.nonlevls_RIG)))
  }

  data %>%
    filter(!Outbound_Connect_Gateway1 %in% unique(prepInput$lkup.list[[9]]$Outbound_Connect_Gateway1)) %>%
    drop_na(Outbound_Connect_Gateway1) %>% # do not want to impute na values here
    pull(Outbound_Connect_Gateway1) %>%
    unique() %>%
    as.character() -> list.nonlevls_OCG1
  if (!(length(list.nonlevls_OCG1) == 0)){
    prepInput$lkup.list[[9]] <- prepInput$lkup.list[[9]] %>%
      add_row(Outbound_Connect_Gateway1 = list.nonlevls_OCG1,
              Categories = rep("bin.1",length(list.nonlevls_OCG1)))
  }

  data %>%
    filter(!Outbound_Connect_Gateway2 %in% unique(prepInput$lkup.list[[10]]$Outbound_Connect_Gateway2)) %>%
    drop_na(Outbound_Connect_Gateway2) %>% # do not want to impute na values here
    pull(Outbound_Connect_Gateway2) %>%
    unique() %>%
    as.character() -> list.nonlevls_OCG2
  if (!(length(list.nonlevls_OCG2) == 0)){
    prepInput$lkup.list[[10]] <- prepInput$lkup.list[[10]] %>%
      add_row(Outbound_Connect_Gateway2 = list.nonlevls_OCG2,
              Categories = rep("bin.1",length(list.nonlevls_OCG2)))
  }

  data %>%
    filter(!Return_Connect_Gateway1 %in% unique(prepInput$lkup.list[[11]]$Return_Connect_Gateway1)) %>%
    drop_na(Return_Connect_Gateway1) %>% # do not want to impute na values here
    pull(Return_Connect_Gateway1) %>%
    unique() %>%
    as.character() -> list.nonlevls_RCG1
  if (!(length(list.nonlevls_RCG1) == 0)){
    prepInput$lkup.list[[11]] <- prepInput$lkup.list[[11]] %>%
      add_row(Return_Connect_Gateway1 = list.nonlevls_RCG1,
              Categories = rep("bin.3",length(list.nonlevls_RCG1)))
  }

  data %>%
    filter(!Return_Connect_Gateway2 %in% unique(prepInput$lkup.list[[12]]$Return_Connect_Gateway2)) %>%
    drop_na(Return_Connect_Gateway2) %>% # do not want to impute na values here
    pull(Return_Connect_Gateway2) %>%
    unique() %>%
    as.character() -> list.nonlevls_RCG2
  if (!(length(list.nonlevls_RCG2) == 0)){
    prepInput$lkup.list[[12]] <- prepInput$lkup.list[[12]] %>%
      add_row(Return_Connect_Gateway2 = list.nonlevls_RCG2,
              Categories = rep("bin.3",length(list.nonlevls_RCG2)))
  }

  # recategorize inputs from lkup.list with all data as data.prep
  data.prep <- tree.bins::bin.oth(list = prepInput$lkup.list,
                                  data = data) %>%
    mutate(State = as.factor(State),
           SourceType = as.factor(SourceType),
           TourCode = as.factor(TourCode),
           TourDate_WeekYear = as.factor(TourDate_WeekYear),
           Outbound_Domestic_Gateway = as.factor(Outbound_Domestic_Gateway),
           Outbound_Intr_Gateway = as.factor(Outbound_Intr_Gateway),
           Return_Domestic_Gateway = as.factor(Return_Domestic_Gateway),
           Return_Intr_Gateway = as.factor(Return_Intr_Gateway),
           Outbound_Connect_Gateway1 = as.factor(Outbound_Connect_Gateway1),
           Outbound_Connect_Gateway2 = as.factor(Outbound_Connect_Gateway2),
           Return_Connect_Gateway1 = as.factor(Return_Connect_Gateway1),
           Return_Connect_Gateway2 = as.factor(Return_Connect_Gateway2))

  # Imputation
  inp.n <- which(colnames(data.prep) %in% c(cleanMeasure$rejected.inputs, "Book_12Mo"))
  indx <- inx(data.prep, inp.n)

  data.imp <- data.prep %>%
    mutate_if(.predicate = (names(.) %in% na.inputs[na.inputs %in% names(data.prep[indx$na])]),
              list(na = ~ as.factor(ifelse(is.na(.),1,0)))) %>%
    mutate_if(.predicate = (names(.) %in% names(data.prep[indx$na])[!(names(data.prep[indx$na]) %in% na.inputs)]),
              ~replace_na(.,names(table(.))[table(.)==max(table(.))]))

  data.imp[which(names(data.imp) %in% names(data.prep[indx$num]))]<-as.data.frame(mapply(z.impute,
                                                                                         x=data.imp[which(names(data.imp) %in% names(data.prep[indx$num]))],
                                                                                         y = imp.parms.num))
  data.imp[which(names(data.imp) %in% names(data.prep[indx$cat]))]<-as.data.frame(mapply(z.impute,
                                                                                         x=data.imp[which(names(data.imp) %in% names(data.prep[indx$cat]))],
                                                                                         y = imp.parms.cat))


  # Gaussian Mixture # names(data.imp_dummies);dim(data.imp_dummies) # 10054   262
  vars.gmm <- which(names(data.imp) %in% attr(GMM.fit$data, "dimnames")[2][[1]])
  data.imp_dummies <- createDummyFeatures(data.imp[vars.gmm]) # dummy vars


  classification <- predict(GMM.fit, data.imp_dummies)      # appy GMM.fit to imp scoring data
  data.imp.GMM <- data.imp %>%                              # pass data
    add_column(cluster = as.factor(classification$classification)) # add cluster input as factor.

  # score
  vars.subset <- which(colnames(data.imp.GMM) %in% attr(RF.Final[["terms"]], "term.labels"))
  predict.class.RF <- predict(RF.Final, newdata=data.imp.GMM[vars.subset], type="class")
  raw[c("predict.class")] <- data.frame(predict.class.RF)
  return(raw[c("EvalID","predict.class")])

  # reg.vars <- which(colnames(data.imp) %in% attr(reg.step$terms, "term.labels"))
  # predict.prob <- predict(reg.step, newdata=data.imp[reg.vars], type = "response")
  # predict.class <- as.factor(ifelse(predict.prob >= regThresh$threshold, 1,0))
}


# loading cleanMeasure from data_clean.RData
print(load("data_clean.RData"))
rm(data.clean)

# loading prepInput, & na.inputs from data_prep.RData
print(load("data_prep.RData"))
na.inputs <- names(data.prep[indx$na]) # names of na flag inputs
rm(inp.n,indx,split,split.valid,split.test,data.prep)

# # loading imp.parms.num, & imp.parms.cat from data_imp.RData
print(load("data_imp.RData"))
rm(inp.n,indx,split,split.valid,split.test,data.imp)

# loading Gaussian Mixture clustering model
print(load("DM_GMM_fit_subj.RData"))
GMM.fit <- fit_subj # pass model
rm(fit_subj)

# loading RF model
print(load("RF_quickGMM.RData"))
RF.Final <- RF.fits[[3]]

# loading scoring data
datScore <- read.csv("scoring_data.csv", header=TRUE, na.strings=c(".", "NA", "", "?"))

# predictions
dfpredict <- funscore(datScore, cleanMeasure, prepInput,
                       imp.parms.num, imp.parms.cat, na.inputs,
                       GMM.fit, RF.Final)


table(dfpredict$predict.class)
# > table(dfpredict$predict.class)
#
# 0    1
# 8217 1837

write.csv(dfpredict, file="LiftMax_Score_RFfinal.csv")
