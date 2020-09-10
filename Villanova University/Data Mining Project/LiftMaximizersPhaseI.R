#   MAT8480 Data Mining: Project Phase I
#   Lift Maximizers: Aaron Siegel, Sebastian Tilson, Dan Preston, and Zack Scalyer.
#   Created 10/31/2019
#
#   Purpose - To understand the nature of the data and relationship
#                among the variables.
#   input
#   %%%%%%%%
#   modeling_data.csv
#
#   output
#   %%%%%%%%
#   data_clean.RData
#     - data.clean: Rdata set with varable measurement corrections and data corrections
#     - cleanMeasure: list of character vectors identifying measurements corrections and
#         rejected inputs by column names to be used in scoring.
#   data_prep.RData
#     - data.prep: Rdata set with clean factor levels, new inputs, and recategorized inputs.
#     - inp.n: column index for non-input
#     - indx: data frame for input column index
#           indx$num: numeric inputs
#           indx$cat: categorical inputs
#           indx$na: inputs with at least one missing
#     - split: index for 50% training data
#     - split.valid: index for 25% vaidation data
#     - split.test: index for 25% testing data
#     - prepInput: list of keys to recode inputs to be used in scoring.
#   data_imp.RData
#     - data.imp: Rdata set with imputed data and na flag inputs.
#     - inp.n: column index for non-input
#     - indx: data frame for input column index
#           indx$num: numeric inputs
#           indx$cat: categorical inputs
#           indx$na: inputs with at least one missing (all FALSE)
#     - split: index for 50% training data
#     - split.valid: index for 25% vaidation data
#     - split.test: index for 25% testing data
#     - imp.parms.num: imputation parameter median for num inputs
#     - imp.parms.cat: imputation parameter mode for cat inputs
#   data_down.RData
#     - data.down: Rdata set with downsampled/balanced training data.
#     - inp.n: column index for non-input
#     - indx: data frame for input column index
#           indx$num: numeric inputs
#           indx$cat: categorical inputs
#           indx$na: inputs with at least one missing
#     - split.down: index for balanced training data
#     - split.down.valid: index for 25% vaidation data (preserved/unbalanced)
#     - split.down.test: index for 25% testing data (preserved/unbalanced)
#   data_down_imp.RData
#     - data.down.in: Rdata set of imputed data.down.
#     - inp.n: column index for non-input
#     - indx: data frame for input column index
#           indx$num: numeric inputs
#           indx$cat: categorical inputs
#           indx$na: inputs with at least one missing
#     - split.down: index for balanced training data
#     - split.down.valid: index for 25% vaidation data (preserved/unbalanced)
#     - split.down.test: index for 25% testing data (preserved/unbalanced)
#     - down.imp.parms.num: imputation parameter median for num inputs
#     - down.imp.parms.cat: imputation parameter mode for cat inputs
#
#   Notes
#   %%%%%%%%
#   1. Currently using the median and mode for imputation of interval and factor inputs,
#        calculated on 50% traning data "split" for data.imp in data_imp.Rdata,
#        and downsampled/balanced traning data "split.down" for data.down.imp in
#        data_down_imp.Rdata.
#   2. Individual Rdata files have been saved so model bilding can be done independent
#        of this scrip file while mataining data integerity.
#   3. As of 11/28/2019 Gaussian Mixture Model clusters where removed from data files
#        and will be included when methods are provied for scoring data.
#   4. All libraries used within this script are:
#        tidyverse, lubridate, caTools, rpart,  tree.bins, & caret.
#
#   Updates
#   %%%%%%%%
#   1. 11-21-2019 data.imp was updated in DMPhaseII_RandomForest.R
#        cluster, a factor input with 11 levels was added. A Gaussian
#        Mixture Model was used to crestet the clusters in data.imp in
#        another scrip made by Sebastian Tilson. File data_imp.RData
#        was updated with data.imp updated data.frame and indx flag.
#   2. 11-22-2019 data.down.imp was updated in DMPhaseII_RandomForest.R
#        cluster, a factor input with 11 levels was added. A Gaussian
#        Mixture Model was used to crestet the clusters in data.imp in
#        another scrip made by Sebastian Tilson. File data_down_imp.RData
#        was updated with data.down.imp updated data.frame and indx flag
#   3. 11-23-2019 Outbound_Connect_Gateway1 & 2 and Return_Connect_Gateway1
#        & 2 where recatergised with tree.bin in the same manor as the other
#        airport code inputs in PhaseI report. They are now no longer
#        rejected inputs.
#   4. 11-23-2019 Due to errors in Outbound_Connections and Return_Connections
#        two new categoric inputs where defined: Outbound_Connect_cat and
#        Return_Connect_cat with levels 0 for no connection flight code,
#        1 for one connection flight code, and 2 for two flight codes in data.
#        This theird level "2" represents two or more connections.
#   5. 11-23-2019 Updated all RData files with recatergised inputs and new
#        inputs in this scrip. Updated data_prep.RData, data_imp.RData, &
#        data_down_imp.RData with Gaussian clusters in DMPhaseII_RandomForest.R
#   6. 11/28/2019 Edits made to address issues found when scoring;
#        -Data Cleaning process with data_clean.RData output had accounted for
#          all measurement corrections and some data corrections. Now output
#          contains all data corrections and cleanMeasure was added for
#          scoring. Promo_Disc was added to the rejected input list due to
#          near zero variance.
#        -Data Prep process now recategorizes TourCode and TourDate_WeekYear
#          with tree.bins. TravelAgain is recoded with defalt value to address
#          conflits in scoring data. prepInput was added to data_prep.RData
#          output for scoring.
#   7. 12/01/2019 Ordenal inputs for survey data that had zero observed where
#        perviously and no longer taken to be missing.
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
#%%%%%%%%%%#
#% Set up %#
#%%%%%%%%%%#
homeDir <- "C:/Users/zack0/Documents/School/Villanova/FA_2019/MAT_8480/Project/LiftMax_DB/Final/Submittion"
setwd(homeDir) # Set the working directory
getwd()        # working directory set by R.Project: LiftMax_DB (edit for submition)

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

impute <- function(x,y) { # impute x with y
  x[is.na(x)]<-y
  x
}

FindMode <- function(x) { # find the first mode of x
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
} # Not used, will return <NA> when NA is the mode

z.mode <- function(x){ # find all modes of x
  return(names(table(x))[table(x)==max(table(x))])
}

## Data Cleaning ----
#% Inport data %#
data <- read.csv("modeling_data.csv", header=TRUE, na.strings=c(".", "NA", "", "?"))

#% measurement and data corrections ----
# character data
measure.char <- data %>%
  select(grep("ID$|Trip_no", names(data))) %>%
  names()
# date data
measure.date <- c("TourDate")
# time data
measure.time <- c("Domestic_Depart_Time", "Intr_Arrival_Time", "Intr_Depart_Time",
                  "Domestic_Arrival_Time")
# categorical inputs
measure.cat <- c("Recommend_GAT", "TravelAgain", "Groups_Interest", "Reference",
                 "Extension", "Insurance", "FltGty", "Complaint_Event", "Voucher_Event",
                 "Book_12Mo")
# ordinal inputs as number
measure.ord <- c("Overall_Impression", "Pre_Departure", "Flight_Itin", "TD_Overall",
                 "Hotels_Avg", "Meals_Avg", "GUSS_Avg", "Optionals_Avg", "Bus_Avg",
                 "Optionals","Eval_Contact_Days")
# ordinal inputs as factor
measure.ordfac <- c("Book_Months", "Age", "DB_Enter_Months", "Past_Trips", "TourPriceCat")

# inputs to drop
rejected.inputs <- c("EvalID", "Cus_ID", "ProdTour_ID", "SalesTourID", "Trip_no", "HH_ID",
                     "TourDate", "Domestic_Depart_Time", "Domestic_Arrival_Time",
                     "Intr_Arrival_Time", "Intr_Depart_Time", "Promo_Disc")

data <- data %>%
  mutate_at(.vars = measure.char, as.character) %>%         # measurement correction
  mutate_at(.vars = measure.cat, as.factor) %>%             # measurement correction
  mutate_at(.vars = measure.date,                           # measurement correction
            ~as.Date(.,format = "%m/%d/%Y")) %>%
  mutate_at(.vars = measure.time,                           # measurement correction
            ~lubridate::hms(., quiet = TRUE)) %>%
  # mutate_at(.vars = measure.ord[c(-10,-11)],                # data correction (no longer implemented)
  #           ~na_if(.,0)) %>%                                #   where Zero is NA (see Update 7)
  mutate_at(.vars = vars(contains("Connections")),          # data correction
            abs) %>%                                        #   for "Outbound_Connections" "Return_Connections"
  mutate_at(.vars = vars(contains("Connect_Time")),         # data correction (step 1)
            ~recode(., `-1` = 0L)) %>%                      #   for "Outbound..." & "Return..." "_Connect_Time_Mins..." "_1" & "_2"
  mutate_at(.vars = vars(contains("Connect_Time")),         # data correction (step 2)
            abs)                                            #   where step 1 (change -1 to zero) then step 2 (take absolute value)

#%% Pass & Save Data: clean ----
data.clean <- data # pass data
cleanMeasure <- list(measure.char=measure.char, measure.date=measure.date,
                     measure.time=measure.time, measure.cat=measure.cat,
                     measure.ord=measure.ord, measure.ordfac=measure.ordfac,
                     rejected.inputs=rejected.inputs) # for scoring
save(cleanMeasure, data.clean,  file="data_clean.RData")

## Data Prep ----

#% cleaning factor levels ----
## adjust Grp_Size_Cat label '43 to 45' to '44 to 45'
key.Grp_Size_Cat <- c(`25 to 29` = "25 to 29",
                      `30 to 34` = "30 to 34",
                      `35 to 39` = "35 to 39",
                      `40 to 43` = "40 to 43",
                      `43 to 45` = "44 to 45", # adjust label '43 to 45' to '44 to 45'
                      `Under 25` = "Under 25") # all values not otherwise matched to be NA
## Give EC its own bucket (12.3%), EU its own bucket (6.2)%, then assign all other unlabeled tour regions to 'All Other' (4.5%)
key.Tour_Region <- c(IT = "IT",
                     BI = "BI",
                     CNE = "CNE",
                     FS = "FS",
                     MD = "MD",
                     AM = "AM",
                     AF = "AF",
                     AU = "AU",
                     EC = "EC",
                     EU = "EU") # all values not otherwise matched to be "All Other"
# ## Assign all unlabeled TourCode to 'All Other' category (81.3%)
key.TourCode <- c(VFR = "VFR",
                  MIT = "MIT",
                  LPR = "LPR") # all values not otherwise matched to be "All Other"
key.Book_Months <- c(`Under 3months` = "Under4",
                     `4-6 months` = "4To6",
                     `7-9 months` = "7To9",
                     `10-12 months` = "10To12",
                     `Over 12months` = "Over12") # all values not otherwise matched to be NA
key.Age <- c(`Under 30` = "Under30",
             `30-39` = "30-39",
             `40-44` = "40-44",
             `40-44` = "40-44",
             `45-49` = "45-49",
             `50-54` = "50-54",
             `55-59` = "55-59",
             `60-69` = "60-69",
             `70-79` = "70-79",
             `Over 80` = "Over80") # all values not otherwise matched (No Age) to be NA
key.DB_Enter_Months <- c(`Under 3 months` = "Under4",
                         `4-6 months` = "4-6",
                         `7-12 months` = "7-12",
                         `Over 12 months` = "Over12") # all values not otherwise matched to be NA
key.Email <- c(Available = "Available",
               Unavailable = "Unavailable",
               Bounced = "Unavailable") # all values not otherwise matched to be NA
key.Past_Trips <- c(`0 Trips` = "0 Trips",
                    `1 Trip` = "1 Trip",
                    `2 Trip` = "2 Or More") # all values not otherwise matched to be "2 Or More"
key.TourPriceCat <- c(`Under 2000` = "Under 2000",
                      `2000 - 2500` = "2000 - 2500",
                      `2501 - 3000` = "2501 - 3000",
                      `3001 - 3500` = "3001 - 3500",
                      `3501 - 4000` = "3501 - 4000",
                      `4001 - 4500` = "4001 - 4500",
                      `4501 - 5000` = "4501 - 5000",
                      `More than 50` = "More than 5000") # all values not otherwise matched to be NA
key.TravelAgain <- c(`0` = "0",
                     `1` = "1",
                     `2` = "2") # all values not otherwise matched to be NA

data <- data %>%
  mutate(Grp_Size_Cat = recode(Grp_Size_Cat, !!!key.Grp_Size_Cat,
                               .default = NA_character_), # all values not otherwise matched to NA
         Tour_Region = recode(Tour_Region, !!!key.Tour_Region,
                              .default = "All Other"), # all values not otherwise matched to "All Other"
         # TourCode = recode(TourCode, !!!key.TourCode,
         #                   .default = "All Other"), # all values not otherwise matched to "All Other"
         Book_Months = recode_factor(Book_Months, !!!key.Book_Months,
                                     .default = NA_character_), # all values not otherwise matched to NA
         Age = recode_factor(Age, !!!key.Age,
                             .default = NA_character_), # all values not otherwise matched to NA
         DB_Enter_Months = recode_factor(DB_Enter_Months, !!!key.DB_Enter_Months,
                                         .default = NA_character_), # all values not otherwise matched to NA
         Email = recode(Email, !!!key.Email,
                        .default = NA_character_), # all values not otherwise matched to NA
         Past_Trips = recode(Past_Trips, !!!key.Past_Trips,
                             .default = "2 Or More"), # all values not otherwise matched to "2 Or More"
         TourPriceCat = recode_factor(TourPriceCat, !!!key.TourPriceCat,
                                      .default = NA_character_), # all values not otherwise matched to NA
         TravelAgain = recode(TravelAgain, !!!key.TravelAgain,
                                      .default = NA_character_) # all values not otherwise matched to NA
         )

#% mutate numerical input ----
## Bucket capacity into 45 (61%), 30 (11.5%), 35 (7.7%), and all other
key.Capacity_Cat <- c(`30` = "30",
                      `35` = "35",
                      `45` = "45") # all values not otherwise matched to be "All Other"

data <- data %>%
  mutate(Capacity_Cat = recode(Capacity, !!!key.Capacity_Cat, # Define new input, Capacity_Cat
                               .default = "All Other"),       #   all values not otherwise matched to "All Other"
         Capacity_Cat = as.factor(Capacity_Cat),              #   Capacity_Cat as factor
         Contact_Event = if_else(Eval_Contact_Days > 0,1, 0), # Define new input, Contact_Event
         Contact_Event = as.factor(Contact_Event),            #   Contact_Event as factor
         Total_Outbound_Connect_Time = Outbound_Connect_Time_Mins_1+Outbound_Connect_Time_Mins_2, # Define new input, Total_Outbound_Connect_Time
         Total_Return_Connect_Time = Return_Connect_Time_Mins_1+Return_Connect_Time_Mins_2, # Define new input, Total_Return_Connect_Time
         Outbound_Connect_cat = as.factor((!is.na(Outbound_Connect_Gateway1)) + (!is.na(Outbound_Connect_Gateway2))), # Define new input see Update: 4
         Return_Connect_cat = as.factor((!is.na(Return_Connect_Gateway1)) + (!is.na(Return_Connect_Gateway2))))       # Define new input see Update: 4

#% mutate date and time inputs ----
data <- data %>%
  mutate(TourDate_WeekYear = strftime(TourDate, format="%V-%Y"), # Aggregate TourDate at the weekly level, returns character vector
         TourDate_WeekYear = as.factor(TourDate_WeekYear)) %>%   #  TourDate_WeekYear as factor
  mutate_at(.vars = measure.time,                   # Define new input,
            list(AM = ~lubridate::am(.))) %>%       #   time measurements as logical flag for AM or not
  mutate_at(.vars = sapply(measure.time, function(x) paste0(x, "_AM"), USE.NAMES = F),
            as.factor) %>%                          #   *flight*_Time_AM to factor
  mutate_at(.vars = sapply(measure.time, function(x) paste0(x, "_AM"), USE.NAMES = F),
            ~recode(.,`FALSE` = "0", `TRUE` = "1")) #   with binary labels

## Derived Inputs
new.inputs <- names(data[,94:length(data)])

#% Split ----
set.seed(77012) # 50-25-25
split <- caTools::sample.split(data$Book_12Mo, SplitRatio = 0.5)
split.valid <- !split
split.test <- !split
split2 <- caTools::sample.split(data$Book_12Mo[!split], SplitRatio = 0.5)
split.valid[split.valid==TRUE] = split2
split.test[split.test==TRUE] = !split2

# checking 50-25-25 split
# prop.table(table(data$Book_12Mo)) # target proportion in data
# prop.table(table(split,data$Book_12Mo),1) # target proportion in training data
# prop.table(table(split.valid,data$Book_12Mo),1) # target proportion in validation data
# prop.table(table(split.test,data$Book_12Mo),1) # target proportion in testing data
# prop.table(table(split.test,split.valid)) # 50-25-25 split training-validation-testing

#% Set Index ----
# column index for all non-inputs
inp.n <- c(which(colnames(data.clean) %in% append(rejected.inputs,"Book_12Mo")))
names(data.clean[inp.n])
# generate num/factor/missing input indices
indx <- inx(data.clean, inp.n)
# names(data.clean[indx$num]) # numeric inputs
# names(data.clean[indx$cat]) # catagric inputs
# names(data.clean[indx$na]) # missing inputs

#% recategorize inputs ----
# exploratory tree: leaves are used to recategorize the current data within tree.bins()
explor.tree_State = rpart::rpart(formula = Book_12Mo ~ State,
                                 data = data[split,], # bins must be defined on training data
                                 method = "anova",
                                 control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01
summary(explor.tree_State) # from this I choose cp = 0.0019

explor.tree_SourceType = rpart::rpart(formula = Book_12Mo ~ SourceType,
                                      data = data[split,],
                                      method = "anova",
                                      control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01
summary(explor.tree_SourceType) # from this I choose cp = 0.0002

explor.tree_TourCode = rpart::rpart(formula = Book_12Mo ~ TourCode,
                                      data = data[split,],
                                      method = "anova",
                                      control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01
summary(explor.tree_TourCode) # from this I choose cp = 0.0015

explor.tree_TDWY = rpart::rpart(formula = Book_12Mo ~ TourDate_WeekYear,
                                    data = data[split,],
                                    method = "anova",
                                    control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01
summary(explor.tree_TDWY) # from this I choose cp = 0.003

explor.tree_ODG = rpart::rpart(formula = Book_12Mo ~ Outbound_Domestic_Gateway,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_ODG) # from this I choose cp = 0.002

explor.tree_OIG = rpart::rpart(formula = Book_12Mo ~ Outbound_Intr_Gateway,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_OIG) # from this I choose cp = 0.003

explor.tree_RDG = rpart::rpart(formula = Book_12Mo ~ Return_Domestic_Gateway,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_RDG) # from this I choose cp = 0.003

explor.tree_RIG = rpart::rpart(formula = Book_12Mo ~ Return_Intr_Gateway,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_RIG) # from this I choose cp = 0.0018

explor.tree_OCG1 = rpart::rpart(formula = Book_12Mo ~ Outbound_Connect_Gateway1,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_OCG1) # from this I choose cp = 0.001 but may add bin with 33 counts to bin with 2693

explor.tree_OCG2 = rpart::rpart(formula = Book_12Mo ~ Outbound_Connect_Gateway2,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_OCG2) # from this I choose cp = 0.007

explor.tree_RCG1 = rpart::rpart(formula = Book_12Mo ~ Return_Connect_Gateway1,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_RCG1) # from this I choose cp = 0.003

explor.tree_RCG2 = rpart::rpart(formula = Book_12Mo ~ Return_Connect_Gateway2,
                               data = data[split,],
                               method = "anova",
                               control = rpart::rpart.control(cp = 0.0001)) # defalt cp=0.01

summary(explor.tree_RCG2) # from this I choose cp = 0.004

# cp data.frame to prun trees by
cpTree.bin <- data.frame(Variables = c("State", "SourceType", "TourCode", "TourDate_WeekYear",
                                       "Outbound_Domestic_Gateway", "Outbound_Intr_Gateway",
                                       "Return_Domestic_Gateway", "Return_Intr_Gateway",
                                       "Outbound_Connect_Gateway1", "Outbound_Connect_Gateway2",
                                       "Return_Connect_Gateway1", "Return_Connect_Gateway2"),
                         CP = c(0.0019, 0.0002, 0.0015, 0.003,
                                0.002, 0.003,
                                0.003, 0.0018,
                                0.001, 0.007,
                                0.003, 0.004))
# tree.bins() for list of the lookup tables to recategorize inputs
lkup.list <- data[split,] %>%
  select(Book_12Mo, State, SourceType, TourCode, TourDate_WeekYear,
         Outbound_Domestic_Gateway, Outbound_Intr_Gateway,
         Return_Domestic_Gateway, Return_Intr_Gateway,
         Outbound_Connect_Gateway1, Outbound_Connect_Gateway2,
         Return_Connect_Gateway1, Return_Connect_Gateway2) %>%
  tree.bins::tree.bins(y = Book_12Mo, bin.nm = "bin.",
                       method = "anova",
                       control = cpTree.bin,
                       return = "lkup.list")
# lkup.list

# new levels of State
# table(lkup.list[[1]]$State, lkup.list[[1]]$Categories)
# table(lkup.list[[1]]$Categories) # bin.2 is largest
# # new levels of SourceType
# table(lkup.list[[2]]$SourceType, lkup.list[[2]]$Categories)
# table(lkup.list[[2]]$Categories) # bin.3 is largest with 6
# # new levels of TourCode
# table(lkup.list[[3]]$TourCode, lkup.list[[3]]$Categories)
# table(lkup.list[[3]]$Categories) # bin.1 is largest with 40
# # new levels of TourCode
# table(lkup.list[[4]]$TourDate_WeekYear, lkup.list[[4]]$Categories)
# table(lkup.list[[4]]$Categories) # bin.1 is largest with 51
# # new levels of Outbound_Domestic_Gateway
# table(lkup.list[[4]]$Outbound_Domestic_Gateway, lkup.list[[4]]$Categories)
# table(lkup.list[[4]]$Categories) # bin.1 is largest with 77
# # new levels of Outbound_Intr_Gateway
# table(lkup.list[[5]]$Outbound_Intr_Gateway, lkup.list[[5]]$Categories)
# table(lkup.list[[5]]$Categories) # bin.1 is largest with 24
# # new levels of Return_Domestic_Gateway
# table(lkup.list[[6]]$Return_Domestic_Gateway, lkup.list[[6]]$Categories)
# table(lkup.list[[6]]$Categories) # bin.1 is largest with 84
# # new levels of Return_Intr_Gateway
# table(lkup.list[[7]]$Return_Intr_Gateway, lkup.list[[7]]$Categories)
# table(lkup.list[[7]]$Categories) # bin.2 is largest with 33
# # new levels of Outbound_Connect_Gateway1
# table(lkup.list[[8]]$Outbound_Connect_Gateway1, lkup.list[[8]]$Categories)
# table(lkup.list[[8]]$Categories) # bin.1 is largest with 22
# # new levels of Outbound_Connect_Gateway2
# table(lkup.list[[9]]$Outbound_Connect_Gateway2, lkup.list[[9]]$Categories)
# table(lkup.list[[9]]$Categories) # bin.1 is largest with 19
# # new levels of Return_Connect_Gateway1
# table(lkup.list[[10]]$Return_Connect_Gateway1, lkup.list[[10]]$Categories)
# table(lkup.list[[10]]$Categories) # bin.3 is largest with 28
# # new levels of Return_Connect_Gateway2
# table(lkup.list[[11]]$Return_Connect_Gateway2, lkup.list[[11]]$Categories)
# table(lkup.list[[11]]$Categories) # bin.3 is largest with 19

## Bins defined on training data may not define all possible levels
# list unknown State lvls in data
data %>%
  filter(!State %in% unique(lkup.list[[1]]$State)) %>%
  drop_na(State) %>% # do not want to impute na values here
  pull(State) %>%
  unique() %>%
  as.character() -> list.nonlevls_State
list.nonlevls_State # to add to largest bin.2

# update lookup table for State bins:
## with unknown lvls in data maped to largest bin.2
lkup.list[[1]] <- lkup.list[[1]] %>%
  add_row(State = list.nonlevls_State,
          Categories = rep("bin.2",length(list.nonlevls_State)))

# list unknown SourceType lvls in data
data %>%
  filter(!SourceType %in% unique(lkup.list[[2]]$SourceType)) %>%
  drop_na(SourceType) %>% # do not want to impute na values here
  pull(SourceType) %>%
  unique() %>%
  as.character() -> list.nonlevls_SourceType
list.nonlevls_SourceType # empty set

# list unknown TourCode lvls in data
data %>%
  filter(!TourCode %in% unique(lkup.list[[3]]$TourCode)) %>%
  drop_na(TourCode) %>% # do not want to impute na values here
  pull(TourCode) %>%
  unique() %>%
  as.character() -> list.nonlevls_TourCode
list.nonlevls_TourCode # empty set

# list unknown TourDate_WeekYear lvls in data
data %>%
  filter(!TourDate_WeekYear %in% unique(lkup.list[[4]]$TourDate_WeekYear)) %>%
  drop_na(TourDate_WeekYear) %>% # do not want to impute na values here
  pull(TourDate_WeekYear) %>%
  unique() %>%
  as.character() -> list.nonlevls_TDWY
list.nonlevls_TDWY # empty set

# list unknown Outbound_Domestic_Gateway lvls in data
data %>%
  filter(!Outbound_Domestic_Gateway %in% unique(lkup.list[[5]]$Outbound_Domestic_Gateway)) %>%
  drop_na(Outbound_Domestic_Gateway) %>% # do not want to impute na values here
  pull(Outbound_Domestic_Gateway) %>%
  unique() %>%
  as.character() -> list.nonlevls_ODG
list.nonlevls_ODG # to add to largest bin.1

# update lookup table for Outbound_Domestic_Gateway bins:
## with unknown lvls in data maped to largest bin.1
lkup.list[[5]] <- lkup.list[[5]] %>%
  add_row(Outbound_Domestic_Gateway = list.nonlevls_ODG,
          Categories = rep("bin.1",length(list.nonlevls_ODG)))

# check lookup table for unknown Outbound_Domestic_Gateway lvls in dat.score
table(lkup.list[[5]]$Categories) # bin.1 is largest with 87

# list unknown Outbound_Intr_Gateway lvls in data
data %>%
  filter(!Outbound_Intr_Gateway %in% unique(lkup.list[[6]]$Outbound_Intr_Gateway)) %>%
  drop_na(Outbound_Intr_Gateway) %>% # do not want to impute na values here
  pull(Outbound_Intr_Gateway) %>%
  unique() %>%
  as.character() -> list.nonlevls_OIG
list.nonlevls_OIG # to add to largest bin.1

# update lookup table for Outbound_Intr_Gateway bins:
## with unknown lvls in data maped to largest bin.1
lkup.list[[6]] <- lkup.list[[6]] %>%
  add_row(Outbound_Intr_Gateway = list.nonlevls_OIG,
          Categories = rep("bin.1",length(list.nonlevls_OIG)))

# check lookup table for unknown Outbound_Domestic_Gateway lvls in dat.score
table(lkup.list[[6]]$Categories) # bin.1 is largest with 25

# list unknown Return_Domestic_Gateway lvls in data
data %>%
  filter(!Return_Domestic_Gateway %in% unique(lkup.list[[7]]$Return_Domestic_Gateway)) %>%
  drop_na(Return_Domestic_Gateway) %>% # do not want to impute na values here
  pull(Return_Domestic_Gateway) %>%
  unique() %>%
  as.character() -> list.nonlevls_RDG
list.nonlevls_RDG # to add to largest bin.1

# update lookup table for Outbound_Intr_Gateway bins:
## with unknown lvls in data maped to largest bin.1
lkup.list[[7]] <- lkup.list[[7]] %>%
  add_row(Return_Domestic_Gateway = list.nonlevls_RDG,
          Categories = rep("bin.1",length(list.nonlevls_RDG)))

# check lookup table for unknown Outbound_Intr_Gateway lvls in dat.score
table(lkup.list[[7]]$Categories) # bin.1 is largest with 93

# list unknown Return_Intr_Gateway lvls in data
data %>%
  filter(!Return_Intr_Gateway %in% unique(lkup.list[[8]]$Return_Intr_Gateway)) %>%
  drop_na(Return_Intr_Gateway) %>% # do not want to impute na values here
  pull(Return_Intr_Gateway) %>%
  unique() %>%
  as.character() -> list.nonlevls_RIG
list.nonlevls_RIG # to add to largest bin.1

# update lookup table for Return_Intr_Gateway bins:
## with unknown lvls in data maped to largest bin.2
lkup.list[[8]] <- lkup.list[[8]] %>%
  add_row(Return_Intr_Gateway = list.nonlevls_RIG,
          Categories = rep("bin.2",length(list.nonlevls_RIG)))

# check lookup table for unknown Outbound_Connect_Gateway1 lvls in dat.score
table(lkup.list[[8]]$Categories) # bin.2 is largest with 34

# list unknown Outbound_Connect_Gateway1 lvls in data
data %>%
  filter(!Outbound_Connect_Gateway1 %in% unique(lkup.list[[9]]$Outbound_Connect_Gateway1)) %>%
  drop_na(Outbound_Connect_Gateway1) %>% # do not want to impute na values here
  pull(Outbound_Connect_Gateway1) %>%
  unique() %>%
  as.character() -> list.nonlevls_OCG1
list.nonlevls_OCG1 # to add to largest bin.1

# update lookup table for Outbound_Connect_Gateway1 bins:
## with unknown lvls in data maped to largest bin.2
lkup.list[[9]] <- lkup.list[[9]] %>%
  add_row(Outbound_Connect_Gateway1 = list.nonlevls_OCG1,
          Categories = rep("bin.1",length(list.nonlevls_OCG1)))

# check lookup table for unknown Outbound_Connect_Gateway1 lvls in dat.score
table(lkup.list[[9]]$Categories) # bin.1 is largest with 24
#$%
# list unknown Outbound_Connect_Gateway2 lvls in data
data %>%
  filter(!Outbound_Connect_Gateway2 %in% unique(lkup.list[[10]]$Outbound_Connect_Gateway2)) %>%
  drop_na(Outbound_Connect_Gateway2) %>% # do not want to impute na values here
  pull(Outbound_Connect_Gateway2) %>%
  unique() %>%
  as.character() -> list.nonlevls_OCG2
list.nonlevls_OCG2 # to add to largest bin.1

# update lookup table for Outbound_Connect_Gateway2 bins:
## with unknown lvls in data maped to largest bin.1
lkup.list[[10]] <- lkup.list[[10]] %>%
  add_row(Outbound_Connect_Gateway2 = list.nonlevls_OCG2,
          Categories = rep("bin.1",length(list.nonlevls_OCG2)))

# check lookup table for unknown Outbound_Connect_Gateway2 lvls in dat.score
table(lkup.list[[10]]$Categories) # bin.1 is largest with 22

# list unknown Return_Connect_Gateway1 lvls in data
data %>%
  filter(!Return_Connect_Gateway1 %in% unique(lkup.list[[11]]$Return_Connect_Gateway1)) %>%
  drop_na(Return_Connect_Gateway1) %>% # do not want to impute na values here
  pull(Return_Connect_Gateway1) %>%
  unique() %>%
  as.character() -> list.nonlevls_RCG1
list.nonlevls_RCG1 # to add to largest bin.3

# update lookup table for Return_Connect_Gateway1 bins:
## with unknown lvls in data maped to largest bin.1
lkup.list[[11]] <- lkup.list[[11]] %>%
  add_row(Return_Connect_Gateway1 = list.nonlevls_RCG1,
          Categories = rep("bin.3",length(list.nonlevls_RCG1)))

# check lookup table for unknown Return_Connect_Gateway1 lvls in dat.score
table(lkup.list[[11]]$Categories) # bin.3 is largest with 34

# list unknown Return_Connect_Gateway2 lvls in data
data %>%
  filter(!Return_Connect_Gateway2 %in% unique(lkup.list[[12]]$Return_Connect_Gateway2)) %>%
  drop_na(Return_Connect_Gateway2) %>% # do not want to impute na values here
  pull(Return_Connect_Gateway2) %>%
  unique() %>%
  as.character() -> list.nonlevls_RCG2
list.nonlevls_RCG2 # to add to largest bin.3

# update lookup table for Return_Connect_Gateway2 bins:
## with unknown lvls in data maped to largest bin.1
lkup.list[[12]] <- lkup.list[[12]] %>%
  add_row(Return_Connect_Gateway2 = list.nonlevls_RCG2,
          Categories = rep("bin.3",length(list.nonlevls_RCG2)))

# check lookup table for unknown Return_Connect_Gateway2 lvls in dat.score
table(lkup.list[[12]]$Categories) # bin.3 is largest with 34

# save for the report
# write.csv(lkup.list[[1]], "bin_State.csv")
# write.csv(lkup.list[[2]], "bin_SourceType.csv")
# write.csv(lkup.list[[4]], "bin_Outbound_Domestic_Gateway.csv")
# write.csv(lkup.list[[5]], "bin_Outbound_Intr_Gateway.csv")
# write.csv(lkup.list[[6]], "bin_Return_Domestic_Gateway.csv")
# write.csv(lkup.list[[7]], "bin_Return_Intr_Gateway.csv")

# post report 11/23/19
# write.csv(lkup.list[[8]], "bin_Outbound_Connect_Gateway1.csv")
# write.csv(lkup.list[[9]], "bin_Outbound_Connect_Gateway2.csv")
# write.csv(lkup.list[[10]], "bin_Return_Connect_Gateway1.csv")
# write.csv(lkup.list[[11]], "bin_Return_Connect_Gateway2.csv")

# post report 11/28/19
# write.csv(lkup.list[[3]], "bin_TourCode.csv")

#% Pass Data: prep ----
## recategorize inputs from lkup.list with all data as data.prep
data.prep <- tree.bins::bin.oth(list = lkup.list,
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

# summary(data.prep)
# prop.table(table(data.prep[["TourDate_WeekYear"]],data.prep[["FY"]]),2)
# prop.table(table(data.prep[["TourDate_WeekYear"]],data.prep[["Book_12Mo"]]),1)
# prop.table(table(data.prep[["TourCode"]]))
# prop.table(table(data.prep[["TourCode"]],data.prep[["Book_12Mo"]]),1)

#% Set Index ----
# non-inputs
inp.n <- c(which(colnames(data.prep) %in% append(rejected.inputs,"Book_12Mo")))
names(data.prep[inp.n])
# generate num/factor/missing input indices
indx <- inx(data.prep, inp.n)
names(data.prep[indx$num]) # numeric inputs
names(data.prep[indx$cat]) # catagric inputs
names(data.prep[indx$na]) # missing inputs

#%% Save Data: prep ----
prepInput <- list( # for scoring
  key.Grp_Size_Cat=key.Grp_Size_Cat, key.Tour_Region=key.Tour_Region,
  key.Tour_Region=key.Tour_Region, key.TourCode=key.TourCode,
  key.Book_Months=key.Book_Months, key.Age=key.Age,
  key.DB_Enter_Months=key.DB_Enter_Months,
  key.Email=key.Email, key.Past_Trips=key.Past_Trips,
  key.TourPriceCat=key.TourPriceCat, key.Capacity_Cat=key.Capacity_Cat,
  key.TravelAgain=key.TravelAgain,
  lkup.list = lkup.list)
save(inp.n, indx, split, split.valid, split.test,
     prepInput, data.prep,  file="data_prep.RData")

# ## Descriptive Stats ----
# ## Numeric data tables on training data
# data.prep[split,indx$num] %>%
#   psych::describe() %>%
#   select(-c(vars,trimmed)) %>%
#   mutate(na = sum(split) - n,
#          names = names(data.prep[indx$num])) %>%
#   select(names,n,na,everything()) -> traning_numeric_descriptives
#
# # traning_numeric_descriptives
# # write.csv(traning_numeric_descriptives, "traning_numeric_descriptives.csv")
#
# ## Numeric data tables on training data by target
# data.prep[split,indx$num] %>%
#   psych::describeBy(data.prep[split,"Book_12Mo"]) -> traning_numeric_Bydescriptives
#
# ## Target = 0
# traning_numeric_Bydescriptives$`0` %>%
#   select(-c(vars,trimmed)) %>%
#   mutate(na = dim(filter(data.prep[split,], Book_12Mo == 0))[1] - n,
#          names = names(data.prep[indx$num])) %>%
#   select(names,n,na,everything()) -> traning_numeric_descriptives_TargetNOT
#
# traning_numeric_descriptives_TargetNOT
# # write.table(traning_numeric_descriptives_TargetNOT, "traning_numeric_descriptives_TargetNOT.csv")
#
# ## Target = 1
# traning_numeric_Bydescriptives$`1` %>%
#   select(-c(vars,trimmed)) %>%
#   mutate(na = dim(filter(data.prep[split,], Book_12Mo == 1))[1] - n,
#          names = names(data.prep[indx$num])) %>%
#   select(names,n,na,everything()) -> traning_numeric_descriptives_Target
#
# traning_numeric_descriptives_Target
# # write.table(traning_numeric_descriptives_Target, "traning_numeric_descriptives_Target.csv")
#
# ## Variance Importance Measure ----
# # nominal input
# chi2p.val<- sapply(data.prep[split,indx$cat],
#                      function(x) chisq.test(x, data.prep[split,"Book_12Mo"])$p.value)
# VIM.cat_1 <- data.frame(chi2p.val = data.frame(chi2p.val)[order(chi2p.val),,drop = F],
#                         rank.chi2 = rank(data.frame(chi2p.val)[order(chi2p.val),]))
#
# chi2p.sim <- sapply(data.prep[split,indx$cat],
#                      function(x) chisq.test(x, data.prep[split,"Book_12Mo"],simulate.p.value = TRUE)$p.value)
# VIM.cat_2 <- data.frame(chi2p.sim = data.frame(chi2p.sim)[order(chi2p.sim),,drop = F],
#                         rank.chi2sim = rank(data.frame(chi2p.sim)[order(chi2p.sim),]))
#
# ROC.cat <- caret::filterVarImp(x = data.prep[split,indx$cat],
#                                y = data.prep[split,"Book_12Mo"])
# VIM.cat_3 <- data.frame(ROC.cat = ROC.cat[order(-ROC.cat$X1),]$X1,
#                         rank.ROC = rank(-ROC.cat[order(-ROC.cat$X1),]$X1),
#                         row.names = row.names(ROC.cat[order(-ROC.cat$X1),,drop = FALSE]))
# VIM.cat <- VIM.cat_1 %>%
#   merge(., VIM.cat_2, by="row.names", all=TRUE, sort = FALSE) %>%
#   merge(., VIM.cat_3, by.x="Row.names", by.y="row.names", all=TRUE, sort = FALSE)
# # write.csv(VIM.cat, "VIMs_cat.csv")
#
# # numeric input
# tTestp.val<- sapply(data.prep[split,indx$num],
#                     function(x) t.test(x ~ data.prep[split,"Book_12Mo"])$p.value)
# VIM.num_1 <- data.frame(tTestp.val = data.frame(tTestp.val)[order(tTestp.val),,drop = F],
#                         rank.t = rank(data.frame(tTestp.val)[order(tTestp.val),]))
#
# ROC.num <- caret::filterVarImp(x = data.prep[split,indx$num],
#                                y = data.prep[split,"Book_12Mo"])
# VIM.num_2 <- data.frame(ROC.num = ROC.num[order(-ROC.num$X1),]$X1,
#                         rank.ROC = rank(-ROC.num[order(-ROC.num$X1),]$X1),
#                         row.names = row.names(ROC.num[order(-ROC.num$X1),,drop = FALSE]))
# VIM.num <- VIM.num_1 %>%
#   merge(., VIM.num_2, by="row.names", all=TRUE, sort = FALSE)
# # write.csv(VIM.num, "VIMs_num.csv")
#
# # all inputs
# ROC_all <- caret::filterVarImp(x = data.prep[split,-inp.n],
#                                y = data.prep[split,"Book_12Mo"])
#
#
# ROC.all <- data.frame(ROC = ROC_all[order(-ROC_all$X1),]$X1,
#                       rank.ROC = rank(-ROC_all[order(-ROC_all$X1),]$X1),
#                       row.names = row.names(ROC_all[order(-ROC_all$X1),,drop = FALSE]))
# # write.csv(ROC.all, "ROC_all.csv")

# Imputation ----
#%%%%%%%%%%%%%%%%%%%%%%#
#% check missing data %#
#%%%%%%%%%%%%%%%%%%%%%%#
# names(data.prep)[indx$na] # names of inputs with missing data
# any(indx$na & indx$cat)   # any cat inputs with missing?
# any(indx$na & indx$num)   # any num inputs with missing?
# names(data.prep)[indx$na & indx$num] # names of num inputs with missing data
# names(data.prep)[indx$na & indx$cat] # names of cat inputs with missing data

imp.parms.num <- data.prep %>%                       # saving imputation parameters
  summarise_at(names(data.prep[indx$num]),           #   for num inputs
               ~median(.[split], na.rm = TRUE))      #     save median
imp.parms.cat <- data.prep %>%                       # saving imputation parameters
  summarise_at(names(data.prep[indx$cat]),           #   for cat inputs
               ~z.mode(.[split]))                    #     save mode

#% Pass Data: imp ----
data.imp <- data.prep %>%                                      # Routine: pass data
  mutate_at(names(data.prep[indx$na]),                         # for vars with missings
            list(na = ~ as.factor(ifelse(is.na(.),1,0)))) %>%  #  create na flag as nominal inputs
  mutate_at(names(data.prep[indx$num]),                        # for numaric vars
            ~ impute(x = .,
                     y = median(.[split], na.rm = TRUE))) %>%  #  impute mean on training data
  mutate_at(names(data.prep[indx$cat]),                        # for cat vars
            ~ impute(x = .,
                     y = z.mode(.[split])))                    #  impute first mode on training data

#%% Save Data: imp ----
#%%%%%%%%%%%%%%%%%%%%%%%%#
#% Routine: indx update %#
#%%%%%%%%%%%%%%%%%%%%%%%%#
# non-inputs
inp.n <- c(which(colnames(data.imp) %in% append(rejected.inputs,"Book_12Mo")))
names(data.imp[inp.n])
# generate num/factor/missing input indices
indx <- inx(data.imp, inp.n)
# names(data.imp[indx$num]) # numeric inputs
# names(data.imp[indx$cat]) # catagric inputs
# names(data.imp[indx$na]) # missing inputs

save(inp.n, indx, split, split.valid, split.test,
     imp.parms.num, imp.parms.cat, data.imp,  file="data_imp.RData")

## Down Sampling ----
vars <- -c(which(colnames(data.prep) %in% append(rejected.inputs[-1],"Book_12Mo"))) # all varables not target or rejected inputs


set.seed(77012)                                                       # reset seed
downSampledTrain <- caret::downSample(x = data.prep[split, vars],     # down sample traning data
                                      y = data.prep$Book_12Mo[split], # target variable with the class membership
                                      yname = "Book_12Mo")            # preserve target name
#%%%%%%%%%%%%%%%%%%%%%#
#% check down sample %#
#%%%%%%%%%%%%%%%%%%%%%#
# table(data.prep[["Book_12Mo"]][split]) # summary of target before downsample
# table(downSampledTrain$Book_12Mo)      # summary of downsample target
# prop.table(table(data.prep[["Book_12Mo"]][split])) # about 70:20 target class unbalence in traing data
# prop.table(table(downSampledTrain[["Book_12Mo"]])) # 50:50 target class balence in downsample traning data
# any(duplicated(data.clean$EvalID)) # Is there any duplicatged EvalID in the data? FALSE

#% Pass Data: down ----
#% Combine downsampled train with validation %#
data.down <- data.prep[(data.prep$EvalID %in% downSampledTrain$EvalID) | (!split),]
# This is the downsampled training data, 25% validation data, and 25% testing data #

#% Indicator for balanced train data %#
split.down <- ifelse(data.down$EvalID %in% downSampledTrain$EvalID, TRUE, FALSE)

#% Indicator for valid and test data in downsample data %#
split.down.valid <- !split.down
split.down.test <- !split.down
split.down.valid[split.down.valid == TRUE] <- split2
split.down.test[split.down.test == TRUE] <- !split2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#% check valid and test target in downsample %#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# table(data.prep[["Book_12Mo"]][split.test],useNA = "always")
# table(data.down[["Book_12Mo"]][split.down.test],useNA = "always")
# table(data.prep[["Book_12Mo"]][split.valid],useNA = "always")
# table(data.down[["Book_12Mo"]][split.down.valid],useNA = "always")

#%% Save Data: down ----
#%%%%%%%%%%%%%%%%%%%%%%%%#
#% Routine: indx update %#
#%%%%%%%%%%%%%%%%%%%%%%%%#
# non-inputs
inp.n <- c(which(colnames(data.down) %in% append(rejected.inputs,"Book_12Mo")))
names(data.down[inp.n])
# generate num/factor/missing input indices
indx <- inx(data.down, inp.n)
# names(data.down[indx$num]) # numeric inputs
# names(data.down[indx$cat]) # catagric inputs
# names(data.down[indx$na]) # missing inputs

save(inp.n, indx, split.down, split.down.valid, split.down.test,
     data.down, file="data_down.RData")

#% Imputation: data.down ----
#%%%%%%%%%%%%%%%%%%%%%%#
#% check missing data %#
#%%%%%%%%%%%%%%%%%%%%%%#
# names(data.down)[indx$na] # names of inputs with missing data
# any(indx$na & indx$cat)   # any cat inputs with missing?
# any(indx$na & indx$num)   # any num inputs with missing?
# names(data.down)[indx$na & indx$num] # names of num inputs with missing data
# names(data.down)[indx$na & indx$cat] # names of cat inputs with missing data


down.imp.parms.num <- data.down %>%                  # saving imputation parameters
  summarise_at(names(data.down[indx$num]),           #   for num inputs
               ~median(.[split.down], na.rm = TRUE)) #     save median
down.imp.parms.cat <- data.down %>%                  # saving imputation parameters
  summarise_at(names(data.down[indx$cat]),           #   for cat inputs
               ~z.mode(.[split.down]))               #     save first mode

#% Pass Data: down.imp ----
data.down.imp <- data.down %>%                                 # Routine: pass data
  mutate_at(names(data.down[indx$na]),                         # for vars with missings
            list(na = ~ as.factor(ifelse(is.na(.),1,0)))) %>%  #  create na flag as nominal inputs
  mutate_at(names(data.down[indx$num]),                        # for numaric vars
            ~ impute(x = .,                                    #  impute median of downsampled training data
                     y = median(.[split.down],na.rm=TRUE))) %>%
  mutate_at(names(data.down[indx$cat]),                        # for numaric vars
            ~ impute(x = .,                                    #  impute mode of downsampled training data
                     y = z.mode(.[split.down])))


#%% Save Data: down.imp ----
#% Routine: indx update %#
# non-inputs
inp.n <- c(which(colnames(data.down.imp) %in% append(rejected.inputs,"Book_12Mo")))
names(data.down.imp[inp.n])
# generate num/factor/missing input indices
indx <- inx(data.down.imp, inp.n)
# names(data.down.imp[indx$num]) # numeric inputs
# names(data.down.imp[indx$cat]) # catagric inputs
# names(data.down.imp[indx$na]) # missing inputs

save(inp.n, indx, split.down, split.down.valid, split.down.test,
     down.imp.parms.num, down.imp.parms.cat, data.down.imp,
     file="data_down_imp.RData")


# distance <- get_dist(data.prep[inp.n])
# factoextra::fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
