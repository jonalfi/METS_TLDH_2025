############################################################################
#  IMPROVING CARDIOPULMONARY EXERCISE TESTING GUIDED PREOPERATIVE RISK
#   PREDICTION IN NON-CARDIAC SURGERY BY ADRESSING BODY-SPECIFIC AND
#     SEX-RELATED DIFFERENCES THROUGH A MACHINE LEARNING MODEL.
#          A SECONDARY POST HOC ANALYSIS OF THE METS STUDY.
#
#                         CODE WRITTEN BY
#                       DR. JONAS ALFITIAN
#     DEPARTMENT FOR ANESTHESIOLOGY AND INTENSIVE CARE MEDICINE,
#           FACULTY OF MEDICINE, UNIVERSITY OF COLOGNE
#
#
# ORIGINAL PATIENT DATA ORIGINATES FROM THE METS STUDY
# AND WILL NOT BE PUBLICLY SHARED
#############################################################################

# ===========================
# LOAD PACKAGES AND SOURCES
# ===========================

library(tidyverse)
library(readstata13)
library(lubridate)
library(boot)
library(caret)
library(shapviz)
library(iml)

source("Functions/spiegelhalter_z.R")
source("Functions/NRIcont.R")
source("Functions/NRIcont_boot.R")

# Original data set will not be shared due to data protection regulations
import_rawdataset <- read.dta13("Data/METS - Export Dataset.dta")

import_rawdataset$cpet.date <- ydm(paste(import_rawdataset$CPET_DT_YYYY, import_rawdataset$CPET_DT_DD, import_rawdataset$CPET_DT_MM))
import_rawdataset$surgery.date <- ydm(paste(import_rawdataset$SURGERY_DT_YYYY, import_rawdataset$SURGERY_DT_DD, import_rawdataset$SURGERY_DT_MM))

# Exclude all cases with incomplete CPET data or missing data on surgery details
import_rawdataset <- import_rawdataset %>%
  filter(!is.na(CPET_DT_YYYY)) %>%
  filter(!is.na(SURGERY_DT_YYYY)) %>%
  filter(surgery.date > cpet.date) %>%
  filter(!is.na(OXYGEN_AT_VO2_PEAK))  %>%
  filter(!is.na(OXYGEN_CONSUM_AT)) %>%
  filter(!is.na(VE_VCO2_AT)) %>%
  filter(!is.na(VE_VO2_AT)) %>%
  filter(!is.na(COMPLICATIONS_NY_STD))

age <- import_rawdataset$AGE
sex <- as.factor(ifelse(import_rawdataset$SEX_STD == 1, "Male", "Female"))
height <- import_rawdataset$HEIGHT1
weight <- import_rawdataset$WEIGHT1

# Binary classification of high-risk surgery as described in the manuscript
proc.type <- as.factor(import_rawdataset$procedure_simple)
highrisk.surgery <- ifelse((import_rawdataset$procedure_simple == "Intra-peritoneal or retro-peritonal" | import_rawdataset$procedure_simple == "Vascular" | import_rawdataset$procedure_simple == "Intra-thoracic"), 1, 0)
highrisk.surgery.fct <- as.factor(ifelse(highrisk.surgery == TRUE, "High-risk", "Non-High-risk"))

# CPET Parameters
# =====================
# Absolute Peak VO2
peakVO2 <- import_rawdataset$OXYGEN_AT_VO2_PEAK * weight
# Absolute Peak VO2
peakVO2.bw <- peakVO2 / weight
# Classification of reduced peak VO2 as described in the manuscript
low.peakVO2.bw <- ifelse(import_rawdataset$OXYGEN_AT_VO2_PEAK < 14, 1, 0)

# Outcome Variables
# =====================
## Complications
complications.stay <- as.factor(import_rawdataset$COMPLICATIONS_NY_STD)
levels(complications.stay) <- c("None", "Mild", "Moderate", "Severe", "Fatal")
complications.dicho <- (complications.stay == "Moderate" | complications.stay == "Severe" | complications.stay == "Fatal")
complications.dicho <- ifelse(complications.dicho, 1, 0)
complications.dicho.fct <- as.factor(complications.dicho)
levels(complications.dicho.fct) <- c("Non.Severe", "Severe")

# Data frame with data relevant for model training
wrangled.data <- data.frame(age, sex, height, weight, proc.type, highrisk.surgery, highrisk.surgery.fct,
                            peakVO2, peakVO2.bw, low.peakVO2.bw,
                            complications.dicho, complications.dicho.fct
)

# Splitting data into training and testing sets
set.seed(125)
caseIDs <- createDataPartition(wrangled.data$complications.dicho, p = 0.7, list = FALSE)
df.training <- wrangled.data[caseIDs, ] %>%
  select(complications.dicho.fct, complications.dicho, age, highrisk.surgery, peakVO2, peakVO2.bw, low.peakVO2.bw, sex, height, weight)

df.test <- wrangled.data[-caseIDs, ] %>%
  select(complications.dicho.fct, complications.dicho, age, highrisk.surgery, peakVO2, peakVO2.bw, low.peakVO2.bw, sex, height, weight)

df.full <- rbind(df.test, df.training)

Y_training <- df.training$complications.dicho.fct
Y_test <- df.test$complications.dicho.fct
Y_full <- df.full$complications.dicho.fct

dmy <- dummyVars(" ~ .", data = df.training)
df.training <- data.frame(predict(dmy, newdata = df.training))

dmy <- dummyVars(" ~ .", data = df.test)
df.test <- data.frame(predict(dmy, newdata = df.test))

dmy <- dummyVars(" ~ .", data = df.full)
df.full <- data.frame(predict(dmy, newdata = df.full))

# Assign model training parameters
# 10-times repeated 5-fold cross-validation with SMOTE resampling due to class imbalance
fit_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  summaryFunction = twoClassSummary,
  sampling = "smote",
  classProbs = TRUE
)

# ===========================
# MODEL TRAINING
# ===========================

# Reference Model (generalized linear logistic regression model (GLM) 
set.seed(125)
cv.glm.ctrl <- caret::train(x = select(df.training, peakVO2.bw), y = Y_training, method = "glm", trControl = fit_control)
probabilities.ctrl <- predict(cv.glm.ctrl, newdata = select(df.test, peakVO2.bw), type = "prob")[, "Severe"]
predictions.ctrl <- ROCR::prediction(probabilities.ctrl, Y_test)
performance.ctrl <- ROCR::performance(predictions.ctrl, "auc")
roc.ctrl <- pROC::roc(Y_test, probabilities.ctrl)

# AbsPeakVO2 Model (gradient boosting machine (GBM))
set.seed(125)
cv.gbm.peakVO2 <- caret::train(x = select(df.training, age, highrisk.surgery, peakVO2), y = Y_training,
                               verbose = FALSE,
                               method = "gbm",
                               metric = "ROC",
                               preProcess = c('center', 'scale'),
                               trControl = fit_control
)
probabilities.peakVO2 <- predict(cv.gbm.peakVO2, newdata = select(df.test, age, highrisk.surgery, peakVO2), type = "prob")[, "Severe"]
predictions.peakVO2 <- ROCR::prediction(probabilities.peakVO2, Y_test)
performance.peakVO2 <- ROCR::performance(predictions.peakVO2, "auc")
roc.peakVO2 <- pROC::roc(Y_test, probabilities.peakVO2)

# NEWPeakVO2 Model (gradient boosting machine (GBM))
set.seed(125) 
cv.gbm.full <- caret::train(x = select(df.training, age, highrisk.surgery, peakVO2, weight, height, sex.Male), y = Y_training,
                            verbose = FALSE,
                            method = "gbm",
                            metric = "ROC",
                            preProcess = c('center', 'scale'),
                            trControl = fit_control
)
probabilities.full <- predict(cv.gbm.full, newdata = select(df.test, age, highrisk.surgery, peakVO2, weight, height, sex.Male), type = "prob")[, "Severe"]
predictions.full <- ROCR::prediction(probabilities.full, Y_test)
performance.full <- ROCR::performance(predictions.full, "auc")
roc.full <- pROC::roc(Y_test, probabilities.full)

# ===========================
# MODEL EVALUATION
# ===========================

#Model Calibration assessed by Brier-Score and Spiegelhalter-Z statistic
#Y_test is a factor and Brier and Spiegelhalter_z functions take numeric vectors as arguments; factor levels start with 1, thus -1 for binary classification
ModelMetrics::brier(as.numeric(Y_test)-1, probabilities.ctrl)
ModelMetrics::brier(as.numeric(Y_test)-1, probabilities.peakVO2)
ModelMetrics::brier(as.numeric(Y_test)-1, probabilities.full)

spiegelhalter_z(as.numeric(Y_test)-1, probabilities.ctrl)
spiegelhalter_z(as.numeric(Y_test)-1, probabilities.peakVO2)
spiegelhalter_z(as.numeric(Y_test)-1, probabilities.full)

#Model discrimination assessed by area under the ROC curve 
pROC::auc(roc.ctrl)
pROC::auc(roc.peakVO2)
pROC::auc(roc.full)

# Model evaluation through computing net reclassification improvement index (NRI)
# Predicted probabilities from the corresponding models are compared for each observation, such that for each observation the direction of change of the predicted probability is computed
# Results are stored in a data frame taken as an argument by the NRI function
nri.df.ctrl_full <- df.test %>%
  select(complications.dicho.fct.Non.Severe, complications.dicho.fct.Severe) %>%
  mutate(id = row_number()) %>%
  bind_cols(fitted.old = probabilities.ctrl, fitted.new = probabilities.full) %>%
  mutate(up = ifelse(fitted.new > fitted.old, TRUE, FALSE)) %>%
  mutate(down = ifelse(fitted.new < fitted.old, TRUE, FALSE)) %>%
  mutate(change = ifelse(fitted.new > fitted.old, "Up", "Down"))

nri.df.peakVO2_full <- df.test %>%
  select(complications.dicho.fct.Non.Severe, complications.dicho.fct.Severe) %>%
  mutate(id = row_number()) %>%
  bind_cols(fitted.old = probabilities.peakVO2, fitted.new = probabilities.full) %>%
  mutate(up = ifelse(fitted.new > fitted.old, TRUE, FALSE)) %>%
  mutate(down = ifelse(fitted.new < fitted.old, TRUE, FALSE)) %>%
  mutate(change = ifelse(fitted.new > fitted.old, "Up", "Down"))

# 500 replicates for bootstrapping approach
boot.R <- 500
#Bootstrapped NRI statistic is returned as a vector containing NRI, event NRI, non-event NRI, 95% confidence interval margins and p-value
nri.ctrl_full <- NRI.cont.boot(nri.df.ctrl_full, boot.R)
nri.peakVO2_full <- NRI.cont.boot(nri.df.peakVO2_full, boot.R)

# MODEL INTERPRETABILITY
# ===========================

#### SHAP VALUES
# Prediction wrapper function for returning vectorized predicted probabilities from as required for computing SHAP values by explain function from fastshap package
pfun <- function(object, newdata) {  
  unname(predict(object, newdata = newdata, type = "prob")[, "Severe"])
}

#Training data with features is extracted for explain function from fastshap package for each model
X_train.full <- df.training %>%
  select(age, highrisk.surgery, peakVO2, weight, height, sex.Male)

X_train.peakVO2 <- df.training %>%
  select(age, highrisk.surgery, peakVO2)

X_train.ctrl <- df.training %>%
  select(peakVO2.bw)

# SHAP objects are computed by Monte-Carlo simulations with 100 repetitions
set.seed(2113)
sh.full <- fastshap::explain(cv.gbm.full, X = X_train.full, pred_wrapper = pfun, nsim = 100, shap_only = FALSE)
sh.peakVO2 <- fastshap::explain(cv.gbm.peakVO2, X = X_train.peakVO2, pred_wrapper = pfun, nsim = 100, shap_only = FALSE)
sh.ctrl <- fastshap::explain(cv.glm.ctrl, X = X_train.ctrl, pred_wrapper = pfun, nsim = 100, shap_only = FALSE)

#### Interaction effects with Friedman's H-statistic computed and conditional inference tree computed
# Predictor R6-object is created to enable quantification of interaction effects 
predictor <- Predictor$new(cv.gbm.full, data = X_train.full, y = Y_training, type = "prob", class = "Severe")
interact.all <- Interaction$new(predictor)
interact.peakVO2 <- Interaction$new(predictor, feature = "peakVO2")
# cif contains a conditional inference tree, where each possible path is represented as a separate string
cif <- TreeSurrogate$new(predictor, maxdepth = 4)

