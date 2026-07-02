# load required packages
library(tidyverse) # version 2.0.0
library(pROC) # version 1.19.0.1
library(caret) # version 7.0-1
# R version 4.4.2


# load data
# The data needs to be in the following format:
# For each patient, the dataset contains several rows; one row per measurement. The time of the
# respective measurement is included as separate variable. For each measurement, the state
# probabilities are included. Further, the start of the therapy administration (in minutes) and the time
# point of evaluation (in minutes) are available for each row. The corresponding AST result is matched to each
# therapy. The dataset is reduced to those values which were measured after the start of therapy administration
# and before the time of therapy reassessment.


####################################################
## Approach (i): Averaging of State Probabilities ##
####################################################

# calculate mean across all probabilities for state one
avg_states <- data %>% 
  group_by(id, drug_name) %>% 
  mutate(est_eff = mean(local_prob_state1))

# calculate AUC, accuracy, sensitivity, specificity, precision, F1-score
avg_states$pred_class <- as.factor(ifelse(avg_states$est_eff >= 0.5, "S", "R"))
confusion_matrix <- table(True = avg_states$susceptibility, Pred = avg_states$pred_class)

TP <- confusion_matrix["S", "S"]
FN <- confusion_matrix["S", "R"]
FP <- confusion_matrix["R", "S"]
TN <- confusion_matrix["R", "R"]

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
precision   <- TP / (TP + FP)
accuracy    <- (TP + TN) / sum(confusion_matrix)
f1          <- 2 * precision * sensitivity / (precision + sensitivity)

roc_obj <- roc(
  response = avg_states$susceptibility,
  predictor = avg_states$est_eff,
  levels = c("R", "S"),
  direction = "<"
)

auc(roc_obj)



##############################################
## Approach (ii): Logistic Regression Model ##
##############################################

# calculate the difference between the probability for state one at the start of therapy administration and the probability for the first state at the time of reassessment
data_LR <- data %>% 
  group_by(id, drug_name) %>% 
  slice(c(1,n())) %>% 
  mutate(count = c(1,2))

# get IDs of patients with a probability available at therapy start and at time of reassessment
IDs <- filter(data_LR, count == 2)
IDs <- IDs$id

# filter data according to stay IDs
data_LR <- data_LR %>%
  subset(id %in% IDs)

# derive differences used as independent variables in LRM
data_LR <- data_LR %>% 
  select(id, local_prob_state1, drug_name, susceptibility) %>%
  group_by(id, drug_name, susceptibility) %>% 
  mutate(diff_prob1 = local_prob_state1 - local_prob_state1[1]) %>%
  slice(n())




## Model Fitting
set.seed(65421)
nfolds <- 10
split <- createFolds(data_LR$susceptibility, k = nfolds)

total_performance <- matrix(nrow = nfolds, ncol = 6)
colnames(total_performance) <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1")

# model fitting per fold
for (i in 1:nfolds) {
  
  rows <- unlist(split[i])
  
  train <- data_LR[-rows,]
  test <- data_LR[rows,]
  test <- ungroup(test)
  
  mean <- mean(c(sum(train$susceptibility == "S"), sum(train$susceptibility == "R")))
  
  # up- and downsampling
  up <- case_when(sum(train$susceptibility == "S") >= mean & sum(train$susceptibility == "R") >= mean ~ "S, R",
                  sum(train$susceptibility == "S") >= mean & sum(train$susceptibility == "R") < mean ~ "S",
                  sum(train$susceptibility == "S") < mean & sum(train$susceptibility == "R") >= mean ~ "R")
  up <- unlist(str_split(up, ", "))
  
  up_data <- train %>% 
    filter(susceptibility %in% up) %>% 
    group_by(susceptibility) %>% 
    sample_n(size = mean, replace = F)
  
  down_data <- train %>% 
    filter(! susceptibility %in% up) %>% 
    group_by(susceptibility) %>% 
    sample_n(size = mean, replace = T)
  
  train_balanced <- rbind(down_data, up_data)
  train_balanced$susceptibility_num <- ifelse(train_balanced$susceptibility == "S", 1, 0)
  
  # define inputs
  X <- subset(test, select = c(diff_prob1))
  
  # fit the model
  logistic <- glm(susceptibility_num ~ diff_prob1, family = binomial(link='logit'), data = train_balanced)
  
  # make predictions
  y_scores <- predict(logistic, newdata = X, type='response')
  y_scores_class <- ifelse(y_scores > 0.5, "S", "R")
  
  # create confusion matrix
  conf <- confusionMatrix(as.factor(y_scores_class), test$susceptibility, positive = "S")
  
  # calculate AUC
  pred <- prediction(y_scores, test$susceptibility)
  auc <- performance(pred, measure = "auc")
  auc <- auc@y.values[[1]]
  
  # save results
  total_performance[i,] <- c(auc, # AUC by class
                             conf$overall[1], # accuracy
                             conf$byClass[c(1:2, 5, 7)]) # sensitivity, specificity, precision, f1
}

total_performance <- as.data.frame(total_performance)
auc_avg <- mean(total_performance$AUC)
acc_avg <- mean(total_performance$Accuracy)
sensitivity_avg <- mean(total_performance$Sensitivity)
specificity_avg <- mean(total_performance$Specificity)
precision_avg <- mean(total_performance$Precision)
f1_avg <- mean(total_performance$F1)
total_performance_avg <- as.data.frame(matrix(c(auc_avg, acc_avg, sensitivity_avg, specificity_avg, precision_avg, f1_avg), nrow = 1, ncol = 6))
colnames(total_performance_avg) <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1")
total_performance_avg




## Logistic Regression Model with Reduced Data ##
# filter differences >= 0.3
data_30 <- data_LR %>% 
  filter(abs(diff_prob1) >= 0.3)

#### split data into training and test data
set.seed(65421)
nfolds <- 10
split <- createFolds(data_30$susceptibility, k = nfolds)

total_performance_red <- matrix(nrow = nfolds, ncol = 6)
colnames(total_performance_red) <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1")

# model fitting per fold
for (i in 1:nfolds) {
  
  rows <- unlist(split[i])
  
  train <- data_30[-rows,]
  test <- data_30[rows,]
  test <- ungroup(test)
  
  mean <- mean(c(sum(train$susceptibility == "S"), sum(train$susceptibility == "R")))
  
  # up- and downsampling
  up <- case_when(sum(train$susceptibility == "S") >= mean & sum(train$susceptibility == "R") >= mean ~ "S, R",
                  sum(train$susceptibility == "S") >= mean & sum(train$susceptibility == "R") < mean ~ "S",
                  sum(train$susceptibility == "S") < mean & sum(train$susceptibility == "R") >= mean ~ "R")
  up <- unlist(str_split(up, ", "))
  
  up_data <- train %>% 
    filter(susceptibility %in% up) %>% 
    group_by(susceptibility) %>% 
    sample_n(size = mean, replace = F)
  
  down_data <- train %>% 
    filter(! susceptibility %in% up) %>% 
    group_by(susceptibility) %>% 
    sample_n(size = mean, replace = T)
  
  train_balanced <- rbind(down_data, up_data)
  train_balanced$susceptibility_num <- ifelse(train_balanced$susceptibility == "S", 1, 0)
  
  # define inputs
  X <- subset(test, select = c(diff_prob1))
  
  # fit the model
  logistic <- glm(susceptibility_num ~ diff_prob1, family = binomial(link='logit'), data = train_balanced)
  
  # make predictions
  y_scores <- predict(logistic, newdata = X, type='response')
  y_scores_class <- ifelse(y_scores > 0.5, "S", "R")
  
  # create confusion matrix
  conf <- confusionMatrix(as.factor(y_scores_class), test$susceptibility, positive = "S")
  
  # calculate AUC
  pred <- prediction(y_scores, test$susceptibility)
  auc <- performance(pred, measure = "auc")
  auc <- auc@y.values[[1]]
  
  # save results
  total_performance_red[i,] <- c(auc, # AUC by class
                                 conf$overall[1], # accuracy
                                 conf$byClass[c(1:2, 5, 7)]) # sensitivity, specificity, precision, f1
}

total_performance_red <- as.data.frame(total_performance_red)
auc_avg <- mean(total_performance_red$AUC)
acc_avg <- mean(total_performance_red$Accuracy)
sensitivity_avg <- mean(total_performance_red$Sensitivity)
specificity_avg <- mean(total_performance_red$Specificity)
precision_avg <- mean(total_performance_red$Precision)
f1_avg <- mean(total_performance_red$F1)
total_performance_red_avg <- as.data.frame(matrix(c(auc_avg, acc_avg, sensitivity_avg, specificity_avg, precision_avg, f1_avg),
                                                  nrow = 1, ncol = 6))
colnames(total_performance_red_avg) <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1")
total_performance_red_avg




