## Libraries
library(haven)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(glmnet)
library(randomForest)

## Data
# Demographics
demo <- read_xpt("DEMO_L.xpt")
# Dietary Interview
diet1 <- read_xpt("DR1IFF_L.xpt")
diet2 <- read_xpt("DR2IFF_L.xpt")
# Food Security Questionarre
food_sec <- read_xpt("FSQ_L.xpt")


## Data Cleaning
# DDS calculation
# Obtain first digit of food code
diet1$DR1IFDCD <- format(diet1$DR1IFDCD,scientific = F)
diet2$DR2IFDCD <- format(diet2$DR2IFDCD,scientific = F)
diet1$food_group <- as.numeric(substr(diet1$DR1IFDCD,1,1))
diet2$food_group <- as.numeric(substr(diet2$DR2IFDCD,1,1))
# Count by number of unique food groups as DDS for both days
diet_dds1 <- diet1 %>% group_by(SEQN) %>%
  summarise(count1 = length(unique(food_group)))
diet_dds2 <- diet2 %>% group_by(SEQN) %>%
  summarise(count2 = length(unique(food_group)))
# Filter participants with complete data
common_ids <- Reduce(intersect, list(diet_dds2$SEQN,diet_dds1$SEQN, demo$SEQN, food_sec$SEQN))
# Combine datasets
data_combined <- demo %>% filter(SEQN %in% common_ids) %>%
  left_join(diet_dds1, by = "SEQN") %>%
  left_join(diet_dds2, by = "SEQN") %>%
  left_join(food_sec, by = "SEQN")
data_unfactored <- data_combined %>% 
  # Caculate day1 and day2 DDS
  mutate(DDS = (count1+count2)/2)%>%
  # Variables of interest
  dplyr::select(DMDEDUC2,RIAGENDR,RIDAGEYR,DMDBORN4,DMDMARTZ,
                DMDHHSIZ,INDFMPIR,RIDRETH3,FSDAD,DDS,SEQN)%>%
  # Remove missing data
  filter(!is.na(DMDEDUC2),!is.na(INDFMPIR),!is.na(FSDAD),!is.na(DMDMARTZ),DMDMARTZ <= 3)
# Factor and scale variables
data_final <- data_unfactored %>%
  mutate(
    DMDEDUC2 = factor(DMDEDUC2,
                      levels = c(1, 2, 3, 4, 5),
                      labels = c("Less than 9th grade", "9–11th grade",
                                 "High school graduate/GED", "Some college or AA",
                                 "College graduate+")),
    
    FSDAD = factor(FSDAD,
                   levels = c(1, 2, 3, 4),
                   labels = c("Full food security", "Marginal", "Low", "Very low")),
    
    RIAGENDR = factor(RIAGENDR,
                      levels = c(1, 2),
                      labels = c("Male", "Female")),
    
    DMDBORN4 = factor(DMDBORN4,
                      levels = c(1, 2),
                      labels = c("US", "Non-US")),
    
    DMDMARTZ = factor(DMDMARTZ,
                      levels = c(1, 2, 3),
                      labels = c("Married/Partner","Widowed/Divorced/Separated",
                                 "Never married")),
    
    RIDRETH3 = factor(RIDRETH3,
                      levels = c(1, 2, 3, 4, 6, 7),
                      labels = c("Mexican American","Other Hispanic",
                                 "Non-Hispanic White","Non-Hispanic Black",
                                 "Non-Hispanic Asian","Other / Multi-racial")),
    RIDAGEYR_s = scale(RIDAGEYR),
    INDFMPIR_s = scale(INDFMPIR),
    DMDHHSIZ_s = scale(DMDHHSIZ)
  )

## Data Summary
# Age and gender distribution
data_final %>%
  mutate(age_group = cut(
    RIDAGEYR,breaks = c(18,28,38,48,58,68,78,81),
    right = FALSE,
    labels = c("18-27", "28-37", "38-47",
               "48-57", "58-67", "68-77","78+"))) %>%
  group_by(age_group, RIAGENDR) %>%
  summarise(count = n(),.groups = "drop") %>%
  ggplot(aes(y= age_group,
             x = ifelse(RIAGENDR == "Male", -count, count),
             fill = RIAGENDR)) +
  geom_col()+
  scale_x_continuous(labels = abs) +
  labs(
    x = "Number of Individuals",
    y = "Age",
    fill = "Gender",
    title = "Age Distribution by Gender") 
# Race Distribution
data_final %>% group_by(RIDRETH3) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))%>%
  mutate(RIDRETH3 = factor(RIDRETH3, levels = RIDRETH3),
         percent = round(count /nrow(data_final) * 100, 1)) %>%
  ggplot(aes(x="",y = count, fill = RIDRETH3)) +
  geom_col() +
  coord_polar(theta = "y",start = 90) +
  geom_text(aes(label = paste0(percent, "%")),
            position = position_stack(vjust = 0.5),
            size = 4) +
  labs(title = "Race Distribution",fill = "Race")+
  theme_void()
# DDS Distribution with mean
hist(data_final$DDS, main = "Distribution of Diet Diversity Score",xlab = "DDS")
abline(v = mean(data_final$DDS), col = "red", lwd = 2, lty = 2)
# Summary Statistics
summary(data_final)

## Data Splitting
set.seed(4893)
# 80-20 split
train_id <- sample(1:nrow(data_final),size = 0.8*nrow(data_final))
train <- data_final[train_id,]
valid <- data_final[-train_id,]

## Linear Regression
# Full linear model removing SEQN and unscaled variables
linear_mod <- lm(DDS~ .-SEQN-RIDAGEYR-INDFMPIR-DMDHHSIZ, data = train)
summary(linear_mod)
# VIF to assess multicollinearity
car::vif(linear_mod)
# Diagnostic plot for assumptions
plot(linear_mod)
# Outlier Test
outlierTest(linear_mod)
# RMSE
lm_pred <- round(predict(linear_mod, newdata = valid)*2)/2
sqrt(mean((lm_pred - valid$DDS)^2))

## Lasso Regression
# Model matrix
train_x <- model.matrix(DDS~ .-SEQN-RIDAGEYR-INDFMPIR-DMDHHSIZ, data = train)[, -1]
train_y <- train$DDS
valid_x <- model.matrix(DDS~ .-SEQN-RIDAGEYR-INDFMPIR-DMDHHSIZ,data =valid)[, -1]
set.seed(4893)
# Cross validation for lambda penalty
lasso_cv <- cv.glmnet(train_x, train_y,alpha = 1)
# Refit with penalty chosen
lasso_mod <- glmnet(train_x, train_y,alpha = 1,lambda = lasso_cv$lambda.min)
lasso_mod
coef(lasso_mod)
# RMSE
lasso_pred <- round(predict(lasso_mod, s = lasso_cv$lambda.min, newx = valid_x)*2)/2
sqrt(mean((lasso_pred - valid$DDS)^2))

## Random Forest
set.seed(4893)
# Random Forest
rf_mod <-randomForest(DDS~ .-SEQN-RIDAGEYR-INDFMPIR-DMDHHSIZ, data = train,mtry = 3,importance = T)
rf_mod
# RMSE
rf_pred <- predict(rf_mod,valid)
pred_round <- round(rf_pred * 2) / 2
sqrt(mean((rf_pred - valid$DDS)^2))
# Variable Importance
importance(rf_mod)
varImpPlot(rf_mod)

## Graphs
# Age vs DDS
data_final %>% group_by(RIDAGEYR) %>%
  summarize(score = mean(DDS)) %>%
  ggplot(aes(x = RIDAGEYR,y = score))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x ="Age", y = "Diet Diversity Score", title = "Average Diet Diveristy Score Per Age Group" )
# Food Security vs DDS
data_final %>% ggplot(aes(x = FSDAD,y = DDS))+
  geom_boxplot()+
  labs(x ="Level of Food Security", y = "Diet Diversity Score", title = "Food Security" )
# Income vs DDS
data_final %>%
  mutate(
    income_group = cut(
      INDFMPIR,
      breaks = c(0, 1, 2, 3, 5),
      labels = c("Below poverty",
                 "Near poverty",
                 "Middle income",
                 "High income"),
      include.lowest = TRUE
    )
  )%>%
  ggplot(aes(x = income_group, y = DDS)) +
  geom_boxplot() +
  labs(x = "Income Groups",
       y = "Diet Diversity Score", title = "Income")
# Education vs DDS
data_final %>% ggplot(aes(x = DMDEDUC2,y = DDS))+
  geom_boxplot()+
  labs(x ="Level of Education", y = "Diet Diversity Score", title = "Education" )
