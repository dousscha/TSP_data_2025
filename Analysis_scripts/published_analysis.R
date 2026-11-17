# ============================================================================
#                           MODELS AND COMPARISONS
# ============================================================================

library(lme4)
library(MuMIn)
library(performance)
library(emmeans)
setwd("~/Shahed_Traplines/Second_experiment_2024/Analysis")
source(file = "glmm_funs.R")
sessionInfo()

# ============================================================================
# 1.  RELATIVE ERROR in EXPERIMENT 2 (TABLE 2; SUPPLEMENT S6)
# ============================================================================

set.seed(123)
GLMM_df_paper <- read.csv("~/Shahed_Traplines/Second_experiment_2024/Combined_data/GLMM_df_paper.csv", header=TRUE)
data <- data.frame(GLMM_df_paper)
data <- na.omit(data)
dim(data)

data$zscore <- (abs(data$error_traj-mean(data$error_traj))/sd(data$error_traj))
data <- subset(data, data$zscore < 3)
dim(data)

data$cond <- as.factor(data$cond)

data$scale <- as.factor(data$scale)
data$scale<- relevel(data$scale, ref = "S3")
data$num_scale <-gsub("S","",as.character(data$scale))
data$num_scale <-as.numeric(data$num_scale)
data$scale_ordi <- ordered(data$num_scale, levels = 1:5, labels = c("XS", "S", "M", "L","XL"))

# allow gamma fit
data[data==0] <- 0.01
data$response <-as.numeric(data$error_traj)
data$session <- as.numeric(data$session_num)
data$scale   <- as.factor(data$scale)
data$Gender <- as.factor(data$Gender)

data$gaming_habits <-as.numeric(data$gaming_habits)
data$trial_num <-as.numeric(data$trial_num)
data$GPS_uage <-as.numeric(data$GPS_uage)
data$Gender <-as.factor(data$Gender)
data$Gender<- relevel(data$Gender, ref = "male")

data_clean <- na.omit(data)  # Remove rows with NA values
data_clean <- data_clean[!apply(data_clean == Inf | data_clean == -Inf, 1, any), ]  # Remove rows with Inf
data <-data_clean



model_base       <- glm(response ~ num_scale *cond + I(num_scale^2) *cond+ trial_num + gaming_habits +
                          GPS_uage + daily_optimisation + session_num + Age + Gender
                        + navigation_rating ,data = data,
                        family=Gamma(link = "log"))
summary(model_base)

model_no_cond <- glm(response ~ num_scale  + I(num_scale^2) + trial_num + gaming_habits +
                       GPS_uage + daily_optimisation + session_num + Age + Gender
                     + navigation_rating ,data = data,
                     family=Gamma(link = "log"))

summary(model_no_cond)

compare_performance(model_no_cond,model_base)

model_no_cond_ppid <- glmer(response ~ num_scale  + I(num_scale^2) + trial_num + gaming_habits +
                       GPS_uage + daily_optimisation+ Age +(1|ppid)
                      ,data = data,
                     family=Gamma(link = "log"))

summary(model_no_cond_ppid)

# manuscript Table 2:
# Gamma(log): response ~ num_scale + I(num_scale^2) + trial_num
#                      + gaming_habits + (1 | ppid).
model_simplified1 <- glmer(response ~ num_scale + I(num_scale^2) + trial_num + gaming_habits
                         +(1|ppid),data = data,
                         family=Gamma(link = "log"),control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

summary(model_simplified1)

# Model comparisons for the reported predictor selection.
aic_tab1 <- MuMIn::model.sel(model_simplified1,model_no_cond_ppid)
aic_table1 <- round(aic_tab1[, c("AICc", "delta", "df")], digits = 2)
aic_table1
r2(model_simplified1)
summary(model_simplified1)


# Additional Gamma candidate models 
model_ppid_simplified2 <- glmer(response ~ num_scale + I(num_scale^2) + trial_num + gaming_habits +
                      GPS_uage +
                       (1|ppid), data = data,
                    family = Gamma(link = "log"),
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

model_ppid_simplified1 <- glmer(response ~ num_scale + I(num_scale^2) + trial_num + gaming_habits +Age+ daily_optimisation+
                     (1 | ppid), data = data,
                    family = Gamma(link = "log"),
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

aic_tab1 <- MuMIn::model.sel(model_no_cond_ppid, model_ppid_simplified1,
                               model_ppid_simplified2, model_simplified1)
aic_table1 <- round(aic_tab1[, c("AICc", "delta", "df")], digits = 2)
aic_table1


# ============================================================================
# 2.  GLMM model fitting on TSP optimal success TABLE 3 and S7 in EXPERIMENT 2  
# ============================================================================

data <- data.frame(GLMM_df_paper)
data <- na.omit(data)

dim(data)
data$cond <- as.factor(data$cond)
data$scale <- as.factor(data$scale)
data$Gender<- as.factor(data$Gender)
data$Gender<- relevel(data$Gender, ref = "male")

data$num_scale <-gsub("S","",as.character(data$scale))
data$num_scale <-as.numeric(data$num_scale)
data$scale<- relevel(data$scale, ref = "S1")

# we expect the relation to be non linear model
model1<- glm(TSP_correct ~ num_scale *cond + I(num_scale^2) *cond
                + trial_num + gaming_habits + GPS_uage + daily_optimisation + session_num + Age  + Gender+
                + navigation_rating ,data = data, family = binomial(link = "logit"))

summary(model1)

model1bis <- glm(TSP_correct ~ num_scale *cond+ I(num_scale^2) *cond
              + gaming_habits +trial_num ,data = data,family = binomial(link = "logit"))
r.squaredGLMM(model1bis)

aic_tab1  <- MuMIn::model.sel(model1,model1bis)
aic_table1 <- round(aic_tab1[, c("AICc", "delta", "df")], digits = 2)
aic_table1
anova(model1,model1bis,test='Chisq')

# MODEL WITH SCALE × CONDITION INTERACTIONS (Supplement S7 model 2).
model2 <- glmer(TSP_correct ~ num_scale*cond + I(num_scale^2)*cond  + gaming_habits + trial_num+
                              (1  | ppid),
                            family = binomial(link = "logit"), data = data,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

table(data$TSP_correct, data$num_scale)

r2(model2)
r2(model1)

# SELECTED TSP MODEL:
# retains condition main effect, removes scale × condition interactions.
model3 <- glmer(TSP_correct ~ num_scale + I(num_scale^2) +cond + gaming_habits + trial_num+
                              (1  | ppid),
                            family = binomial(link = "logit"), data = data)

# comparison model2 VS model3.
summary(model2)
summary(model3)
aic_tab1 <- MuMIn::model.sel(model1, model2, model3)
aic_table1 <- round(aic_tab1[, c("AICc", "delta", "df")], digits = 2)
aic_table1
r2(model2)
r2(model3)

# ============================================================================
# 3. Fig1 & S5: GLMM model fitting on Relative Error EXP 1: 
# ============================================================================

set.seed(123)
GLMM_df <- read.csv("~/Shahed_Traplines/Second_experiment_2024/Combined_data/GLMM_Exp1.csv", header=TRUE)
data <- data.frame(GLMM_df)
data <- na.omit(data)

dim(data)
data$Scale <- as.factor(data$Scale)
data$Scale<- relevel(data$Scale, ref = "S2")

data$num_scale <-gsub("S","",as.character(data$Scale))
data$num_scale <-as.numeric(data$num_scale)
data$num_scale <- data$num_scale - 2

data$zscore <- (abs(data$perc_error-mean(data$perc_error))/sd(data$perc_error))
data <- subset(data, data$zscore < 3)

data[data==0]<- 0.001
data$response <-as.numeric(data$perc_error)
data$session <- as.numeric(data$session_num)
data$Scale   <- as.factor(data$Scale)
data$gender   <- as.factor(data$gender)
data$Configuration   <- as.factor(data$Configuration)
data$Configuration<- relevel(data$Configuration, ref = "C1")

data$gaming_habits <-as.numeric(data$gaming_habits)
data$trial_num <-as.numeric(data$trial)
data$GPS_usage <-as.numeric(data$GPS_usage)
data$gender <-as.factor(data$gender)
data$gender<- relevel(data$gender, ref = "male")

model_base  <- glm(response ~ num_scale *Configuration + I(num_scale^2) *Configuration ,
                   data = data,family=Gamma(link = "log"))
summary(model_base)
r2(model_base)

slopes <- emtrends(model_base, ~ Configuration, var = "num_scale")
slopes
plot(slopes, comparisons = TRUE)
pairs(slopes, adjust = "tukey")  # Tukey correction for multiple comparisons

model_base2  <- glm(response ~ num_scale *Configuration + I(num_scale^2) *Configuration 
                    + trial + session_num, data = data,family=Gamma(link = "log"))
summary(model_base2)
r2(model_base2)

model_base_ppid <- glmer(response ~ num_scale *Configuration + I(num_scale^2) *Configuration
                          +(1|user),data = data,
                        family=Gamma(link = "log"),control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e9)))
summary(model_base_ppid)

r2(model_base_ppid)

model_base2_ppid <- glmer(response ~ num_scale *Configuration
                         + session_num + trial+(1|user),data = data,
                         family=Gamma(link = "log"),control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

summary(model_base2_ppid)

r2(model_base2_ppid)


# ============================================================================
# 4. EXPERIMENT 2 — CHEST VISIBILITY
# ============================================================================

set.seed(123)
df_chests <- read.csv("~/Shahed_Traplines/Second_experiment_2024/Combined_data/GLMM_chest_visibility.csv", header=TRUE)
data <- data.frame(df_chests)
data <- na.omit(df_chests)
dim(data)
data$num_scale <-gsub("S","",as.character(data$scale))
data$num_scale <-as.numeric(data$num_scale)

model_base_poly <- glm(total_in_sight ~ poly(num_scale,2) *cond , data = data)
summary(model_base_poly)

model_base <- glm(total_in_sight ~num_scale *cond , data = data)
summary(model_base)

