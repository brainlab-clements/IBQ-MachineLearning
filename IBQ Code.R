# Load Libraries
library(tidyverse)
library(mice)
library(glmnet)
library(caret)
library(doParallel)
library(purrr)
library(patchwork)
library(ggplot2)

# Create functions for bootstrapping
boot_lasso <- function(data, threshold=5) {
  
  model_data <- as.matrix(data)
  
  cv_model <- cv.glmnet(
    x = model_data[,c(1:91)], 
    y = model_data[,92],
    family = "binomial",
    nfolds = 5
  )
  
  lambda <- cv_model$lambda
  no_lambda <- length(lambda)
  
  for (l in seq_len(no_lambda)) {
    current_lambda <- lambda[l]
    model <- glmnet(x=model_data[,c(1:91)], y=model_data[,92], alpha = 1, family = "binomial",
                    lambda = current_lambda)
    item.names <- rownames(model$beta)
    lasso.coef <- data.frame("item"=item.names, "beta"=round(model$beta[,1],3), "odds_rat"=round(exp(model$beta[,1]),3))
    
    if (sum(lasso.coef$beta == 0) < 91-threshold){
      break
    }
  }
  
  return(lasso.coef)
}

boot_ASD <- function(data, colname){
  autism <- data[,colname]
  ASD <- which(autism==1)
  TD <- which(autism==0)
  
  n_ASD <- length(ASD)
  n_TD <- length(TD)
  
  ASD_boot_index <- sample(1:n_ASD, n_ASD, replace = TRUE)
  ASD_boot <- data[ASD_boot_index,]
  
  TD_boot_index <- sample(1:n_TD, n_TD, replace = TRUE)
  TD_boot <- data[TD_boot_index,]
  
  boot_sample <- rbind(ASD_boot, TD_boot)
  return(boot_sample)
}

# Load and clean data

## Load main dataset
data <- read.csv("ENTER DATA")

## Load additional EMO data
emo_asd <- read.csv("ENTER DATA")

## Load all demographics
emo_asd_dem <- read.csv("ENTER DATA")
bsp_demo <- read.csv("ENTER DATA")
emo_demo <- read.csv("ENTER DATA")
isp_demo <- read.csv("ENTER DATA")

## Combine demographics and IBQ data
bsp_ids <- bsp_demo$studyid
emo_ids <- emo_demo$studyid
isp_ids <- isp_demo$studyid

bsp_dem_reduced <- data.frame("studyid"=bsp_demo$studyid, "sex"=bsp_demo$sex, "income"=bsp_demo$income, "race"=bsp_demo$childrace)
emo_dem_reduced1 <- data.frame("studyid"=emo_demo$studyid, "sex"=emo_demo$Sex, "income"=emo_demo$Family_Income, "race"=emo_demo$Race)
emo_dem_reduced1 <- as.data.frame(apply(emo_dem_reduced1, 2, as.character))
emo_dem_reduced2 <- data.frame("studyid"=emo_asd_dem$subj, "sex"=emo_asd_dem$Sex, "income"=emo_asd_dem$Family_Income, "race"=emo_asd_dem$Race)
emo_dem_reduced2 <- as.data.frame(apply(emo_dem_reduced2, 2, as.character))
emo_dem_all <- bind_rows(emo_dem_reduced1, emo_dem_reduced2)

isp_dem_reduced <- data.frame("studyid"=isp_demo$studyid, "sex"=isp_demo$dem_participant_sex, "income"=isp_demo$dem_12_month_income, "race"=isp_demo$dem_child_race)

reduced_data <- data[,c(1:94, 96)]
combined_bsp <- reduced_data[reduced_data$study=="bsp",] %>% merge(bsp_dem_reduced, by="studyid", all.y=T)


ibq_emo_names <- colnames(emo_asd)[4:94]
ibq_data_names <- colnames(data)[4:94]

## Reverse code 13 EMO participants
colnames(emo_asd)[4:94] <- ibq_data_names
colnames(emo_asd)[2] <- "age_group"
emo_asd$ibq33 <- 8 - emo_asd$ibq33
emo_asd$ibq3 <- 8 - emo_asd$ibq3
emo_asd$ibq82 <- 8 - emo_asd$ibq82
emo_asd$ibq84 <- 8 - emo_asd$ibq84
emo_asd$ibq86 <- 8 - emo_asd$ibq86
emo_asd$ibq91 <- 8 - emo_asd$ibq91
emo_asd$ibq37 <- 8 - emo_asd$ibq37
emo_asd$ibq72 <- 8 - emo_asd$ibq72
emo_asd$ibq61 <- 8 - emo_asd$ibq61
emo_asd$ibq50 <- 8 - emo_asd$ibq50
emo_asd$ibq60 <- 8 - emo_asd$ibq60
emo_asd$ibq73 <- 8 - emo_asd$ibq73

## Combine all data
reduced_emo_asd <- emo_asd[,c(1:94, 126)]
asd_td_emo <- bind_rows(reduced_data[reduced_data$study=="emo",], reduced_emo_asd)
combined_emo <- asd_td_emo %>% merge(emo_dem_all, by="studyid", all.y=T)

combined_isp <- reduced_data[reduced_data$study=="isp",] %>% merge(isp_dem_reduced, by="studyid", all.y=T)

all_data <- bind_rows(combined_bsp, combined_emo, combined_isp)

## Recode sex variable
all_data$sex <- all_data$sex %>% recode("0"="Male",
                                        "1"="Female",
                                        "female"="Female",
                                        "male"="Male")

## Exclude developmental concern participants
TD <- all_data[all_data$dx==0 & !is.na(all_data$dx),]
ASD <- all_data[all_data$dx==1 & !is.na(all_data$dx),]

### isolate TD and ASD
TD_ASD <- rbind(TD, ASD)

### get rid of duplicates
TD_ASD_nodup <- NULL
all_ids <- unique(TD_ASD$studyid)
for (i in 1:length(all_ids)) {
  current <- all_ids[i]
  case <- TD_ASD[TD_ASD$studyid==current,]
  len <- nrow(case)
  
  if (len > 1){
    case <- case[nrow(case),]
  }
  
  TD_ASD_nodup <- rbind(TD_ASD_nodup, case)
}

ibq <- TD_ASD_nodup %>% select(starts_with("ibq"))
ibq <- cbind("studyid"=TD_ASD_nodup$studyid, ibq)

ibq_ids <- ibq$studyid
ibq_dem <- TD_ASD_nodup[TD_ASD_nodup$studyid %in% ibq_ids,c(1:3,95:98)]

## Impute IBQ items using age and gender as auxiliary variables
df <- ibq[,-1]
df$sex <- ibq_dem$sex
df$age_group <- ibq_dem$age_group
ini <- mice(df, maxit = 0, printFlag = FALSE)
pred <- ini$predictorMatrix
pred[, "age_group"] <- 1
pred[, "sex"] <- 1
pred_optimized <- quickpred(
  df, 
  mincor = 0.1, 
  include = c("age_group", "sex")
)

imputed_data <- mice(
  df,
  m = 5,                         # Number of multiply imputed datasets
  maxit = 10,                    # Number of iterations per dataset
  predictorMatrix = pred_optimized, # Use the predictor matrix forcing age & sex
  seed = 123,                    # Set seed for reproducibility
  printFlag = TRUE
)

df_completed <- complete(imputed_data, action = 1)
df_completed <- df_completed[,c(-92:-93)]

imputed <- cbind("studyid"=ibq_ids, df_completed)
imputed <- imputed %>% merge(ibq_dem, by="studyid")

# Fit LASSO model to data

## set reproducible seed
set.seed(250624)

## create 5-fold cross validated samples stratified by age, gender, and autism
imputed$strata <- factor(paste(imputed$dx, imputed$age_group, imputed$sex, sep = "_"))
folds <- createFolds(imputed$strata, k = 5, list = TRUE, returnTrain = FALSE)
fold1_test_indices  <- imputed$studyid[folds[[1]]]
fold1_train_indices <- imputed$studyid[-folds[[1]]]
imputed$sex_num <- ifelse(imputed$sex == "Male", 0, 1)
test_fold1 <- imputed[folds[[1]], ]
imputed$sex_num <- ifelse(imputed$sex == "Male", 0, 1)
imputed <- imputed %>%
  group_by(dx, age_group, sex) %>%
  mutate(fold = sample(rep(1:5, length.out = n()))) %>%
  ungroup()

# look at balance
# by(imputed, imputed$fold, function(sub) {
#   c(
#     n = nrow(sub),
#     mean_age = mean(sub$age_group),
#     pct_female = mean(sub$sex == "Female"),
#     pct_dx1 = mean(sub$dx == "1")
#   )
# })

imputed <- as.matrix(imputed)

### isolate variables used in analysis
model_data <- apply(imputed[,c(2:92,95,101)], 2, as.numeric)

cv_model <- cv.glmnet(
  x = model_data[,c(1:91)], 
  y = model_data[,92],
  family = "binomial",
  foldid = model_data[,93],
  nfolds = 5
)

plot(cv_model)

lambda <- cv_model$lambda
no_lambda <- length(lambda)

### search for model that satisfies item threshold (1)
for (l in seq_len(no_lambda)) {
  current_lambda <- lambda[l]
  model <- glmnet(x=model_data[,c(1:91)], y=model_data[,92], alpha = 1, family = "binomial",
                  lambda = current_lambda)
  item.names <- rownames(model$beta)
  lasso.coef <- data.frame("item"=item.names, "beta"=round(model$beta[,1],3), "odds_rat"=round(exp(model$beta[,1]),3))
  
  ### stop search when condition is satisfied (at least 1 item is non-zero)
  if (sum(lasso.coef$beta == 0) < 90){
    break
  }
}

remaining <- lasso.coef[lasso.coef$beta!=0,]
coef <- remaining[order(remaining$odds_rat, decreasing = T),]

# Conduct bootstrap LASSO

## set number of bootstrap samples (100,000)
B <- 100000
## specify number of cores to use for parallel
cores <- 12

## register cluster and begin parallel
### WARNING BOOTSTRAP TAKES A LONG TIME
cl <- makeCluster(cores)
registerDoParallel(cl)

all_boot <- foreach(b = seq_len(B), .packages = c("glmnet", "caret"), 
                    .export = c("model_data", "boot_ASD", "boot_lasso")) %dopar% {
                      boot_sample <- boot_ASD(data = model_data, colname = 92)
                      boot_lasso(boot_sample, threshold=1)
                    }

## end parallel
stopCluster(cl)

# Obtain bootstrap CI and SE
boot_matrix <- do.call(rbind, lapply(all_boot, function(x) x$beta))

## Compute standard error
boot_se <- sqrt((colSums((boot_matrix - mean(boot_matrix))^2))/(B-1))

## Compute CI (.025 to .975)
ci_bounds <- apply(boot_matrix, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

ci_df <- data.frame(
  Coefficient = if (!is.null(colnames(boot_matrix))) colnames(boot_matrix) else paste0("ibq", 1:ncol(boot_matrix)),
  Lower_CI    = ci_bounds[1, ],
  Upper_CI    = ci_bounds[2, ],
  SE          = boot_se
)

## Store estimates from original model (no bootstrap)
ci_df$Estimate <- lasso.coef$beta

## Store number of times a coefficient is non-zero across bootstraps
counts <- rep(NaN, 91)
for (i in seq_len(91)) {
  var <- boot_matrix[,i]
  counts[i] <- sum(var!=0)
}

ci_df$Counts <- counts
ci_df$Counts_percent <- counts/B

## Store all results in csv
write.csv(ci_df, "lasso_ci.csv", row.names=T)

# Create forest plot
ci_data <- ci_df
ci_sort <- ci_data[order(ci_data$Counts_percent, decreasing=T), ]
ci_data$Counts_percent

ci_10 <- ci_sort[1:10,]

## Preserve order as in the dataset (top item at the top of the graph)
ci_10$X <- factor(ci_10$X, levels = rev(ci_10$X))

## if using valence instead of subscale color scheme
ci_10$Group <- factor(c(names(custom_colors)[1],
                        names(custom_colors)[2],
                        names(custom_colors)[1],
                        names(custom_colors)[2],
                        names(custom_colors)[1],
                        names(custom_colors)[1],
                        names(custom_colors)[1],
                        names(custom_colors)[2],
                        names(custom_colors)[2],
                        names(custom_colors)[1])
)

ci_10$Significant <- with(ci_10, (Lower_CI > 0 & Upper_CI > 0) | (Lower_CI < 0 & Upper_CI < 0))


## Shared Theme Setup
theme_clean <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "#777777", size = 0.5),
    axis.ticks = element_line(color = "#777777"),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  )


## Panel A: Selection Stability (Left Plot)
p1 <- ggplot(ci_10, aes(y = X, x = Counts_percent, color = Group)) +
  # Lollipop stems starting from 0
  geom_segment(aes(x = 0, xend = Counts_percent, yend = X), size = 0.8) +
  # Lollipop points
  geom_point(size = 3.5) +
  # Vertical dashed reference line (e.g., permutation null at 0.18)
  # geom_vline(xintercept = 0.18, linetype = "dashed", color = "gray50") +
  # Text label for dashed line
  # annotate("text", x = 0.22, y = length(ci_10$X) + 0.6,
  #          label = "Permutation null (95th pct)", hjust = 0, size = 3.8, color = "gray30") +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, 0.1),
    labels = function(x) paste0(x * 100, "%"),
    expand = c(0, 0)
  ) +
  labs(x = "Selection Rate") +
  theme_clean +
  theme(legend.position = "none")


## Panel B: Effect Size (Right Plot)
p2 <- ggplot(ci_10, aes(y = X, x = Estimate, color = Group)) +
  # Reference line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # Confidence Interval horizontal lines
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, size = 0.8) +
  # Points: Filled if significant, Open (white background) if crossing zero
  geom_point(aes(fill = ifelse(Significant, as.character(Group), "white")),
             shape = 21, size = 3, stroke = 1.2) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = c(custom_colors, "white" = "white")) +
  scale_x_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.1, 0.1, .1)) +
  labs(x = "Beta Estimate (95% interval)") +
  theme_clean +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),  # Hide Y-axis text for panel B
    axis.ticks.y = element_blank(), # Hide Y-axis ticks for panel B
    axis.line.y = element_blank()   # Remove left axis line on plot B
  )


## Combine and save plots
combined_plot <- p1 + p2 + plot_layout(widths = c(1.2, 1))

ggsave("ibqfigure.png", plot = combined_plot, width = 6, height = 4, units = "in", dpi = 300)















