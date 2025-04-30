library(dplyr)
library(pROC)
library(glmnet)
library(doParallel)
library(PRROC)
library(ggplot2)
library(gridExtra)

source("~/projects/eICU_TL/utils.R")


eICU = read.csv("~/projects/eICU_TL/alldata_cleaned.csv")


eICU$age_group = as.factor(eICU$age_group)
eICU$apachedxgroup = as.factor(eICU$apachedxgroup)

eICU = eICU %>%
  filter(apachedxgroup == "Sepsis",
         ethnicity != "Other/Unknown") %>%
  select(-apachedxgroup)

# Run this once to save the glm.source model and then comment it out
# fit.source.model(eICU, ethnicity = "Caucasian", nfold = 5, seed = 12121)
# fit.source.model(eICU, ethnicity = NULL, nfold = 5, seed = 12121)

glm.source = readRDS(file = "~/projects/eICU_TL/glm_source.rds")

all.measures = data.frame()
all.beta.target = data.frame()


minorities = c("Hispanic", "African American", "Asian", "Native American")
nfold = 3
# minority = "Native American"

for(i.seed in c(12121, 67287, 72839))
{
  for(minority in minorities){
    
    set.seed(i.seed)
    eICU.min = eICU %>% filter(ethnicity == minority) %>% select(-ethnicity)
    y = as.numeric(eICU.min$actualhospitalmortality == "EXPIRED")
    fold0 <- sample.int(sum(y==0)) %% nfold
    fold1 <- sample.int(sum(y==1)) %% nfold
    foldid <- numeric(length(y))
    foldid[y==0] <- fold0
    foldid[y==1] <- fold1
    foldid <- foldid + 1
    foldid.min = foldid
    
    for(id in 1:nfold){
      cat("running minority", minority,"id", id, "seed", i.seed, "...\n")
      eICU.min.train = eICU.min[foldid != id, ]
      eICU.min.test = eICU.min[foldid == id, ]
      
      measure = "auc"
      
      ### Only target
      m = only.target(eICU.min.train, eICU.min.test, nfold = 5, verbose = F, measure = measure)
      all.measures = rbind(all.measures, c(i.seed, id, minority, "only.target", m[[1]]))
      all.beta.target = rbind(all.beta.target, m[[2]])
      
      ## two step only main effects
      m = only.main.effect(eICU.min.train, eICU.min.test, glm.source, nfold = 5, verbose = F, measure = measure)
      all.measures = rbind(all.measures, c(i.seed, id, minority, "only.main", m[[1]]))
      all.beta.target = rbind(all.beta.target, m[[2]])
      
      ## two step all effects
      m = all.effect(eICU.min.train, eICU.min.test, glm.source, nfold = 5, verbose = F, measure = measure)
      all.measures = rbind(all.measures, c(i.seed, id, minority, "all.effects", m[[1]]))
      all.beta.target = rbind(all.beta.target, m[[2]])
      
      ## Adaptive
      m = adaptive(eICU.min.train, eICU.min.test, glm.source, nfold = 5, verbose = F, measure = measure)
      all.measures = rbind(all.measures, c(i.seed, id, minority, "adaptive", m[[1]]))
      all.beta.target = rbind(all.beta.target, m[[2]])
      
      m = adaptive.main(eICU.min.train, eICU.min.test, glm.source, nfold = 5, verbose = F, measure = measure)
      all.measures = rbind(all.measures, c(i.seed, id, minority, "adaptive.main", m[[1]]))
      all.beta.target = rbind(all.beta.target, m[[2]])
    }
    
  }
}


colnames(all.measures) <- c("seed", "foldid", "target", "method", "auc",
                            "auprc", "sensitivity", "specificity", "gamma")
write.csv(all.measures, file = "~/projects/eICU_TL/all_measures_auc.csv", row.names = FALSE)
write.csv(all.beta.target, file = "~/projects/eICU_TL/all_coef.csv", row.names = FALSE)


all.measures = read.csv(file = "~/projects/eICU_TL/all_measures_auc.csv")

all.measures$auc = as.numeric(all.measures$auc)
all.measures$auprc = as.numeric(all.measures$auprc)
all.measures$sensitivity = as.numeric(all.measures$sensitivity)
all.measures$specificity = as.numeric(all.measures$specificity)
all.measures$balanced_acc = 0.5 * (all.measures$sensitivity + all.measures$specificity)


p_auc = ggplot(all.measures, aes(x = target, y = auc, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8))
p_auprc = ggplot(all.measures, aes(x = target, y = auprc, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8))
p_TPR = ggplot(all.measures, aes(x = target, y = sensitivity, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8))
p_TNR = ggplot(all.measures, aes(x = target, y = specificity, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8))
p_bal_acc = ggplot(all.measures, aes(x = target, y = balanced_acc, fill=method)) +
  geom_boxplot(position = position_dodge(width = 0.8))


grid.arrange(p_auc, p_auprc, p_TPR, p_TNR, p_bal_acc, ncol = 1)

# dist.df = data.frame()
# 
# for(minority in minorities){
#   for(method in unique(all.measures$method)){
#     idx = as.numeric(all.measures$target == minority) * as.numeric(all.measures$method == method)
#     idx = which(idx == 1)
#     coeffs = all.beta.target[idx, ]
#     distances = as.vector(dist(coeffs))
#     df = data.frame(pairwise.dist = distances)
#     df$target = minority
#     df$method = method
#     dist.df = rbind(dist.df, df)
#   }
# }
# 
# ggplot(dist.df, aes(x = target, y = pairwise.dist, fill=method)) +
#   geom_boxplot(position = position_dodge(width = 0.8))
