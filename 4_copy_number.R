#copy number estimation and statistical test
source("0_setup.R")

otu = read.delim("resample_zotutab_new_21567.txt",row.names = 1, header = T,check.names = TRUE, sep = "\t")
otu = otu[-which(rowSums(otu)== 0),]
otu = as.data.frame(t(otu))

copynum = read.table("datasets/copy number/copy_number_output.txt",row.names = 1, header = T,check.names = FALSE)
#copynum$copy_number[is.na(copynum$copy_number)] <- 0
#which(is.na(copynum$copy_number))
copynum = copynum[match(colnames(otu),row.names(copynum)),]
#copynum[is.na(copynum)] = 0
#copynum = copynum

#install.packages("matrixStats")
library("matrixStats") 

xm = c()
for (i in 1:96){
  a = sum(as.numeric(otu[i,])/copynum)
  a = 21567/a
  xm = c(xm,a)
}

write.csv(xm,"datasets/copy number/copy_number_per_sample_right.csv")

#linear mixed model
copynum = read.csv("datasets/copy number/copy_num_R.csv",header = T,check.names = FALSE)
copynum$Month = as.factor(copynum$Month)
copynum$Warming = as.factor(copynum$Warming)
copynum$block = as.factor(copynum$block)
copynum$Season = as.factor(copynum$Season)
copynum$Period = as.factor(copynum$Period)
copynum$Period_trt = as.factor(copynum$Period_trt)
copynum$Growth = as.factor(Group$Growing_season)
copynum$G_trt = as.factor(paste(Group$Growing_season,copynum$Warming,sep = "_"))

df <- copynum %>% 
  group_by(Period_trt) %>% 
  dplyr::summarise(avg = mean(copy_number),
                   se = std.error(copy_number))

df <- copynum %>% 
  group_by(Month_trt) %>% 
  dplyr::summarise(avg = mean(copy_number),
                   se = std.error(copy_number))


df <- copynum %>% 
  group_by(G_trt) %>% 
  dplyr::summarise(avg = mean(copy_number),
                   se = std.error(copy_number))

copy_number = copynum$copy_number
Period = copynum$Period
Warming = copynum$Warming
Month = copynum$Month
Growth = copynum$Growth


library(nlme)
library(lme4)

fit = anova(lm(copy_number ~ Month*Warming + block,copynum))
fit
fit = anova(lm(copy_number ~ Period*Warming + block:Month,copynum))
fit

library(lmerTest)
lmer <- lmer(copy_number ~ Month*Warming + (1|block), data = copynum)
summary(lmer)
anova(lmer)

lmer <- lmer(copy_number ~ Period*Warming + (1|block), data = copynum)
summary(lmer)
anova(lmer)

lmer <- lme(fixed = copy_number ~ Month*Warming, data = copynum, random = ~1|block)
summary(aov(lmer))

lmer <- lme(fixed = copy_number ~ Period*Warming, data = copynum, random = ~ block)
summary(aov(lmer))

lmer <- lme(copy_number ~ Period:Warming+Period+Warming+block, data = copynum, 
            fixed = Period:Warming +Period+Warming,
            random = block ~ 1 )
summary(aov(lmer))

#linear mixed model -tNST
tNST = read.csv("datasets/NST/NST_R.csv",header = T,check.names = FALSE)
tNST$Period = as.factor(tNST$Period)
tNST$NST = as.numeric(tNST$NST.ij.bray)
tNST$Warming = as.character(tNST$Warming)
tNST$Month = as.character(tNST$Month)

#copynum$Season = as.factor(copynum$Season)

library(lmerTest)

fit = anova(lm(NST.ij.bray ~ Warming*Period, data = tNST))
fit
fit = anova(lm(tNST$NST ~ tNST$Warming*tNST$Month))
fit

#linear mixed model - beta dispersion
dis = read.csv("datasets/dissimilarity/asv_within_for_anova.csv",header = T,check.names = FALSE)
dis$Period = as.factor(dis$Period)
dis$dis= as.numeric(dis$dis)
dis$Warming = as.character(dis$trt1)
dis$Month = as.character(dis$month1)

#copynum$Season = as.factor(copynum$Season)

library(lmerTest)

fit = anova(lm(dis ~ Warming*Period, data = dis))
fit
fit = anova(lm(dis ~ Warming*Month, data = dis))
fit

