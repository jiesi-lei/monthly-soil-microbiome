# SEM for both Rh and Rs
library(dplyr)
library(lavaan)
library(AICcmodavg)
library(tidyverse)
library(semPlot)
library(vegan)

source("0_Setup.R")
source("0_Setup_gene.R")

# source("lavaan.modavg.R")
#group
Group<-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 
Group = Group[match(rownames(gene),row.names(Group)),]
Group$Season = ordered(Group$Season, c("sp","sm","fl","wt"))
Group$Growing_season = ordered(Group$Growing_season, c("G","NG"))
Group$Period = ordered(Group$Period, c("Cold season 1","Warm season","Cold season 2"))
head(Group)
id_cool = which(Group$Period %in% c("Cold season 1","Cold season 2"))
id_warm = which(Group$Period == "Warm season" )

#preparation-----------------------------------
#functional profile
dat = read.csv("cut.treat.2.log.RA.csv",row.names = 1)
gene = as.data.frame(t(dat[,-c(1:9)]))
which((colSums(gene)) == 0)
info = dat[,c(1:9)]
set.seed(2000)
nmds1 <- metaMDS(gene, distance = 'bray', k = 4)
nmds1
summary(nmds1)
nmds_plot <- nmds1
func.struc = as.data.frame(scores(nmds1))[,1]*10
func.richness = specnumber(gene)/20000

#taxonomic profile: shannon index
otu = read.delim("resample_zotutab_new_21567.txt",row.names = 1, header = T,check.names = TRUE, sep = "\t")
otu = otu[-which(rowSums(otu)== 0),]
otu = as.data.frame(t(otu))
bactr.richness = specnumber(otu)/1000
nmds1 <- metaMDS(log(otu+1), distance = 'bray', k = 4)
nmds1
summary(nmds1)
nmds_plot <- nmds1
bactr.struc = as.data.frame(scores(nmds1))[,1]*30

# c.deg gene absolute abundance
c.deg.select = c("amyA","cda","glucoamylase","nplT","pula",#Starch
                 "ara","mannanase","xyla","xylanase",#Hemicellulose
                 "cellobiase","endoglucanase","exoglucanase",#Cellulous
                 "chitinase",#Chitin
                 "pectinase (pectate_lyase)",#Pectin
                 "glx","mnp","phenol_oxidase",#lignin
                 "cutinase",#cutin
                 "cdh")#Terpene
subselect = which(info[,"gene"] %in% c.deg.select)
subgene <- get(paste0("gene"))[,subselect]
C_deg_gene_abun = rowSums(subgene)

#copy number
copynum = read.csv("datasets/copy number/copy_num_R.csv",header = T,check.names = FALSE)
copy_number = copynum$copy_number

#Env
env <- read.csv("monthly_ENV.csv",row.names = 1)
env = env[match(rownames(gene),row.names(env)),]

warming = as.numeric(rep(rep(c("1","2"),each = 4), times = 12))
moisture = env$GWC
pH = env$pH
TC = env$Tcpercent*100
TN = env$Tnpercent*100
NO3N = env$NO3N
NH4N = env$NH4N
temp = env$X3dayTemperatureAt7.5cm
precipitation = env$Rain/10
doy = env$DayOfSample

RA = env$RootResp
RH = env$HeteroResp
NEE = env$NEE
GPP = env$GPP

RS = env$SoilResp

semdat = data.frame(moisture, pH, temp, TN, TC, NH4N, NO3N, precipitation, 
                    RA, NEE, GPP, RS, RH,
                    func.struc, func.richness,
                    bactr.struc, bactr.richness,
                    copy_number, doy, warming)
#write.table(semdat,"datasets/SEM/semdata.dat")
#semdat %>% {scale(.)} %>% {data.frame(.,warming)}-> semdat.scaled

semdat.scaled.cool = scale(semdat[id_cool,])
semdat.scaled.warm = scale(semdat[id_warm,])
save(semdat.scaled.cool,semdat.scaled.warm,file = "SEM_RH/SEM_dat.RData")

library(psych)
semdat.scaled.cool1 = semdat.scaled.cool[-8,]
cor(semdat.scaled.cool1)
pairs.panels(semdat.scaled.cool)
pairs.panels(semdat.scaled.warm)


# hypothesized conceptual model
mod.8a <- 'temp ~ warming
           moisture ~ warming
           NO3N ~ warming
           GPP ~ warming + pH
           pH ~ warming + NO3N
           func.struc ~ moisture + NO3N + GPP + temp
           bactr.richness ~ moisture + NO3N + GPP + temp + pH
           copy_number ~ moisture + temp + NO3N + func.struc + bactr.richness + pHv
           RH ~ moisture + temp + NO3N + func.struc + bactr.richness + copy_number
           RS ~ moisture + temp + NO3N + func.struc + bactr.richness + copy_number + GPP
'

# final model
# Cool season ----
mod.8a <- 'temp ~ warming
           moisture ~ warming
           NO3N ~ warming
           GPP ~ moisture + pH + warming
           pH ~ warming + NO3N
           func.struc ~ temp + moisture + NO3N + GPP
           bactr.richness ~ pH + temp + moisture
           copy_number ~ moisture + temp + NO3N + func.struc
           RH ~ bactr.richness + temp + pH + copy_number 
           RS ~ func.struc + GPP + moisture + NO3N + temp
           bactr.richness ~~ RS
           func.struc ~~ bactr.richness
'

fit.8a <- sem(mod.8a, data=semdat.scaled.cool1, estimator = "GLS")#bactr.richness + NEE + temp + precipitation + moisture
semPaths(fit.8a, "std.all", edge.label.cex = 0.5, exoVar = FALSE, 
         exoCov = FALSE)
summary(fit.8a, rsq=T, standardize = T)

fitMeasures(fit.8a)

resid(fit.8a, type="standardized")
modindices(fit.8a)
mi1<-modindices(fit.8a);print(mi1[mi1$mi>1.0,])

Beta <- lavInspect(fit.8a, "std.all")$beta # matrix of path coefficients
# view effects of each order
Beta # direct effects
Beta %*% Beta # indirect effects involving 2 path coefficents
Beta %*% Beta %*% Beta # indirect effects involving 3 path coefficients
Beta %*% Beta %*% Beta %*% Beta # indirect effects involving 4 path coefficients
Beta %*% Beta %*% Beta %*% Beta %*% Beta
Beta %*% Beta %*% Beta %*% Beta %*% Beta %*% Beta

toteff = as.data.frame(Beta + Beta %*% Beta + Beta %*% Beta %*% Beta +
                         Beta %*% Beta %*% Beta %*% Beta + 
                         Beta %*% Beta %*% Beta %*% Beta %*% Beta)
RS_eff = toteff["RS",]
RH_eff = toteff["RH",]
write.csv(toteff, "SEM_RH/SEM_total_effects_cool_season.csv", row.names = T)

#Warm season ----
# initial
mod.8a <- 'temp ~ warming
           moisture ~ warming
           NO3N ~ warming
           GPP ~ warming + pH
           pH ~ warming + NO3N
           func.struc ~ moisture + NO3N + temp
           bactr.richness ~ moisture + NO3N + temp + pH
           copy_number ~ func.struc + temp
           RH ~ moisture + func.struc
           RS ~ moisture + temp + func.struc + GPP
           NO3N ~~ bactr.richness
'

fit.8a <- sem(mod.8a, data=semdat.scaled.warm, estimator = "GLS")
semPaths(fit.8a)
summary(fit.8a, rsq=T, standardize = T)
fitMeasures(fit.8a)

resid(fit.8a, type="standardized")
modindices(fit.8a)
mi1<-modindices(fit.8a);print(mi1[mi1$mi>2.0,])

Beta <- lavInspect(fit.8a, "std.all")$beta # matrix of path coefficients
# view effects of each order
Beta # direct effects
Beta %*% Beta # indirect effects involving 2 path coefficents
Beta %*% Beta %*% Beta # indirect effects involving 3 path coefficients
Beta %*% Beta %*% Beta %*% Beta # indirect effects involving 4 path coefficients
Beta %*% Beta %*% Beta %*% Beta %*% Beta
Beta %*% Beta %*% Beta %*% Beta %*% Beta %*% Beta

toteff = as.data.frame(Beta + Beta %*% Beta + Beta %*% Beta %*% Beta +
                         Beta %*% Beta %*% Beta %*% Beta + 
                         Beta %*% Beta %*% Beta %*% Beta %*% Beta)
RS_eff = toteff["RS",]
RH_eff = toteff["RH",]
write.csv(toteff, "SEM_RH/SEM_total_effects_warm_season.csv", row.names = T)

save.image("SEM_RH/SEM_revised_RHRS.RData")
load("SEM_RH/SEM_revised_RHRS.RData")
