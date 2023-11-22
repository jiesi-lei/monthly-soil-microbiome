library(lavaan)
library(AICcmodavg)
library(tidyverse)
library(semPlot)
library(vegan)

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
semdat %>% {scale(.)} %>% {data.frame(.,warming)}-> semdat.scaled

semdat.scaled.cool = semdat[id_cool,]
semdat.scaled.warm = semdat[id_warm,]

# write.csv(semdat.scaled.cool,"SEM/cool.csv",row.names = F)
# write.csv(semdat.scaled.warm,"SEM/warm.csv",row.names = F )

semdat.scaled.cool = read.csv("SEM/cool.csv")
semdat.scaled.warm = read.csv("SEM/warm.csv")

#initial model
mod.8a <- 'func.struc ~ moisture + pH + temp + TN + TC + NH4N + NO3N + precipitation + GPP
           func.richness ~ moisture + pH + temp + TN + TC + NH4N + NO3N + precipitation + GPP
           bactr.struc ~ moisture + pH + temp + TN + TC + NH4N + NO3N + precipitation + GPP
           bactr.richness ~ moisture + pH + temp + TN + TC + NH4N + NO3N + precipitation + GPP
           copy_number ~ moisture + pH + temp + TN + TC + NH4N + NO3N + precipitation + GPP
           GPP ~ moisture + pH + temp + TN + TC + NH4N + NO3N + precipitation
           RS ~ func.struc + func.richness  + bactr.struc + bactr.richness + copy_number + GPP'
fit.8a <- sem(mod.8a, data = semdat)
summary(fit.8a, rsq=T, standardize = T)
#The Std.lv column shows results that are standardized so that the latent variables 
#have a variance of one. The Std.all column shows results that are standardized 
#so that both the latent variables and the observed variables have a variance of one.
#The latter is reported more often and is akin to standardized regression coefficients.
fitMeasures(fit.8a)
resid(fit.8a, type="standardized")
modindices(fit.8a)
mi1<-modindices(fit.8a);print(mi1[mi1$mi>4.0,])
#Insignificant. P = 0.000

#improved model temp + NH4N + precipitation + func.richness
mod.8b <- 'func.struc ~ moisture + + NH4N + temp + TN + RA + copy_number + func.richness
           func.richness ~ moisture + temp + RA + copy_number + TN + bactr.richness
           bactr.struc ~ pH + TN + TC + copy_number
           bactr.richness ~ pH + temp + TN + NO3N
           copy_number ~ moisture + TC
           GPP ~ moisture + TC + precipitation + RA
           RS ~ func.struc + func.richness + GPP + TC + moisture + pH + NO3N'
fit.8b <- sem(mod.8b, data=semdat)#bactr.richness + NEE + temp + precipitation + moisture
fit.8b@SampleStats@cov
summary(fit.8b, rsq=T, standardize = T)
fitMeasures(fit.8b)

resid(fit.8b, type="standardized")
modindices(fit.8b)
mi1<-modindices(fit.8b);print(mi1[mi1$mi>4.0,])

anova(fit.8a,fit.8b)
#p = 0.000


# run diagnosis, remove insignificant links and add back suggested links

sink("SEM/log_0123.txt")
# Warm season SEM
mod.8a <- 'temp ~ warming
           moisture ~ warming
           NO3N ~ warming + pH + temp
           GPP ~ warming + pH
           pH ~ warming + temp
           copy_number ~ temp + moisture + func.struc + GPP
           func.struc ~ temp + moisture + NO3N
           RS ~ func.struc + GPP + moisture + temp + NO3N
           temp ~~ RS
           temp ~~ NO3N'

fit.8a <- sem(mod.8a, data=semdat.scaled.warm)
summary(fit.8a, rsq=T, standardize = T,fit.measures=TRUE)
fitMeasures(fit.8a)
resid(fit.8a, type="cor")
modindices(fit.8a)
mi1<-modindices(fit.8a);print(mi1[mi1$mi>3.5,])

#raw SEM plot
library(semPlot)
semPaths(fit.8a)
semPaths(fit.8a,"std",residuals=FALSE,nCharNodes=0,layout="tree2",edge.label.cex = 1,sizeMan =12)

# Cool season SEM
mod.8d <- 'temp ~ warming
           moisture ~ warming
           NO3N ~ warming + pH + temp
           GPP ~ warming + pH
           copy_number ~ temp + moisture + NO3N + func.struc
           func.struc ~ temp + moisture + GPP
           RS ~ func.struc + GPP + moisture + NO3N + temp
           temp ~~ func.struc
           temp ~~ NO3N
           moisture ~~ GPP'
#temp ~~ RS
fit.8d <- sem(mod.8d, data=semdat.scaled.cool)
summary(fit.8d, rsq=T, standardize = T,fit.measures=TRUE)
fitMeasures(fit.8d)
resid(fit.8d, type="cor")
modindices(fit.8d)
mi1<-modindices(fit.8d);print(mi1[mi1$mi>2,])

#raw SEM plot
library(semPlot)
semPaths(fit.8d)
semPaths(fit.8d,"std",residuals=FALSE,nCharNodes=0,layout="tree2",
         edge.label.cex = 1,sizeMan =12, nDigits = 3)


