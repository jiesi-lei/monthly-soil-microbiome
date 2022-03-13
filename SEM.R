source("basic_setting.R")
source("load_datasets.R")

library(lavaan)
library(AICcmodavg)
library(tidyverse)
library(semPlot)
library(vegan)

# source("lavaan.modavg.R")
# group
Group<-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 
Group = Group[match(rownames(gene),row.names(Group)),]
Group$Season = ordered(Group$Season, c("sp","sm","fl","wt"))
Group$Growing_season = ordered(Group$Growing_season, c("G","NG"))
Group$Period = ordered(Group$Period, c("Cold season 1","Warm season","Cold season 2"))
head(Group)
id_cool = which(Group$Period %in% c("Cold season 1","Cold season 2"))
id_warm = which(Group$Period == "Warm season" )

# preparation----
# functional profile
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

# copy number
copynum = read.csv("datasets/copy number/copy_num_R.csv",header = T,check.names = FALSE)
copy_number = copynum$copy_number

# Env
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
# write.table(semdat,"datasets/SEM/semdata.dat")
semdat %>% {scale(.)} %>% {data.frame(.,warming)}-> semdat.scaled

semdat.scaled.cool = semdat[id_cool,]
semdat.scaled.warm = semdat[id_warm,]

# write.csv(semdat.scaled.cool,"SEM/cool.csv",row.names = F)
# write.csv(semdat.scaled.warm,"SEM/warm.csv",row.names = F )

semdat.scaled.cool = read.csv("SEM/cool.csv")
semdat.scaled.warm = read.csv("SEM/warm.csv")

# SEM construction ----
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
mi1<-modindices(fit.8a);print(mi1[mi1$mi>3.5,]) # check potential links, adjust the above model structure iteratively

# raw SEM plot
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
fit.8d <- sem(mod.8d, data=semdat.scaled.cool)
summary(fit.8d, rsq=T, standardize = T,fit.measures=TRUE)
fitMeasures(fit.8d)
resid(fit.8d, type="cor")
modindices(fit.8d)
mi1<-modindices(fit.8d);print(mi1[mi1$mi>2,])# check potential links, adjust the above model structure iteratively

# raw SEM plot
library(semPlot)
semPaths(fit.8d)
semPaths(fit.8d,"std",residuals=FALSE,nCharNodes=0,layout="tree2",
         edge.label.cex = 1,sizeMan =12, nDigits = 3)


# partial correlation ----
library(ppcor)
pcor <- pcor(semdat, method = "spearman")
M = pcor$estimate
p.mat <- pcor$p.value
head(p.mat[, 1:5])

corrplot(M, type="upper", order="hclust", tl.col="black", tl.srt=45)
col <- colorRampPalette(c("#4477AA",  "#77AADD", "#FFFFFF","#EE9988", "#BB4444"))

pdf(file = "1.Monthly_warming_geochip_2012_manuscript/correlation.pdf")
corrplot(M, method="color", col=col(200), number.cex = .7,
         type="upper", order="original",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         addgrid.col = "grey",
         number.digits = 3
)
dev.off()
