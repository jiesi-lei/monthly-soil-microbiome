source("0_setup.R")
source("0_setup_gene.R")

c.deg.select = c("amyA","apu","cda","glucoamylase","nplT","pula",#Starch
                 "ara","mannanase","xyla","xylanase",#Hemicellulose
                 "cellobiase","endoglucanase","exoglucanase",#Cellulous
                 "chitinase","chitin_deacetylase_fungi",#Chitin
                 "pectinase (pectate_lyase)","pectin_lyase_Oomycetes",#Pectin
                 "glx","mnp","ligninase","phenol_oxidase",#lignin
                 "cutinase",#cutin
                 "cdh","lmo")#Terpene
N.select = c("gdh","glnA_fungi","urec",#Ammonification
             "amoa","amoa_quasi","hao",#Nitrification
             "narg","nirk","nirs","norb","nosz",#Denitrification
             "napa","nrfa",#Dissimilatory N reduction
             "nasa","narb","NiR",#assimilatory N reduction
             "nifh",#N fixation
             "hzsa"#anammox, gene hzo close to 0
)
Methane.select = c("hdrB","Mch_methane","mcra","mrtH","MT2",#Methanogenesis
                   "mmox","pmoa")#Methane oxidation
C.fixation.select = c("rubisco","pcc","codh","aclb")

subselect = gene[,info$Gene.name%in%c(c.deg.select,Methane.select,C.fixation.select)]
subselect = gene[,info$Gene.name%in%c(c.deg.select)]#,N.select)]
subselect = gene[,5000:10000]
  
#for GeoChip data - ADONIS
ado1<-adonis(gene ~ Warming*Month+Block, data=Group, distance = "bray")
ado1
Group1 = filter(Group, Period == "Cold season 1")
ado1<-adonis(gene1 ~ Warming+Block, data=Group1, distance = "bray")
ado1
Group2 = filter(Group, Period == "Warm season")
ado1<-adonis(gene2 ~ Warming+Block, data=Group2, distance = "bray")
ado1
Group3 = filter(Group, Period == "Cold season 2")
ado1<-adonis(gene3 ~ Warming+Block, data=Group3, distance = "bray")
ado1

mrpp1 <- mrpp(gene, grouping = Group$Warming)
mrpp1
mrpp1 <- mrpp(gene, grouping = Group$Month)
mrpp1
mrpp1 <- mrpp(gene, grouping = Group$Period:Group$Warming)
mrpp1

anosim1 <- anosim(gene, grouping = Group$Warming)
anosim1
anosim1 <- anosim(gene, grouping = Group$Month)
anosim1


#For 16s data - ADONIS
dat = read.delim("resample_zotutab_new_21567.txt",row.names = 1, header = T,check.names = TRUE, sep = "\t")
dat = dat[-which(rowSums(dat)== 0),]
dat = data.frame(t(dat))
Group<-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 
Group = Group[match(rownames(dat),row.names(Group)),]
Group$Season = ordered(Group$Season, c("sp","sm","fl","wt"))

ado1<-adonis(dat ~ Warming*Month+Block, data=Group, distance = "bray")
ado1
Group1 = filter(Group, Period == "Cold season 1")
ado1<-adonis(asv1 ~ Warming*Month+Block, data=Group1, distance = "bray")
ado1
Group2 = filter(Group, Period == "Warm season")
ado1<-adonis(asv2 ~ Warming+Block, data=Group2, distance = "bray")
ado1
Group3 = filter(Group, Period == "Cold season 2")
ado1<-adonis(asv3 ~ Warming+Block, data=Group3, distance = "bray")
ado1

mrpp1 <- mrpp(dat, grouping = Group$Warming)
mrpp1
mrpp1 <- mrpp(dat, grouping = Group$Month)
mrpp1

anosim1 <- anosim(dat, grouping = Group$Warming)
anosim1
anosim1 <- anosim(dat, grouping = Group$Month)
anosim1

