source("basic_setting.R")
source("load_datasets.R")

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

