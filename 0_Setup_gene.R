#gene data
dat = read.csv("cut.treat.2.log.RA.csv",row.names = 1)
gene = as.data.frame(t(dat[,-c(1:9)]))
which((colSums(gene)) == 0)
info = dat[,c(1:9)]

dat1 = read.csv("cut.treat.2.log.RA.M1to2.csv",row.names = 1)
gene1 = as.data.frame(t(dat1[,-c(1:9)]))
#gene1 = gene1[,-which((colSums(gene1)) == 0)]

dat2 = read.csv("cut.treat.2.log.RA.M3to9.csv",row.names = 1)
gene2 = as.data.frame(t(dat2[,-c(1:9)]))
#gene2 = gene2[,-which((colSums(gene2)) == 0)]

dat3 = read.csv("cut.treat.2.log.RA.M10to12.csv",row.names = 1)
gene3 = as.data.frame(t(dat3[,-c(1:9)]))
#gene3 = gene3[,-which((colSums(gene3)) == 0)]

rm(dat, dat1, dat2, dat3)

dat = read.delim("resample_zotutab_new_21567.txt",row.names = 1, header = T,check.names = TRUE, sep = "\t")
dat = dat[-which(rowSums(dat)== 0),]
dat = as.data.frame(t(dat))
asv1 = dat[1:16,]
asv2 = dat[17:72,]
asv3 = dat[73:96,]
write.csv(as.data.frame(rowSums(dat)),"datasets/copy number/otu_abundance.csv",row.names = T)

#group
Group<-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 
Group = Group[match(rownames(gene),row.names(Group)),]
Group$Season = ordered(Group$Season, c("sp","sm","fl","wt"))
Group$Growing_season = ordered(Group$Growing_season, c("G","NG"))
Group$Period = ordered(Group$Period, c("Cold season 1","Warm season","Cold season 2"))
head(Group)

#Env variables
env <- read.csv("monthly_ENV.csv",row.names = 1)
env = env[match(rownames(gene),row.names(env)),]
env = env[match(rownames(dat),row.names(env)),]
