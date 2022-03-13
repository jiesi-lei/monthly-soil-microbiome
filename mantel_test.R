#Mantel linking C deg genes with respiration
source("Basic_setting.R")
source("Load_datasets.R")

index = c(1:3)
season = c("cold1","warm","cold2")
env_dat = env[-c(8,14,18,19)]

# Geochip
# For different season and treatment (6 in total)
for (j in 1:3){
  
fg.dat = as.data.frame(t(get(paste0("gene",index[j]))))

fg.dat[is.na(fg.dat)] = 0

# match env and fg datasets
samp.fg = colnames(fg.dat)
samp.env= rownames(env_dat)
my.fg = match(samp.fg, samp.env)
env_dat2 = env_dat[my.fg,]

if(is.na(which(is.na(my.fg))[1] > 0)){	
  fg.dat2 = fg.dat 
}else{
  a <- which(is.na(my.fg))
  fg.dat2 <- fg.dat[,-a]
  env_dat2 = env_dat2[-a,]
}

BC.beta = vegdist(t(fg.dat2), method="bray") 
JC.beta = vegdist(t(fg.dat2), method="jaccard",binary=T)
EU.beta = vegdist(t(fg.dat2), method="euclidean",binary=T)

env.select = env_dat2
env.std = decostand(env.select, method = "standardize", MARGIN=2)
envdis = dist(env.std)

#calculate the relationship between community and each env_factor
report =c()
for(i in 1:ncol(env.std)){
  envdis =vegdist(env.std[,i],method = "euclidean", na.rm=T)
  mantel.BC = mantel(envdis, BC.beta, na.rm=T) 
  mantel.JC = mantel(envdis, JC.beta, na.rm=T)
  report = rbind(report,c(colnames(env.select)[i], mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic, mantel.JC$signif))
}
colnames(report)<- c("Envs",paste(rep(c("r","p"),2),rep(c("BC","JC"),each=2),sep = "."))
write.csv(report,file=paste0("datasets/mantel/mantel_2012_",season[j],".csv"),row.names=FALSE)

#calculate partial mantel between community and each env_factor
report = c()
for(i in 1:ncol(env.std)){
  envdis = dist(env.std[,i])  
  envdis2 = dist(env.std[,-i])  
  pmantel.BC = mantel.partial(BC.beta, envdis, envdis2, na.rm=T)
  pmantel.JC = mantel.partial(JC.beta, envdis, envdis2, na.rm=T)  
  pmantel.EU = mantel.partial(EU.beta, envdis, envdis2, na.rm=T)  
  report = rbind(report,c(colnames(env.select)[i], pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic, pmantel.JC$signif,
                          pmantel.EU$statistic, pmantel.EU$signif))
}
colnames(report)<- c("Envs",paste(rep(c("r","p"),3),rep(c("BC","JC","EU"),each=2),sep = "."))
write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_",season[j],".csv"),row.names=FALSE)
}

#All timepoints
fg.dat = as.data.frame(t(gene))

fg.dat[is.na(fg.dat)] = 0

# match env and fg datasets
samp.fg = colnames(fg.dat)
samp.env= rownames(env_dat)
my.fg = match(samp.fg, samp.env)
env_dat2 = env_dat[my.fg,]

if(is.na(which(is.na(my.fg))[1] > 0)){	
  fg.dat2 = fg.dat 
}else{
  a <- which(is.na(my.fg))
  fg.dat2 <- fg.dat[,-a]
  env_dat2 = env_dat2[-a,]
}

BC.beta = vegdist(t(fg.dat2), method="bray") 
JC.beta = vegdist(t(fg.dat2), method="jaccard",binary=T)
EU.beta = vegdist(t(fg.dat2), method="euclidean",binary=T)

env.select = env_dat2
env.std = decostand(env.select, method = "standardize", MARGIN=2)
envdis = dist(env.std)

#calculate the relationship between community and each env_factor
report =c()
for(i in 1:ncol(env.std)){
  envdis =vegdist(env.std[,i],method = "euclidean", na.rm=T)
  mantel.BC = mantel(envdis, BC.beta, na.rm=T) 
  mantel.JC = mantel(envdis, JC.beta, na.rm=T)
  mantel.EU = mantel(envdis, EU.beta, na.rm=T)
  report = rbind(report,c(colnames(env.select)[i], mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic, mantel.JC$signif,
                          mantel.EU$statistic, mantel.EU$signif))
}
colnames(report)<- c("Envs",paste(rep(c("r","p"),3),rep(c("BC","JC","EU"),each=2),sep = "."))
write.csv(report,file=paste0("datasets/mantel/mantel_2012_all.csv"),row.names=FALSE)

#calculate partial mantel between community and each env_factor
report = c()
for(i in 1:ncol(env.std)){
  envdis = dist(env.std[,i])  
  envdis2 = dist(env.std[,-i])  
  pmantel.BC = mantel.partial(BC.beta, envdis, envdis2, na.rm=T)
  pmantel.JC = mantel.partial(JC.beta, envdis, envdis2, na.rm=T)
  pmantel.EU = mantel.partial(EU.beta,envdis, envdis2, na.rm=T)
  report = rbind(report,c(colnames(env.select)[i], pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic, pmantel.JC$signif,
                          pmantel.EU$statistic, pmantel.EU$signif))
}
colnames(report)<- c("Envs",paste(rep(c("r","p"),3),rep(c("BC","JC","EU"),each=2),sep = "."))
write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_all.csv"),row.names=FALSE)

