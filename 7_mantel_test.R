#Mantel code from pipeline
source("0_Setup_gene.R")
source("0_Setup.R")

index = c(1:3)
season = c("cold1","warm","cold2")
env_dat = env[-c(8,14,18,19)]

#Geochip
#For different seasons (6 in total)
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

#ASV
#For different seasons (3 in total)
for (j in 1:3){
  
  fg.dat = as.data.frame(t(get(paste0("asv",index[j]))))
  
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
  write.csv(report,file=paste0("datasets/mantel/mantel_2012_ASV_",season[j],".csv"),row.names=FALSE)
  
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
  write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_ASV_",season[j],".csv"),row.names=FALSE)
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
write.csv(report,file=paste0("datasets/mantel/mantel_2012_all_ASV.csv"),row.names=FALSE)

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
write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_all_ASV.csv"),row.names=FALSE)


#For warming and control - all C degradation genes combined
subselect = which(info[,1]%in% c.deg.select)
subinformation<-info[subselect,]
fg.dat<-as.data.frame(t(gene))[subselect,]
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
write.csv(report,file=paste0("datasets/mantel/mantel_2012_C_deg_all.csv"),row.names=FALSE)

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
write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_C_deg_all.csv"),row.names=FALSE)

#==============================
#Gene category & env variables
index = c(1:4)
season = c("sp","sm","fl","wt")
env_dat = env[-c(8,14,18,19)]
control = c(1:4,9:12,17:20)
warming = c(5:8,13:16,21:24)
report.mantel =list()
report.pmantel =list()

#For different seasons (4 in total)
for (j in 1:4){
  for (h in 1:length(c.deg.select)){
  print(c.deg.select[h])
  subselect=which(info[,"gene"]==c.deg.select[h])
  subinformation<-info[subselect,]
  fg.dat<-as.data.frame(t(get(paste0("gene",index[j]))))[subselect,]
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
  report.mantel[[h]] = data.frame()
  for(i in 1:ncol(env.std)){
    envdis =vegdist(env.std[,i],method = "euclidean", na.rm=T)
    mantel.BC = mantel(envdis, BC.beta, na.rm=T) 
    mantel.JC = mantel(envdis, JC.beta, na.rm=T)
    report.mantel[[h]] = rbind(report.mantel[[h]],c(colnames(env.select)[i], mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic, mantel.JC$signif))
  }
  colnames(report.mantel[[h]])<- c("Envs",paste(rep(c("r","p"),2),rep(c("BC","JC"),each=2),sep = "."))
  #write.csv(report,file=paste0("datasets/mantel/mantel_2012_",season[j],c.deg.select[h],"_control.csv"),row.names=FALSE)

  
  #calculate partial mantel between community and each env_factor
  report.pmantel[[h]] = data.frame()
  for(i in 1:ncol(env.std)){
    envdis = dist(env.std[,i])  
    envdis2 = dist(env.std[,-i])  
    pmantel.BC = mantel.partial(BC.beta, envdis, envdis2, na.rm=T)
    pmantel.JC = mantel.partial(JC.beta, envdis, envdis2, na.rm=T)  
    pmantel.EU = mantel.partial(EU.beta, envdis, envdis2, na.rm=T)  
    report.pmantel[[h]] = rbind(report.pmantel[[h]],c(colnames(env.select)[i], pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic, pmantel.JC$signif,
                            pmantel.EU$statistic, pmantel.EU$signif))
  }
  colnames(report.pmantel[[h]])<- c("Envs",paste(rep(c("r","p"),3),rep(c("BC","JC","EU"),each=2),sep = "."))
  #write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_",season[j],".csv"),row.names=FALSE)
  }
  save(report.mantel,file=paste0("datasets/mantel/mantel_2012_",season[j],"_all.Rdata"))
  save(report.pmantel,file=paste0("datasets/mantel/p_mantel_2012_",season[j],"_all.Rdata"))
}


#All timepoints - ASV counts ===================================================
fg.dat = read.delim("Galaxy61-[resample_ASV_counts_minsize3_n.txt].txt",row.names = 1, header = T,check.names = TRUE, sep = "\t")
fg.dat = as.data.frame(t(fg.dat))

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
write.csv(report,file=paste0("datasets/mantel/mantel_2012_all_ASV.csv"),row.names=FALSE)

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
write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_all_ASV.csv"),row.names=FALSE)

#----------------------------
#mantel的深入分析
load("datasets/mantel/mantel_2012_sp_control.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_sp_control.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_sp_control.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_sp_control.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_sm_control.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_sm_control.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_sm_control.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_sm_control.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_fl_control.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_fl_control.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_fl_control.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_fl_control.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_wt_control.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_wt_control.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_wt_control.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_wt_control.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_sp_warming.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_sp_warming.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_sp_warming.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_sp_warming.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_sm_warming.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_sm_warming.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_sm_warming.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_sm_warming.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_fl_warming.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_fl_warming.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_fl_warming.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_fl_warming.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_wt_warming.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_wt_warming.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_wt_warming.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_wt_warming.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_sp_all.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_sp_all.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_sp_all.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_sp_all.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_sm_all.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_sm_all.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_sm_all.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_sm_all.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_fl_all.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_fl_all.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_fl_all.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_fl_all.xlsx",row.names = T)

load("datasets/mantel/mantel_2012_wt_all.Rdata")
write.xlsx(report.mantel, "datasets/mantel/mantel_2012_wt_all.xlsx",row.names = T)
load("datasets/mantel/p_mantel_2012_wt_all.Rdata")
write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_wt_all.xlsx",row.names = T)


#For different treatment

report.mantel =list()
report.pmantel =list()

  for (h in 1:length(c.deg.select)){
    print(c.deg.select[h])
    subselect=which(info[,"gene"]==c.deg.select[h])
    subinformation<-info[subselect,]
    fg.dat<-as.data.frame(t(gene))[subselect,warming]
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
    report.mantel[[h]] = data.frame()
    for(i in 1:ncol(env.std)){
      envdis =vegdist(env.std[,i],method = "euclidean", na.rm=T)
      mantel.BC = mantel(envdis, BC.beta, na.rm=T) 
      mantel.JC = mantel(envdis, JC.beta, na.rm=T)
      report.mantel[[h]] = rbind(report.mantel[[h]],c(colnames(env.select)[i], mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic, mantel.JC$signif))
    }
    colnames(report.mantel[[h]])<- c("Envs",paste(rep(c("r","p"),2),rep(c("BC","JC"),each=2),sep = "."))
    #write.csv(report,file=paste0("datasets/mantel/mantel_2012_",season[j],c.deg.select[h],"_control.csv"),row.names=FALSE)
    
    
    #calculate partial mantel between community and each env_factor
    report.pmantel[[h]] = data.frame()
    for(i in 1:ncol(env.std)){
      envdis = dist(env.std[,i])  
      envdis2 = dist(env.std[,-i])  
      pmantel.BC = mantel.partial(BC.beta, envdis, envdis2, na.rm=T)
      pmantel.JC = mantel.partial(JC.beta, envdis, envdis2, na.rm=T)  
      pmantel.EU = mantel.partial(EU.beta, envdis, envdis2, na.rm=T)  
      report.pmantel[[h]] = rbind(report.pmantel[[h]],c(colnames(env.select)[i], pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic, pmantel.JC$signif,
                                                        pmantel.EU$statistic, pmantel.EU$signif))
    }
    colnames(report.pmantel[[h]])<- c("Envs",paste(rep(c("r","p"),3),rep(c("BC","JC","EU"),each=2),sep = "."))
    #write.csv(report,file=paste0("datasets/mantel/p_mantel_2012_",season[j],".csv"),row.names=FALSE)
  }
  save(report.mantel,file=paste0("datasets/mantel/mantel_2012_warming.Rdata"))
  save(report.pmantel,file=paste0("datasets/mantel/p_mantel_2012_warming.Rdata"))

  load("datasets/mantel/mantel_2012_warming.Rdata")
  write.xlsx(report.mantel, "datasets/mantel/mantel_2012_warming.xlsx",row.names = T)
  load("datasets/mantel/p_mantel_2012_warming.Rdata")
  write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_warming.xlsx",row.names = T)
  
  load("datasets/mantel/mantel_2012_control.Rdata")
  write.xlsx(report.mantel, "datasets/mantel/mantel_2012_control.xlsx",row.names = T)
  load("datasets/mantel/p_mantel_2012_control.Rdata")
  write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_control.xlsx",row.names = T)
  
  load("datasets/mantel/mantel_2012_all.Rdata")
  write.xlsx(report.mantel, "datasets/mantel/mantel_2012_all.xlsx",row.names = T)
  load("datasets/mantel/p_mantel_2012_all.Rdata")
  write.xlsx(report.pmantel, "datasets/mantel/p_mantel_2012_all.xlsx",row.names = T)

  #p.adjust
  p = c(0.008,0.011,0.017,0.002,0.006,0.007,0.007,0.007,0.032,0.005,0.024,0.012,0.009,
        0.007,0.004,0.005,0.030,0.034,0.061,0.008)
  p.adjust(p, method = "BH")
  p.adjust(p, method = "bonferroni")
  p.adjust(p, method = "fdr")
  p.adjust(p, method = "holm")
  p.adjust(p, method = "hochberg")
  p.adjust(p, method = "hommel")
  #library(fdrtool)
  #FDR = fdrtool(p, statistic = "pvalue",plot = F)$qval
  
  p = c(0.013,0.013,0.006,0.011,0.007,0.003,0.008,0.007,0.005,0.008,
        0.019,0.012,0.008,0.010,0.022,0.013,0.011,0.013,0.005,0.013)
  p.adjust(p, method = "fdr")
  p.adjust(p, method = "bonferroni")
  
  p = c(0.042, 0.03, 0.021, 0.033, 0.016, 0.035, 0.045, 0.034, 0.074, 0.056, 0.078, 0.043, 0.044, 0.059, 0.006, 0.056, 0.046, 0.081, 0.089,0.042)
  p.adjust(p, method = "fdr")
  p.adjust(p, method = "bonferroni")
  
  p = c(0.013,0.003,0.003,0.003,0.004,0.002,0.003,0.002,0.002,0.006,0.001,0.001,0.001,0.001,0.003,0.005,0.002,0.004,0.003,0.013)
  p.adjust(p, method = "fdr")
  p.adjust(p, method = "bonferroni")
  
  p = c(0.405,0.423,0.477,0.272,0.374,0.371,0.375,0.406,0.469,0.355,0.382,0.438,0.421,0.374,0.392,0.403,0.46,0.411,0.505,0.405)
  p.adjust(p, method = "fdr")
  p.adjust(p, method = "bonferroni")
  
  p = c(0.517,0.532,0.442,0.381,0.448,0.476,0.55,0.533,0.471,0.408,0.521,0.572,0.536,0.501,0.541,0.642,0.542,0.52,0.406,0.517)
  p.adjust(p, method = "fdr")
  p.adjust(p, method = "bonferroni")
  