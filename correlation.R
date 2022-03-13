# Correlation between rrn copy number & C deg genes
source("basic_setting.R")
source("load_datasets.R")

library(basicTrendline)
Group <-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 
Group = Group[match(rownames(gene),row.names(Group)),]
Group$Season = ordered(Group$Season, c("sp","sm","fl","wt"))
Group$Growing_season = ordered(Group$Growing_season, c("G","NG"))
Group$Period = ordered(Group$Period, c("Cold_season_1","Hot_season","Cold_season_2"))
head(Group)

# Env variables
env <- read.csv("monthly_ENV.csv",row.names = 1)
env = env[match(rownames(gene),row.names(env)),]

# Copy_number
copynum = read.csv("datasets/copy number/copy_num_R.csv",header = T,check.names = FALSE)


# C deg gene abundance
dat = read.csv("cut.treat.2.log.RA.csv",row.names = 1)
gene = as.data.frame(t(dat[,-c(1:9)]))
which((colSums(gene)) == 0)
info = dat[,c(1:9)]
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

env.all = data.frame(env,Group,copynum, C_deg_gene_abun)

env.all.C = filter(env.all, Warming == "C")
correl(env.all.C$copy_number,env.all.C$X3dayTemperatureAt7.5cm)
summary(lm(env.all.C$copy_number~env.all.C$X3dayTemperatureAt7.5cm))
env.all.W = filter(env.all, Warming == "W")
correl(env.all.W$copy_number,env.all.W$X3dayTemperatureAt7.5cm)
summary(lm(env.all.W$copy_number~env.all.W$X3dayTemperatureAt7.5cm))
correl(env.all$copy_number,env.all$X3dayTemperatureAt7.5cm)

trendline(env.all$X3dayTemperatureAt7.5cm, env.all$copy_number, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)

(ggplot(data = env.all, aes(X3dayTemperatureAt7.5cm,copy_number, fill = Warming)) 
  + geom_point(aes(fill = Warming, color = Warming))
  + geom_smooth(method = "lm", aes(color = Warming), alpha = 0.1)
  + theme_test()
  + theme(axis.text = element_text(color = "black",family = "Arial", size = 12),
          axis.title.y = element_text(color = "black",family = "Arial", size = 12),
          legend.position = "none",
          panel.border = element_rect(size = 1),
          axis.text.x = element_text(vjust = -1.6))
  + scale_y_continuous(expand = c(0, 0),breaks = c(0,0.5,1,1.5,2,2.5) ,limits = c(1,2.5))
  + scale_fill_manual(values = c("steelblue3","red4"))
  + scale_color_manual(values = c("steelblue3","red4"))
  + xlab("") + ylab(""))

# ggsave("figures/copynum/rrn_temp.pdf", dpi = 500, width = 3.8, height = 3)

df <- env.all %>% 
  group_by(Group) %>% 
  dplyr::summarise(GWC.avg = mean(GWC),
            GWC.se = std.error(GWC),
            pH.avg = mean(pH),
            pH.se = std.error(pH),
            NO3N.avg = mean(NO3N),
            NO3N.se = std.error(NO3N),
            NH4N.avg = mean(NH4N),
            NH4N.se = std.error(NH4N),
            TC.avg = mean(Tcpercent),
            TC.se = std.error(Tcpercent),
            TN.avg = mean(Tnpercent),
            TN.se = std.error(Tnpercent),
            RS.avg = mean(SoilResp),
            RS.se = std.error(SoilResp),
            RA.avg = mean(RootResp),
            RA.se = std.error(RootResp),
            RH.avg = mean(HeteroResp),
            RH.se = std.error(HeteroResp),
            ER.avg = mean(ER),
            ER.se = std.error(ER),
            NEE.avg = mean(NEE),
            NEE.se = std.error(NEE),
            GPP.avg = mean(GPP),
            GPP.se = std.error(GPP),
            Rain.avg = mean(Rain),
            Rain.se = std.error(Rain),
            temp.avg = mean(X3dayTemperatureAt7.5cm),
            temp.se = std.error(X3dayTemperatureAt7.5cm),
            copynum.avg = mean(copy_number),
            copynum.se = std.error(copy_number),
            C.gene.abun.avg = mean(C_deg_gene_abun),
            C.gene.abun.se = std.error(C_deg_gene_abun))
write.csv(df, "datasets/ENV/Env_all_summary.csv")



t.test(env.all[which(env.all$Warming.1=="C"),"copy_number"],
       env.all[which(env.all$Warming.1=="W"),"copy_number"], paired = T)
mean(env.all[which(env.all$Warming.1=="C"),"copy_number"])
mean(env.all[which(env.all$Warming.1=="W"),"copy_number"])

(ggplot(data = env.all, aes(Warming,copy_number, fill = Warming)) + geom_boxplot(outlier.shape=20,outlier.size=2, size = 0.5, color = "black")
  + theme_test()
  + theme(axis.text = element_text(color = "black",family = "Arial", size = 12),
        axis.title.y = element_text(color = "black",family = "Arial", size = 12),
        legend.position = "none",
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(vjust = -1.6))
  + scale_y_continuous(expand = c(0, 0),breaks = c(0,0.5,1,1.5,2,2.5) ,limits = c(1,2.5))
  + scale_fill_manual(values = c("steelblue3","red4"))
  + xlab("") + ylab(""))

# ggsave("figures/copynum/CW.pdf", dpi = 500, width = 3.8, height = 3)

(ggplot(data = env.all, aes(Period.1,copy_number, fill = Warming)) + geom_boxplot(outlier.size=1, size = 0.5, color = "black",position=position_dodge(0.85))
  + theme_test()
  + scale_x_discrete(limits=c("Cold season 1", "Warm season", "Cold season 2"))
  + theme(axis.text = element_text(color = "black",family = "Arial", size = 12),
          axis.title.y = element_text(color = "black",family = "Arial", size = 12),
          legend.position = "none",
          panel.border = element_rect(size = 1),
          axis.text.x = element_text(vjust = -1.6))
  + scale_y_continuous(expand = c(0, 0),breaks = c(0,0.5,1,1.5,2,2.5,3) ,limits = c(1,2.75))
  + scale_fill_manual(values = c("steelblue3","darkred"))
  + xlab("") + ylab(""))
# ggsave("figures/copynum/CW_period.pdf", dpi = 500, width = 3.8, height = 3)

copynum_rr = data.frame()
for (i in (1:12)){
  subcopynum = env.all$copy_number[c((8*i-7):(8*i))]
  for (j in (1:4)){
  copynum_rr[i,j] = log(subcopynum[j+4]/subcopynum[j],base = exp(1))
}
}
copynum_rr$mean = apply(copynum_rr,1,mean)
copynum_rr$se = apply(copynum_rr,1,std.error)
write.csv(copynum_rr, "datasets/ENV/copynum_rr.csv")

c_gene_rr = data.frame()
for (i in (1:12)){
  sub_c_gene = env.all$C_deg_gene_abun[c((8*i-7):(8*i))]
  for (j in (1:4)){
    c_gene_rr[i,j] = log(sub_c_gene[j+4]/sub_c_gene[j],base = exp(1))
  }
}
c_gene_rr$mean = apply(c_gene_rr,1,mean)
c_gene_rr$se = apply(c_gene_rr,1,std.error)
write.csv(c_gene_rr, "datasets/ENV/c_gene_rr.csv")

library("Hmisc")
corr = rcorr(as.matrix(env.all[,c(1:17,20,22)]))
corr$P
library(corrplot)
library(plotrix)
corr = cor(as.matrix(env.all[,c(1:17,20,22)]))
corrplot(corr$r, p.mat = corr$P, insig ="blank")

env.all.table = data.table::as.data.table(env.all)

a = as.data.frame(aggregate(env.all[,c(1:17,20,22)], by = list(env.all$Month), mean))
write.csv(a, "datasets/ENV/monthly_mean.csv")
corr = Hmisc::rcorr(as.matrix(a[,2:20]))
corrplot(corr$r, p.mat = corr$P, insig ="blank")

a = read.csv("datasets/ENV/correlation.csv", header = T, row.names = 1)
corr = Hmisc::rcorr(as.matrix(a), type = "pearson")
corrplot(corr$r, p.mat = corr$P, insig ="blank")
library(basicTrendline)
trendline(a$temp.avg, a$copy_num_rr, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)

library(ggthemes)
trendline(a$copy_num_rr, a$cate.rr.mean, model="line3P", ePos.x = "topleft", summary=TRUE, eDigit=5)
correl(a$copy_num_rr, a$cate.rr.mean, method = "spearman")
correl(a$copy_num_rr, a$temp.avg, method = "pearson")

# power test ----
library(pwr)
pwr.r.test(n = 12, r = 0.5454545, sig.level = 0.06661188, power = NULL,
           alternative = c("two.sided"))

trendline(a$copy_num_rr, a$probe.rr.mean, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
cor(a$copy_num_rr, a$probe.rr.mean)
(ggplot(a, aes(copy_num_rr, cate.rr.mean))
  #+ geom_smooth(method = "lm", color = "lightskyblue2", fill = "lightskyblue1", alpha = 0.3)
  + theme_test() +
    theme(axis.text = element_text(color = "black",family = "Arial", size = 12),
          axis.title.y = element_text(color = "black",family = "Arial", size = 12),
          panel.border = element_rect(size = 1),
          legend.position = "none",
          axis.text.x = element_text(vjust = -1.6))
  + scale_y_continuous(expand = c(0, 0),breaks = c(-0.2,-0.1,0, 0.1, 0.2, 0.3) ,limits = c(-0.15,0.3))
  + scale_x_continuous(expand = c(0, 0),breaks = c(-0.2, -0.1, 0, 0.1, 0.2) ,limits = c(-0.25,0.27))
  + xlab("") + ylab("") 
  + geom_errorbar(aes(x=copy_num_rr, ymin=cate.rr.mean-cate.rr.se, ymax=cate.rr.mean+cate.rr.se),width = 0.015, size = 0.2)
  + geom_errorbarh(aes(y=cate.rr.mean, xmin=copy_num_rr-copy_num_se, xmax=copy_num_rr+copy_num_se), height = 0.015, size = 0.2)
  + geom_point(shape = 21, size = 2.5, fill = "black", stroke = 0.2 ))

# ggsave("figures/correlation/rr_copynum_avg_1231.rr.pdf", dpi = 500, width = 3.8, height = 3)

trendline(a$temp.avg, a$copy_num_rr, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
(ggplot(a, aes(temp.avg, copy_num_rr))
  + geom_smooth(method = "lm", color = "lightskyblue2", fill = "lightskyblue1", alpha = 0.3)
  + theme_test() +
    theme(axis.text = element_text(color = "black",family = "Arial", size = 12),
          axis.title.y = element_text(color = "black",family = "Arial", size = 12),
          panel.border = element_rect(size = 1),
          legend.position = "none",
          axis.text.x = element_text(vjust = -1.6))
  + scale_y_continuous(expand = c(0, 0),breaks = c(-0.2,-0.1,0, 0.1, 0.2, 0.3) ,limits = c(-0.2,0.3))
  + scale_x_continuous(expand = c(0, 0),breaks = c(0, 5, 10, 15, 20, 25, 30) ,limits = c(3,30))
  + xlab("") + ylab("") 
  + geom_errorbar(aes(x=temp.avg, ymin=copy_num_rr-copy_num_se, ymax=copy_num_rr+copy_num_se),width = 0.5, size = 0.4)
  + geom_errorbarh(aes(y=copy_num_rr, xmin=temp.avg-temp.se, xmax=temp.avg+temp.se), height = 0.015, size = 0.4)
  + geom_point(shape = 21, size = 2.5, fill = "black", stroke = 0.2 ))

# ggsave("figures/correlation/temperature_rr_copynum.pdf", dpi = 500, width = 3.8, height = 3)

trendline(a$copy_number, a$avg.rr, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
trendline(a$copy_number, a$C_gene, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
trendline(a$rr_copynum, a$rr_C_abundance, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)

pairwise = read.csv("datasets/ENV/pairwise.csv", header = T)
trendline(pairwise$copy_num_rr, pairwise$C_deg_rr, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
trendline(pairwise$copy_num_rr, pairwise$X3dayTemperatureAt7.5cm, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)
trendline(a$rr_copynum, a$rr_C_abundance, model="line2P", ePos.x = "topleft", summary=TRUE, eDigit=5)


(ggplot(pairwise, aes(X3dayTemperatureAt7.5cm, copy_num_rr))
  + geom_smooth(method = "lm", color = "lightskyblue2", fill = "lightskyblue1", alpha = 0.3)
  + theme_test() +
    theme(axis.text = element_text(color = "black",family = "Arial", size = 10),
          axis.title.y = element_text(color = "black",family = "Arial", size = 10),
          panel.border = element_rect(size = 1),
          legend.position = "none",
          axis.text.x = element_text(vjust = -1.6))
  + scale_y_continuous(expand = c(0, 0),breaks = c(-0.4,-0.3,-0.2,-0.1,0, 0.1, 0.2, 0.3, 0.4, 0.5) ,limits = c(-0.4,0.5))
  + scale_x_continuous(expand = c(0, 0),breaks = c(0, 5, 10, 15, 20, 25, 30) ,limits = c(1.4,30))
  + xlab("") + ylab("") 
  + geom_point(shape = 21, size = 2.5, fill = "black", stroke = 0.2))
correl(pairwise$copy_num_rr, pairwise$X3dayTemperatureAt7.5cm, method = "pearson")
correl(pairwise$copy_num_rr, pairwise$C_deg_rr, method = "spearman")
# ggsave("figures/correlation/temperature_rr_copynum_0918.pdf", dpi = 500, width = 3.5, height = 2.8)



