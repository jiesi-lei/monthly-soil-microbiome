source("0_setup.R")
source("0_Setup_gene.R")

env = read.csv("datasets/ENV/Env_all_summary.csv", header = T, row.names = 1)
env$GWC.avg = env$GWC.avg*100
env$GWC.se = env$GWC.se*100
env$Month = rep(1:12, each = 2)
env$Warming = rep(c("C","W"), times = 12)

var.all = c("GWC","pH","NO3N","NH4N",
            "TC","TN","RA","RS","RH",
            "ER","NEE",
            "GPP","Rain","temp",
            "DNA_yield_PG","MQ_PG")

for (i in 1:16){
  var = var.all[i]
  var.p = paste0(var,".avg")
  ymin.p = paste0(var.p,"-",var,".se")
  ymax.p = paste0(var.p,"+",var,".se")
  (ggplot(env, aes_string(x = "Month", y = var.p ,color = "black", fill = "Warming"))+
      #geom_vline(xintercept = c(2.5,5.5,8.5,11.5), linetype = "dashed")+
      geom_errorbar(aes_string(ymin = ymin.p, ymax = ymax.p ),width=.2,color = "black")+
      geom_line(aes_string(color = "Warming"))+
      geom_point(shape = 21, size=4, color = "black")+
      scale_fill_manual(values = c("royalblue4","firebrick4"))+
      scale_color_manual(values = c("royalblue4","firebrick4"))+
      scale_x_continuous(limits = c(0.75,12.25), breaks = c(1:12), expand = c(0,0))+
      scale_y_continuous(expand = c(0,0), limits = c(0,0.13)) +  
      expand_limits(x = 0) +
      ylab(var) +
      theme_test()+
      theme(axis.text = element_text(color = "black",family = "Arial", size = 20),
            axis.title.y = element_text(color = "black",family = "Arial", size = 20),
            legend.position = "none",
            axis.title.x = element_text(color = "black",family = "Arial", size = 20),
            panel.border = element_rect(size = 1.5)))
  
  filename = paste0("figures/ENV/2012/",var,".month.pdf")
  ggsave(filename, dpi = 500, width = 8, height = 4)
  i = i+1
}

#Env variables
env <- read.csv("monthly_ENV_Y2209.csv",row.names = 1)
env = env[,-c(2,8,9,15,19,20)]
Group<-read.csv("monthly_SampleList.csv",row.names = 1,stringsAsFactors = TRUE, check.names = FALSE) 
env.all = data.frame(env,Group)
env.all$Period_trt = paste0(Group$Period,"_",Group$Warming)
df <- env.all %>% 
  group_by(Period_trt) %>% 
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
                   DNA_yield_ND.avg = mean(DNA_yield_ND),
                   DNA_yield_ND.se = std.error(DNA_yield_ND),
                   DNA_yield_PG.avg = mean(DNA_yield_PG),
                   DNA_yield_PG.se = std.error(DNA_yield_PG),
                   MQ_ND.avg = mean(MQ_ND),
                   MQ_ND.se = std.error(MQ_ND),
                   MQ_PG.avg = mean(MQ_PG),
                   MQ_PG.se = std.error(MQ_PG)
                   )
df <- env.all %>% 
  group_by(Warming) %>% 
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
                   DNA_yield_ND.avg = mean(DNA_yield_ND),
                   DNA_yield_ND.se = std.error(DNA_yield_ND),
                   DNA_yield_PG.avg = mean(DNA_yield_PG),
                   DNA_yield_PG.se = std.error(DNA_yield_PG),
                   MQ_ND.avg = mean(MQ_ND),
                   MQ_ND.se = std.error(MQ_ND),
                   MQ_PG.avg = mean(MQ_PG),
                   MQ_PG.se = std.error(MQ_PG)
  )

#write.csv(df,"datasets/ENV/period_summary.csv",row.names = F)

# names(env.all)
for (i in (1:19)){
print((colnames(env.all))[i])
df.1 <- env.all[,i]
df = data.frame(df.1,env.all$Warming, env.all$Period)
cold.1_control = df[which(df$env.all.Period == "Cold season 1"& df$env.all.Warming == "C"),1]
cold.1_warming = df[which(df$env.all.Period == "Cold season 1"& df$env.all.Warming == "W"),1]
t.test(cold.1_control,cold.1_warming, paired = T)
warm_control = df[which(df$env.all.Period == "Warm season"& df$env.all.Warming == "C"),1]
warm_warming = df[which(df$env.all.Period == "Warm season"& df$env.all.Warming == "W"),1]
t.test(warm_control,warm_warming, paired = T)
cold.2_control = df[which(df$env.all.Period == "Cold season 2"& df$env.all.Warming == "C"),1]
cold.2_warming = df[which(df$env.all.Period == "Cold season 2"& df$env.all.Warming == "W"),1]
t.test(cold.2_control,cold.2_warming, paired = T)
control = df[which(df$env.all.Warming == "C"),1]
warming = df[which(df$env.all.Warming == "W"),1]
t.test(control,warming, paired = T)
i = i+1
}

sink("datasets/ENV/monthly_t_test_Y220912.txt")
env.all$Month = rep(1:12, each = 8)
for (i in (1:19)){
  print((colnames(env.all))[i])
  for (j in 1:12){
  print(paste("Month",j))
  df.1 <- env.all[,i]
  df = data.frame(df.1,env.all$Warming, env.all$Month)
  control = df[which(df$env.all.Month == j & df$env.all.Warming == "C"),1]
  warming = df[which(df$env.all.Month == j & df$env.all.Warming == "W"),1]
  print(t.test(warming,control, paired = T))
}
}
sink()

Result.S <- c()
for (i in (1:length(colnames(env)))){
print(colnames(env)[i])
aov = aov(env[,i]~Warming*Period+Block, data=Group)
print(summary(aov))
df = df.residual(aov)
MSerror<-deviance(aov)/df
com1.S=LSD.test(aov, "Period", df,MSerror, p.adj="holm", group=TRUE,alpha = 0.05)
print(com1.S$groups)

a.S = data.frame(com1.S$means[,1:2],rownames(com1.S$means))
b.S = data.frame(com1.S$groups)
c.S = merge(a.S,b.S)
names(c.S)[1] <- "value"
c.S$env <- colnames(env)[i]
Result.S = rbind(Result.S,c.S)}

write.csv(Result.S,"Env.season.LSD.2012.csv",row.names = F)

