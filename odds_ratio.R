# install.packages("oddsratio")
library(oddsratio)

major_odds_p<-data.frame()


for (cell in c("Exc","Inh",'Ast',"Oli","OPC",'Mic_Immune',"Vasc_Epithelia") ){

  for (region in c("All",'EC',"HC","TH","AG","MTC","PFC")){
    tmp<-Major_count5[Major_count5$Celltype==cell & Major_count5$BrainRegion==region,]
    tmp$Percentage2<-tmp$Percentage/100

    fit_glm <- glm(Percentage2 ~ Epigenetic_Information ,
                   data = tmp,
                   family = "quasibinomial")

    summ<-summary(fit_glm)
    pvalues<-as.data.frame(summ$coefficients)$"Pr(>|t|)" [2]

    coef_Epigenetic_Information <- coef(fit_glm)["Epigenetic_Information"]

    odds_ratio_Epigenetic_Information <- exp(coef_Epigenetic_Information)

    tmp_data<-data.frame(
      pathology=c('Epigenetic_Information'),
      celltype=cell,
      region=region,
      Odds_ratio=odds_ratio_Epigenetic_Information,
      pvalues=pvalues)

    major_odds_p<-rbind(major_odds_p,tmp_data)

  }
}

write.table(major_odds_p,"ATAC.Class.Epigenetic_Information.OddsRatio.tsv",sep="\t",quote = F,row.names = F)
