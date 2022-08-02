
#' Identifies which variants are statistically significant in x number of pools.
#'
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param altFreq1 numeric. The alternate allele frequency for bulk A.
#' @param altFreq2 numeric. The alternate allele frequency for bulk B.
#' @param replicates integer. The number of bootstrap replications.
#' @param filter numeric. an optional minimum SNP-index filter
#'
#' @return Returns a vector of length replicates delta SNP-indeces
#'
plotMakerBSA=function(datum, xplot=datum$binMiddle/1000, yplot=datum$All_pools){
  ggplot(data=datum,aes(x=xplot, y=yplot))+
    geom_pointrange(data=temp,
                    aes(ymin=pool_min,
                        ymax=pool_max ,
                        x=binMiddle/1000,
                        y=All_pools,
                        size=combine,
                        colour=combine,
                        alpha=combine,
                        order="combine"))+
    scale_size_manual(name="Condition and significance",
                      labels=c("Lung, not significant","YPD, not significant","Lung, significant",  "YPD, significant"),
                      values = c(0.1,0.1,0.4, 0.4), drop=F)+
    scale_colour_manual(name="Condition and significance",
                        labels=c("Lung, not significant","YPD, not significant","Lung, significant",  "YPD, significant"),
                        values = c('FALSE_Lungs'="red",'FALSE_YPD'="black",
                                   'TRUE_Lungs'="red",'TRUE_YPD'="black"), drop=F)+
    scale_alpha_manual(name="Condition and significance",
                       labels=c("Lung, not significant","YPD, not significant","Lung, significant",  "YPD, significant"),
                       values = c(0.95,0.95,1, 1), drop=F)+
    #geom_vline(xintercept  = IR_list[[paste0("chr",i)]],color="black", size=0.3, linetype="dashed")+
    #geom_vline(xintercept = 0, size=0.3)+
    labs(title=paste0("Chromosome ", i), x=paste0("Position on chromosome ",i," (kb)"),y="Change in allele frequency (BSA)", color="Pool")+
    geom_hline(yintercept = c(0), color="black", size=0.6)+
    scale_x_continuous(breaks = seq(from=0, to=2400, by=300), limits=c(0, 2400),expand = c(0, 0))+
    scale_y_continuous(breaks = c(-0.5, -0.25, 0,0.25,0.5), limits = c(-0.55,0.55), expand = c(0,0))+
    theme(#legend.title = element_text(size=18),
      legend.position = "right",
      plot.title = element_blank(),
      legend.text = element_text(size=16),
      axis.title=element_text(size=32),
      axis.title.y =element_text(margin = margin(r=25, t=0, l=5, b=0)),
      axis.title.x =element_text(margin = margin(r=0, t=25, l=0, b=5)),
      axis.text=element_text(size = 28),
      axis.ticks.length = unit(0.3,"cm"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_rect(fill = "white",colour = "black", size=1),
      text=element_text(family="sans"))

}

