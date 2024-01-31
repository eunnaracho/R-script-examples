#ggplot2 example

ggsave(filename = paste0(chemical,"_Simple mut spectra_proportion.jpeg"),
       
       ggplot(data = MF_subtypes, 
              mapping = aes(x = subtype,
                            y = proportion,
                            fill = dose,
                            pattern = sample)) + 
         geom_bar(stat="identity", 
                  width=0.75, 
                  position=position_dodge(width=0.9),
                  colour="black") +
         scale_fill_manual(values=c("darkslategray1",
                                    "deepskyblue",
                                    "blue2",
                                    "darkblue")) +
         labs(x = "Mutation Subtype", 
              y = "Proportion (%)",
              fill = paste0(chemical, " Concentration\n","(", conc_unit, ")")) +
         scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
         theme_classic() +
         theme(legend.position = c(0.85, 0.75),
               axis.text.x = element_text(colour = "black", size=12),
               axis.text.y = element_text(colour = "black", size=12),
               axis.title.x = element_text(margin=margin(t=10), size=14),
               axis.title.y = element_text(margin=margin(r=10), size=14)),
       
       width = 18, height = 12.5, units = "cm")
