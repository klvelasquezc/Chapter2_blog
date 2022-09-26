##%######################################################%##
#                                                          #
####                Panel de correlación                ####
#                                                          #
##%######################################################%##



# Performance Analyticis --------------------------------------------------

library(PerformanceAnalytics)

chart.Correlation(log(airpoll+1),
                  method="pearson",
                  histogram=TRUE,
                  pch=20)

?chart.Correlation


# Corrplot ----------------------------------------------------------------

library(corrplot)

corr <- round(cor(log(airpoll+1), method = "spearman"),2)
cor.mat <- cor.mtest(log(airpoll+1), conf.level = 0.95)


p1 <- corrplot(corr, method="color",  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat= cor.mat$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)$corrPos      
text(p1$x, p1$y, round(p1$corr, 2))
       


# ggcorrplot --------------------------------------------------------------

library(ggcorrplot)
         
ggcorrplot(corr, 
           type = "lower",
           lab = T, show.diag = F, 
           legend.title = " Pearson\nCorrelation", 
           colors= c("#BB4444", "#FFFFFF", "#4477AA"), 
           hc.order = T, 
           sig.level = 0.05, insig = "pch", pch=8, pch.cex = 2,  
           p.mat= cor.mat$p, 
           ggtheme = ggplot2::theme(
            panel.background = element_blank()))      



# GGally ------------------------------------------------------------------

library(GGally)

ggpairs(log(airpoll+1),
        upper = list(continuous= wrap("cor", method= "pearson", digits=2)),
        lower = list( continuous= "smooth"))
+theme_classic()
