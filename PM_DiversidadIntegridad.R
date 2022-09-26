### Análisis Diversidad
##  Pequeños mamiferos y áreas con cambio
##  de integridad ecológica


## Librerias

library(iNEXT)
library(AICcmodavg)
library(ggeffects)
library(DHARMa)
library(performance)
library(patchwork)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(ggpubr)
library(ggdist)
library(parameters)
library(hillR)

# Preparación datos -----------------------------------------------------------------

PM <- read.csv("Roedores_doctorado110922.csv", sep=",", header = T) %>% 
  filter(!grepl('Recaptura', Obs)) %>% 
  count(Nodo, Epiteto) %>% 
  mutate(Especies= case_when(
    Epiteto== "irroratus"~ "L. irroratus",
    Epiteto== "hispidus" ~ "S. hispidus",
    Epiteto== "parva" ~ "C. parva",
    Epiteto== "fulvescens" ~ "R. fulvescens",
    Epiteto== "morfo 3" ~ "P. morfo 3",
    Epiteto== "morfo 1" ~ "P. morfo 1",
    Epiteto== "morfo 2" ~ "P. morfo 2",
    Epiteto== "morfo 3" ~ "P. morfo 3",
    Epiteto== "" ~ "Peromyscus sp.")) %>% 
  pivot_wider(names_from = Nodo, values_from = n) %>% 
  select(-Epiteto,) %>%
  replace(is.na(.),0) %>%  
  column_to_rownames(var = "Especies")

# Taxonomic diversity -----------------------------------------------------

q0 <- hill_taxa(PM, q=0, MARGIN = 2)
q1 <- hill_taxa(PM, q=1, MARGIN = 2)
q2 <- hill_taxa(PM, q=2, MARGIN = 2)

hill_div <- data.frame(q0=q0, q1=q1, q2=q2)


X <- iNEXT(PM, q=c(0,1,2), datatype="abundance", )

Nodo57 <- X$iNextEst$`57` %>% 
  filter(method== "observed")

Nodo58 <- X$iNextEst$`58` %>% 
  filter(method== "observed")

Nodo60 <- X$iNextEst$`60`%>% 
  filter(method== "observed")

Nodo61 <- X$iNextEst$`61`%>% 
  filter(method== "observed")

Nodo62 <- X$iNextEst$`62`%>% 
  filter(method== "observed")

Nodo63 <- X$iNextEst$`63`%>% 
  filter(method== "observed")

Nodo64 <- X$iNextEst$`64`%>% 
  filter(method== "observed")

Nodo65 <- X$iNextEst$`65`%>% 
  filter(method== "observed")

Nodo66 <- X$iNextEst$`66`%>% 
  filter(method== "observed")

div <- bind_rows(Nodo57, Nodo58, Nodo60, Nodo61, Nodo62, Nodo63, Nodo64,
                 Nodo65, Nodo66) %>% 
  mutate(nodo= c(rep("Nodo57", 3),
                 rep("Nodo58", 3),
                 rep("Nodo60", 3),
                 rep("Nodo61", 3),
                 rep("Nodo62", 3),
                 rep("Nodo63", 3),
                 rep("Nodo64", 3),
                 rep("Nodo65", 3),
                 rep("Nodo66", 3))) %>% 
  select(nodo, order, qD, qD.LCL, qD.UCL) %>% 
  pivot_wider(names_from = order, 
              values_from = c(qD,qD.LCL, qD.UCL )) %>% 
  mutate(Integridad= c("0", "0", "0", "0", "1", "1", "1", "1", "1"))




# Gráficas cada qD --------------------------------------------------------

p1 <- ggplot(div, aes(x= nodo , y=qD_0, colour=nodo))+
  geom_pointrange(aes(shape=Integridad, ymin= qD.LCL_0, ymax=qD.UCL_0), size=.8, alpha=0.7)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x=NULL, y="Número de especies",
       title="Riqueza por nodo",
       subtitle = "")+
  theme(axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 16))+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.8)



p2 <- ggplot(div, aes(x= nodo , y=qD_1, colour=nodo))+
  geom_pointrange(aes(shape=Integridad, ymin= qD.LCL_1, ymax=qD.UCL_1), size=.8, alpha=0.7)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x=NULL, y="Número de especies igualmente abundantes",
       title="qD1",
       subtitle = "")+
  theme(axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 16))+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.8)

p3 <- ggplot(div, aes(x= nodo , y=qD_2, colour=nodo))+
  geom_pointrange(aes(shape=Integridad, ymin= qD.LCL_2, ymax=qD.UCL_2), size=.8, alpha=0.7)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x=NULL, y="Número de especies dominantes",
       title="qD2",
       subtitle = "")+
  theme(axis.text.x=element_text( size = 14),
        axis.text.y=element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        plot.title = element_text(size = 15),
        plot.subtitle = element_text(size = 16))+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.8)

plot_qs <- p1+p2+p3

ggsave(filename = "figures/Qs.png", plot = plot_qs, width = 20, height = 7)

# Gráfica perfiles por nodo -----------------------------------------------

perf <- bind_rows(Nodo57, Nodo58, Nodo60, Nodo61, Nodo62, Nodo63, Nodo64,
                 Nodo65, Nodo66) %>% 
  mutate(nodo= c(rep("Nodo57", 3),
                 rep("Nodo58", 3),
                 rep("Nodo60", 3),
                 rep("Nodo61", 3),
                 rep("Nodo62", 3),
                 rep("Nodo63", 3),
                 rep("Nodo64", 3),
                 rep("Nodo65", 3),
                 rep("Nodo66", 3))) %>% 
  select(nodo, order, qD, qD.LCL, qD.UCL) %>% 
  mutate(Integridad= rep(c("0", "0", "0", "0", "1", "1", "1", "1", "1"), each=3))


PerfPlot <- ggplot(perf, aes(x=nodo, y=qD, 
                             group= as.factor(order), 
                             color=as.factor(order)))+
            geom_pointrange(aes(shape=Integridad, ymin=qD.LCL, ymax= qD.UCL),
                  position = position_dodge(width = 0.5),
                  size=1,
                  alpha=0.7)+
            labs(x= "Nodo",color= " Diversidad")+
            theme_ipsum_tw(grid="Y", axis=T)+
            labs(x=NULL, y="Número de especies dominantes",
                title="Perfiles de diversidad por nodo",
                subtitle = "")+
            theme(axis.text.x=element_text( size = 18),
                  axis.text.y=element_text(size = 18),
                 axis.title.y = element_text(size = 18),
                 plot.title = element_text(size = 19),
                 plot.subtitle = element_text(size = 16),
                 legend.text = element_text(size = 16),
                 legend.title = element_text(size = 18),)+
            scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.8)

ggsave(filename = "figures/perfPlot.png", plot = PerfPlot, width = 20, height = 7)

# Modelos -----------------------------------------------------------------

adicionales <- data.frame(nodo=c("Nodo59"),
                          qD_0= rep(0,1),
                          qD_1= rep(0,1),
                          qD_2= rep(0,1))

div_mod <- div %>% 
  select(nodo, qD_0, qD_1, qD_2) %>% 
  bind_rows(adicionales) %>% 
  mutate(Integridad= c("0", "0", "0", "0", "1", "1", "1", "1", "1", "0"))


# Modelos qD0

summary(m0_qD0 <- glm(qD_0 ~ 1, data = div_mod, family = "poisson")) 
summary(m1_qD0 <- glm(qD_0 ~ Integridad, data = div_mod,  family = "poisson"))

cand_q0 <- list(m0_qD0, m1_qD0)
names_q0 <- c("Nulo",
              "qD0 ~ Integridad")
AICc_q0 <- aictab(cand_q0,
                  modnames = names_q0,
                  second.ord = T,
                  sort = T)

q0_predict <- ggpredict(m1_qD0, terms = "Integridad")

par_m0qD0 <- model_parameters(m0_qD0)
par_m1qD0 <- model_parameters(m1_qD0)

write.csv(par_m0qD0, file = "results/par_m0qD0.csv")
write.csv(par_m1qD0, file = "results/par_m1qD0.csv")

  # Plot modelo qD0 

qD0Plot <- ggplot(q0_predict, aes(x= as.factor(x), y=predicted, colour=x))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width= 0.11,
                size=0.8)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x= NULL, y="qD0",
       title="Diversidad estimada por categoria de integridad",
       subtitle = "")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none")+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.6)
qD0Plot 

# Modelos qD1

summary(m0_qD1 <- glm(qD_1 ~ 1, data = div_mod, family = "gaussian")) 
summary(m1_qD1 <- lm(qD_1 ~ Integridad, data = div_mod))

cand_q1 <- list(m0_qD1, m1_qD1)
names_q1 <- c("Nulo",
              "qD1 ~ Integridad")
AICc_q1 <- aictab(cand_q1,
                  modnames = names_q1,
                  second.ord = T,
                  sort = T)

q1_predict <- ggpredict(m1_qD1, terms = "Integridad")

par_m0qD1 <- model_parameters(m0_qD1)
par_m1qD1 <- model_parameters(m1_qD1)

write.csv(par_m0qD1, file = "results/par_m0qD1.csv")
write.csv(par_m1qD1, file = "results/par_m1qD1.csv")

## Plot modelo qD1

qD1Plot <- ggplot(q1_predict, aes(x= as.factor(x), y=predicted, colour=x))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width= 0.11,
                size=0.8)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x= "Categoria de integridad", y="qD1")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        axis.title.x = element_text(size = 19),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none")+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.6)

# Modelos qD2

summary(m0_qD2 <- glm(qD_2 ~ 1, data = div_mod, family = "gaussian")) 
summary(m1_qD2 <- lm(qD_2 ~ Integridad, data = div_mod))

cand_q2 <- list(m0_qD2, m1_qD2)
names_q2 <- c("Nulo",
              "qD2 ~ Integridad")
AICc_q2 <- aictab(cand_q2,
                  modnames = names_q2,
                  second.ord = T,
                  sort = T)

q2_predict <- ggpredict(m1_qD2, terms = "Integridad")

par_m0qD2 <- model_parameters(m0_qD2)
par_m1qD2 <- model_parameters(m1_qD2)

write.csv(par_m0qD2, file = "results/par_m0qD2.csv")
write.csv(par_m1qD2, file = "results/par_m1qD2.csv")

## Plot modelos qD2

qD2Plot <- ggplot(q2_predict, aes(x= as.factor(x), y=predicted, colour=x))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width= 0.11,
                size=0.8)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x= NULL, y="qD2")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = "none")+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.6)

## Plot completo modelos diversidad especies

modQsPlot <- qD0Plot + qD1Plot + qD2Plot
ggsave(filename = "figures/modQsPlot.png", plot = modQsPlot, width = 20, height = 7)

# Goftest -----------------------------------------------------------------

testDispersion(m1_qD0)
check_overdispersion(m1_qD0)


check_model(m1_qD1, check = c("qq","linearity", "normality"))
check_normality(m1_qD1)
plot(check_normality(m1_qD1), type= "qq")

windows()
check_model(m1_qD2)
check_normality(m1_qD2)
plot(check_normality(m1_qD2), type= "qq")


# Biomasa -----------------------------------------------------------------
str(PM)

biomass <- PM %>% 
  filter(Edad %in% c("Ad", "SAd")) %>% 
  drop_na(Epiteto) %>% 
  filter(!grepl('Recaptura', Obs)) %>% 
  group_by(Nodo, Epiteto) %>% 
  summarise(m_bio= mean(masa, na.rm=T)) %>% 
  mutate(Especies= case_when(
    Epiteto== "irroratus"~ "L. irroratus",
    Epiteto== "hispidus" ~ "S. hispidus",
    Epiteto== "fulvescens" ~ "R. fulvescens",
    Epiteto== "morfo 1" ~ "P. morfo 1",
    Epiteto== "morfo 2" ~ "P. morfo 2",
    Epiteto== "morfo 3" ~ "P. morfo 3"))


#%>% 
biomass_div <- PM %>% 
    filter(Edad %in% c("Ad", "SAd")) %>% 
    drop_na(Epiteto) %>% 
    filter(!grepl('Recaptura', Obs)) %>% 
    group_by(Nodo, Epiteto) %>% 
    summarise(m_bio= mean(masa, na.rm=T)) %>%  
    pivot_wider(names_from = Nodo,
              values_from = m_bio) %>% 
    replace(is.na(.), 0) %>% 
    column_to_rownames("Epiteto")
    

# Visualizaci?n promedios
  
plot_biomass <- ggplot(biomass, aes(x=Especies, y=m_bio))+
  geom_point(aes(color=as.factor(Nodo)), size=2.5)+
  theme_ipsum_tw(axis=T, grid = "Y")+
  labs(x= NULL, y="Masa",
       title="Promedio biomasa por especie y nodo")+
  theme(axis.text.x=element_text(face="italic", size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 16))+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.2, end= 0.8)
plot_biomass 

#write.csv(biomass, "div_biomass.csv")

# iNEXT biomasa
# 
X_biom <- iNEXT(biomass_div, q=c(0,1,2), datatype="abundance")

Nodo58b <- X_biom$iNextEst$`58` %>% 
  filter(method== "observed")

Nodo61b <- X_biom $iNextEst$`61`%>% 
  filter(method== "observed")


Nodo62b <- X_biom $iNextEst$`62`%>% 
  filter(method== "observed")

Nodo63b <- X_biom $iNextEst$`63`%>% 
  filter(method== "observed")

Nodo64b <- X_biom $iNextEst$`64`%>% 
  filter(method== "observed")

Nodo65b <- X_biom $iNextEst$`65`%>% 
  filter(method== "observed")

Nodo66b <- X_biom $iNextEst$`66`%>% 
  filter(method== "observed")

divb <- bind_rows(Nodo58b, Nodo61b, Nodo62b, Nodo63b, Nodo64b,
                 Nodo65b, Nodo66b) %>% 
  mutate(nodo= c(rep("Nodo58", 3),
                 rep("Nodo61", 3),
                 rep("Nodo62", 3),
                 rep("Nodo63", 3),
                 rep("Nodo64", 3),
                 rep("Nodo65", 3),
                 rep("Nodo66", 3))) %>% 
  select(nodo, order, qD, qD.LCL, qD.UCL) %>% 
  pivot_wider(names_from = order, 
              values_from = c(qD,qD.LCL, qD.UCL ))

ggplot(divb, aes(x= nodo , y=qD_0))+
  geom_pointrange(aes(ymin= qD.LCL_0, ymax=qD.UCL_0))+
  theme_bw()

ggplot(divb, aes(x= nodo , y=qD_1))+
  geom_pointrange(aes(ymin= qD.LCL_1, ymax=qD.UCL_1))+
  theme_bw()

ggplot(divb, aes(x= nodo , y=qD_2))+
  geom_pointrange(aes(ymin= qD.LCL_2, ymax=qD.UCL_2))+
  theme_bw()

perfb <- bind_rows(Nodo58b, Nodo61b, Nodo62b, Nodo63b, Nodo64b,
                  Nodo65b, Nodo66b) %>% 
  mutate(nodo= c(rep("Nodo58", 3),
                 rep("Nodo61", 3),
                 rep("Nodo62", 3),
                 rep("Nodo63", 3),
                 rep("Nodo64", 3),
                 rep("Nodo65", 3),
                 rep("Nodo66", 3))) %>% 
  select(nodo, order, qD, qD.LCL, qD.UCL)

## Gráfico perfil diversidad biomasa

porfbPlot<- ggplot(perfb, aes(x=nodo, y=qD, 
                 group= as.factor(order), 
                 color=as.factor(order)))+
  geom_pointrange(aes(ymin=qD.LCL, ymax= qD.UCL),
                  position = position_dodge(width = 0.5),
                  size=1, alpha=0.7)+
                  theme_ipsum_tw(grid="Y", axis=T)+
                  labs(x=NULL, y="Biomasa",
                  title="Perfiles de diversidad/biomasa por nodo",
                  subtitle = "",
                  color = "Diversidad")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        plot.title = element_text(size = 19),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),)+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.8)

ggsave(filename = "figures/perfbPlot.png", plot = porfbPlot, width = 20, height = 7)




# Modelos biomasa ---------------------------------------------------------

div_modb <- divb %>% 
  select(nodo, qD_0, qD_1, qD_2) %>% 
  bind_rows(adicionales) %>% 
  mutate(Integridad= c("0", "0", "1", "1", "1", "1", "1", "0", "0", "0"))

# Modelos qD0

summary(m0_qD0b <- glm(qD_0 ~ 1, data = div_modb, family = "poisson")) 
summary(m1_qD0b <- glm(qD_0 ~ Integridad, data = div_modb,  family = "poisson"))

cand_q0b <- list(m0_qD0b, m1_qD0b)
names_q0b <- c("Nulo",
              "qD0_b ~ Integridad")
(AICc_q0b <- aictab(cand_q0b,
                  modnames = names_q0b,
                  second.ord = T,
                  sort = T))

q0b_predict <- ggpredict(m1_qD0b, terms = "Integridad")

par_m0qD0b <- model_parameters(m0_qD0b)
par_m1qD0b <- model_parameters(m1_qD0b)

write.csv(par_m0qD0b, file = "results/par_m0qD0b.csv")
write.csv(par_m1qD0b, file = "results/par_m1qD0b.csv")

## Plot qDo biomasa

qD0Bio <- ggplot(q0b_predict, aes(x= as.factor(x), y=predicted, colour=factor(x)))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width= 0.11,
                size=0.8)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x= NULL, y="qD0 Biomasa",
       title="Biomasa estimada por categoria de integridad",
       subtitle = "")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.position = "none")+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.6)



# Modelos qD1

summary(m0_qD1b <- lm(qD_1 ~ 1, data = div_modb)) 
summary(m1_qD1b <- lm(qD_1 ~ Integridad, data = div_modb))

cand_q1b <- list(m0_qD1b, m1_qD1b)
names_q1b <- c("Nulo",
              "qD1b ~ Integridad")
(AICc_q1b <- aictab(cand_q1b,
                  modnames = names_q1b,
                  second.ord = T,
                  sort = T))

q1b_predict <- ggpredict(m1_qD1b, terms = "Integridad")

par_m0qD1b <- model_parameters(m0_qD1b)
par_m1qD1b <- model_parameters(m1_qD1b)

write.csv(par_m0qD1b, file = "results/par_m0qD1b.csv")
write.csv(par_m1qD1b, file = "results/par_m1qD1b.csv")

## Plot modelo qD1 biomasa

qD1Bio <- ggplot(q1b_predict, aes(x= as.factor(x), y=predicted, colour=factor(x)))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width= 0.11,
                size=0.8)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x= "Categoria integridad", y="qD1 Biomasa")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        axis.title.x = element_text(size = 19),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none")+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.6) 

# Modelos qD2

summary(m0_qD2b <- lm(qD_2 ~ 1, data = div_modb)) 
summary(m1_qD2b <- lm(qD_2 ~ Integridad, data = div_modb))

cand_q2b <- list(m0_qD2b, m1_qD2b)
names_q2b <- c("Nulo",
              "qD2b ~ Integridad")
(AICc_q2b <- aictab(cand_q2b,
                  modnames = names_q2b,
                  second.ord = T,
                  sort = T))

q2b_predict <- ggpredict(m1_qD2b, terms = "Integridad")

par_m0qD2b <- model_parameters(m0_qD2b)
par_m1qD2b <- model_parameters(m1_qD2b)

write.csv(par_m0qD2b, file = "results/par_m0qD2b.csv")
write.csv(par_m1qD2b, file = "results/par_m1qD2b.csv")

## Plot modelo qD2 biomasa

qD2Bio <- ggplot(q2b_predict, aes(x= as.factor(x), y=predicted, colour=factor(x)))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin= conf.low, ymax= conf.high),
                width= 0.11,
                size=0.8)+
  theme_ipsum_tw(grid="Y", axis=T)+
  labs(x= NULL, y="qD2 Biomasa")+
  theme(axis.text.x=element_text( size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.y = element_text(size = 19),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = "none")+
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.3, end= 0.6)

modQsBioPlot <- qD0Bio + qD1Bio + qD2Bio
ggsave(filename = "figures/modQsBioPlot.png", plot = modQsBioPlot, width = 20, height = 7)

# Goftest -----------------------------------------------------------------

testDispersion(m1_qD0b)
check_overdispersion(m1_qD0b)

windows()
check_model(m1_qD1b)
check_normality(m1_qD1b)
plot(check_normality(m1_qD1b), type= "qq")


check_model(m1_qD2b)
check_normality(m1_qD2b)
plot(check_normality(m1_qD2b), type= "qq")

