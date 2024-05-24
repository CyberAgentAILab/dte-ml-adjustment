rm(list = ls())

##------------------------------------------------------------------------------
## Load libraries and source files
##------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(foreign)
library(xgboost)
library(fastglm)
library(ggpubr)
library(RColorBrewer)

# Call source files
source( "functions.R" )
##------------------------------------------------------------------------------
## Setup
##------------------------------------------------------------------------------
# Set seed
set.seed(12345)

# Color palette for figures
cb_colors = brewer.pal(n = 8, name = "Dark2") # discrete colorblind palette

##------------------------------------------------------------------------------
## Load data
##------------------------------------------------------------------------------
df = read.csv('./data/data_ferraroprice.csv') %>% 
  na.omit() %>%
  rename(D = W) %>% 
  filter(D==3 | D==4) %>%
  mutate( D = as.factor( 1 * (D==3) )) %>%    ## treatment effect dummy
  data.frame()

## select variables 
vec.y = as.numeric( as.matrix(df$Y) )
vec.d = as.numeric( as.matrix(df$D) )
mat.x = as.matrix(df[,3:ncol(df)])
num_obs = length(vec.y) ## sample size (for whole dataset)

##------------------------------------------------------------------------------
## Estimation setup
##------------------------------------------------------------------------------
# Locations for DTE
vec.loc = seq( min(df$Y), 200, by = 1) 
# h for PTE
h.pte   = 10  ## gap length for PTE
# Number of folds for cross-fitting
F = 10 

##------------------------------------------------------------------------------
## Regression-adjusted DTE and PTE
##------------------------------------------------------------------------------
# F-fold cross-fitting setup
set.seed(123)
folds = sample(c(1:F), nrow(mat.x), replace=TRUE) 

start_time = Sys.time()
res.ml.est = DTE.ML.estimation(vec.y, vec.d, mat.x, vec.loc, h.pte, "gradient_boosting", 1)
end_time = Sys.time()
print(paste("Time spent:", end_time-start_time))

##------------------------------------------------------------------------------
## Plot PTE (simple)
##------------------------------------------------------------------------------
y.max = max(max(res.ml.est$pte$cv975), max(res.ml.est$pte$cv975)) + 1e-5
y.min = min(min(res.ml.est$pte$cv025), min(res.ml.est$pte$cv025)) - 1e-5

ggplot(res.ml.est$pte, aes(vec.loc, pte) ) + 
  geom_bar( stat = "identity", color= cb_colors[1], fill=cb_colors[1]) +
  geom_errorbar(aes(ymin = cv025, 
                    ymax = cv975),
                width= 5,                    # Width of the error bars
                #position=position_dodge(.9)
  ) +
  ylim(y.min, y.max) + 
  geom_hline(yintercept=0, color="black", size=0.01, alpha = .7) +
  theme_bw() + 
  labs(x= "Water Consumption", y="PTE")  +
  scale_x_continuous(breaks = seq(0,200,by=10), limit=c(-5,205)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("./result/FP_PTE_simple.pdf", width=5, height =3)  ## save 


##------------------------------------------------------------------------------
## Plot PTE (adjusted)
##------------------------------------------------------------------------------
y.max = max(max(res.ml.est$pte.ra$cv975), max(res.ml.est$pte.ra$cv975)) + 1e-5
y.min = min(min(res.ml.est$pte.ra$cv025), min(res.ml.est$pte.ra$cv025)) - 1e-5

ggplot(res.ml.est$pte.ra, aes(vec.loc, pte) ) + 
  geom_bar( stat = "identity", color= cb_colors[3], fill=cb_colors[3]) +
  geom_errorbar(aes(ymin = cv025, 
                    ymax = cv975),
                width= 5,                    # Width of the error bars
                #position=position_dodge(.9)
  ) +
  ylim(y.min, y.max) + 
  geom_hline(yintercept=0, color="black", size=0.01, alpha = .7) +
  theme_bw() + 
  labs(x= "Water Consumption", y="PTE")  +
  scale_x_continuous(breaks = seq(0,200,by=10), limit=c(-5,205))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("./result/FP_PTE_adj.pdf", width=5, height =3) ## save

##------------------------------------------------------------------------------
## Plot DTE (simple)
##------------------------------------------------------------------------------
ggplot() + 
  geom_line( aes(vec.loc, res.ml.est$dte$cv025), color= cb_colors[1], linetype=2) +
  geom_line( aes(vec.loc, res.ml.est$dte$cv975), color= cb_colors[1], linetype=2) +
  geom_ribbon(aes(x    = vec.loc, 
                  ymin = res.ml.est$dte$cv025, 
                  ymax = res.ml.est$dte$cv975), 
              fill = cb_colors[1], alpha = .3) +
  geom_line( aes(vec.loc, res.ml.est$dte$dte), color = cb_colors[1]) +
  theme_bw() + 
  #xlim(0, 200) +
  scale_x_continuous(breaks = seq(0,200,by=10), limit=c(0,200)) +
  geom_hline(yintercept=0, color="black", size=0.01, alpha = .3) +
  labs(title = "", 
       x= "Water Consumption", y="DTE") +
  theme(text=element_text(size=17))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("./result/FP_DTE_simple.pdf", width=5, height =3)  ## save 

##------------------------------------------------------------------------------
## Plot DTE (adjusted)
##------------------------------------------------------------------------------
ggplot() + 
  geom_line( aes(vec.loc, res.ml.est$dte.ra$cv025_moment), color= cb_colors[3], linetype=2) +
  geom_line( aes(vec.loc, res.ml.est$dte.ra$cv975_moment), color= cb_colors[3], linetype=2) +
  geom_ribbon(aes(x = vec.loc, 
                  ymin = res.ml.est$dte.ra$cv025_moment, 
                  ymax = res.ml.est$dte.ra$cv975_moment), 
              fill = cb_colors[3], alpha = .4)+
  geom_line( aes(vec.loc, res.ml.est$dte.ra$dte), color = cb_colors[3]) +
  theme_bw() + 
  #xlim(0, 200) +
  scale_x_continuous(breaks = seq(0,200,by=10), limit=c(0,200)) +
  geom_hline(yintercept=0, color="black", size=0.01, alpha = .3) +
  labs(title = "", 
       x= "Water Consumption", y="DTE") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("./result/FP_DTE_adj.pdf", width=5, height =3)  ## save 

