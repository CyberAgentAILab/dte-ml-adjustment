rm(list = ls())

##------------------------------------------------------------------------------
## Load libraries and source files
##------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
require(gridExtra)
library(ggpubr)
library(RColorBrewer)

cb_colors = brewer.pal(n = 8, name = "Dark2") # discrete colorblind palette

##------------------------------------------------------------------------------
##  Setup
##------------------------------------------------------------------------------
vec.n = c(500, 1000, 5000)
v = c("Simple", "OLS", "Lasso")   
quantiles = rep(seq(0.1, 0.9, 0.1), 3)

##------------------------------------------------------------------------------
## Simulation for DTE
##------------------------------------------------------------------------------
df.bias = c()
df.ratio = c()
df.rmse = c()

# True DTE
dte.true = rep(readRDS('./result/dte.true.rds'),3)

for(j in 1:length(vec.n)) { ## loop: sample size 
  n = vec.n[j] 
  file.info = paste0("dte_n", n)
  
  df.bias.temp = read.csv(paste0("./result/", file.info, "/", file.info, "_bias.csv"))
  df.bias = rbind(df.bias, df.bias.temp)
  
  df.rmse.temp = read.csv(paste0("./result/", file.info, "/", file.info, "_rmse.csv"))
  df.rmse = rbind(df.rmse, df.rmse.temp)
  
  df.ratio.temp = read.csv(paste0("./result/", file.info, "/", file.info, "_ratio.csv"))
  df.ratio = rbind(df.ratio, df.ratio.temp)
  
}

df.bias = cbind(df.bias, dte.true) %>%
  mutate(quantiles = rep(seq(0.1, 0.9, 0.1), 3)) %>%
  select( c(v, dte.true, n,  quantiles)) %>%
  mutate(Simple = 100*Simple/dte.true, OLS = 100*OLS/dte.true, Lasso = 100*Lasso/dte.true) %>%
  pivot_longer(cols = v,
               names_to = "est", 
               values_to = "val") %>%
  mutate(est = as.factor(est)) %>%   
  data.frame() 
df.bias = transform(df.bias, est= factor(est, levels = v))

df.ratio = df.ratio %>%
  mutate(quantiles = rep(seq(0.1, 0.9, 0.1), 3)) %>%
  select(c("OLS", "Lasso", n, quantiles)) %>%
  pivot_longer(cols = c("OLS", "Lasso"),
               names_to = "est", 
               values_to = "val") %>%
  mutate(est = as.factor(est)) %>%   
  data.frame() 
df.ratio = transform(df.ratio, est= factor(est, levels = c("OLS", "Lasso")))

##------------------------------------------------------------------------------
## DTE Bias figure 
##------------------------------------------------------------------------------
ggplot(df.bias, aes(x = quantiles, y = val, group = est, color = est, shape=est, linetype=est)) +
  geom_line(size=0.5) +
  geom_hline(yintercept = 0, size=0.3)+
  geom_point(size=2)+
  facet_grid(cols=vars(n),scales = "free_y", labeller = labeller(n = c("500" = "n=500", "1000" = "n=1000", "5000" = "n=5000")))+
  labs(title = " ",
       x = "Quantiles",
       y = "Bias Ratio (%)",
       color= "Estimator",
       shape = "Estimator",
       linetype="Estimator") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0.1,0.9,by=0.1)) +
  theme(panel.spacing = unit(1, "cm"))+
  scale_color_brewer(palette="Dark2") +
  scale_color_manual(values = c("Simple" = cb_colors[1], "OLS" = cb_colors[2], "Lasso" = cb_colors[3]),
                     labels = c("Simple", "Linear adjustment", "ML adjustment")) +
  scale_shape_manual(values = c("Simple" = 16, "OLS" = 17, "Lasso" = 15),
                     labels = c("Simple", "Linear adjustment", "ML adjustment"))+
  scale_linetype_manual(values = c("Simple" = "dashed", "OLS" = "dashed", "Lasso" = "solid"),
                        labels = c("Simple", "Linear adjustment", "ML adjustment"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face="bold"),
        plot.title = element_text(size = 16, face="bold"),
        legend.text =element_text(size = 14, face="bold"),
        legend.title =element_text(size = 16, face="bold"),
        strip.text = element_text(size = 14, face = "bold"))+
  ylim(-3,50)
ggsave("fig_dte_bias.pdf", width = 11, height=3.5)  ## save 

##------------------------------------------------------------------------------
## DTE RMSE reduction figure
##------------------------------------------------------------------------------
ggplot(df.ratio, aes(x = quantiles, y = val, group = est, color = est, shape=est, linetype=est)) +
  geom_line(size=0.5) +
  geom_hline(yintercept = 0, size=0.3)+
  geom_point(size = 2)+
  facet_grid(cols=vars(n), scales = "free_y", labeller = labeller(n = c("500"= "n=500", "1000" = "n=1000", "5000" = "n=5000"))) +
  labs(title = " ",
       x = "Quantiles",
       y = "RMSE Reduction (%)",
       color= "Estimator",
       shape= "Estimator",
       linetype= "Estimator") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0.10,0.90,by=0.10)) +
  theme(panel.spacing = unit(1, "cm"))+
  scale_color_brewer(palette="Dark2")+
  scale_color_manual(values = c("OLS" = cb_colors[2], "Lasso" = cb_colors[3]),
                     labels = c("Linear adjustment", "ML adjustment")) +
  scale_shape_manual(values = c("OLS" = 17, "Lasso" = 15),
                     labels = c("Linear adjustment", "ML adjustment"))+
  scale_linetype_manual(values = c("OLS" = "dashed", "Lasso" = "solid"),
                        labels = c("Linear adjustment", "ML adjustment"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face="bold"),
        plot.title = element_text(size = 16, face="bold"),
        legend.text =element_text(size = 14, face="bold"),
        legend.title =element_text(size = 16, face="bold"),
        strip.text = element_text(size = 14, face = "bold")) 
ggsave("fig_dte_ratio.pdf", width = 11, height=3.5)  ## save 


##------------------------------------------------------------------------------
## Simulation for QTE
##------------------------------------------------------------------------------
df.bias = c()
df.ratio = c()
df.rmse = c()


for(j in 1:length(vec.n)) { ## loop: sample size 
  n = vec.n[j] 
  file.info = paste0("qte_n", n)
  df.bias.temp = read.csv(paste0("./result/", file.info, "/", file.info, "_bias.csv"))
  df.bias = rbind(df.bias, df.bias.temp)
  
  df.rmse.temp = read.csv(paste0("./result/", file.info, "/", file.info, "_rmse.csv"))
  df.rmse = rbind(df.rmse, df.rmse.temp)
  
  df.ratio.temp = read.csv(paste0("./result/", file.info, "/", file.info, "_ratio.csv"))
  df.ratio = rbind(df.ratio, df.ratio.temp)
  
}

df.bias = df.bias %>%
  mutate(quantiles = rep(seq(0.10, 0.90, 0.10), length(vec.n))) %>%
  select( c(v, n,  quantiles)) %>%
  pivot_longer(cols = v,
               names_to = "est", 
               values_to = "val") %>%
  mutate(est = as.factor(est)) %>%   
  data.frame() 
df.bias = transform(df.bias, est= factor(est, levels = v))

df.ratio = df.ratio %>%
  mutate(quantiles = rep(seq(0.1, 0.9, 0.1), 3)) %>%
  select(c("OLS", "Lasso", n, quantiles)) %>%
  pivot_longer(cols = c("OLS", "Lasso"),
               names_to = "est", 
               values_to = "val") %>%
  mutate(est = as.factor(est)) %>%   
  data.frame() 
df.ratio = transform(df.ratio, est= factor(est, levels = c("OLS", "Lasso")))

##------------------------------------------------------------------------------
## QTE Bias figure 
##------------------------------------------------------------------------------
ggplot(df.bias, aes(x = quantiles, y = 100*val, group = est, color = est, shape=est, linetype=est)) +
  geom_line(size=0.5) +
  geom_hline(yintercept = 0, size=0.3)+
  geom_point(size=2)+
  facet_grid(cols=vars(n),scales = "free_y", labeller = labeller(n = c("500" = "n=500", "1000" = "n=1000", "5000" = "n=5000")))+
  labs(title = " ",
       x = "Quantiles",
       y = "Bias Ratio (%)",
       color= "Estimator",
       shape = "Estimator",
       linetype="Estimator") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0.1,0.9,by=0.1)) +
  theme(panel.spacing = unit(1, "cm"))+
  scale_color_brewer(palette="Dark2") +
  scale_color_manual(values = c("Simple" = cb_colors[1], "OLS" = cb_colors[2], "Lasso" = cb_colors[3]),
                     labels = c("Simple", "Linear adjustment", "ML adjustment")) +
  scale_shape_manual(values = c("Simple" = 16, "OLS" = 17, "Lasso" = 15),
                     labels = c("Simple", "Linear adjustment", "ML adjustment"))+
  scale_linetype_manual(values = c("Simple" = "dashed", "OLS" = "dashed", "Lasso" = "solid"),
                        labels = c("Simple", "Linear adjustment", "ML adjustment"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face="bold"),
        plot.title = element_text(size = 16, face="bold"),
        legend.text =element_text(size = 14, face="bold"),
        legend.title =element_text(size = 16, face="bold"),
        strip.text = element_text(size = 14, face = "bold"))+
  ylim(-4,52)
ggsave("fig_qte_bias.pdf", width = 11, height=3.5)  ## save 

##------------------------------------------------------------------------------
## QTE RMSE Reduction figure 
##------------------------------------------------------------------------------
ggplot(df.ratio, aes(x = quantiles, y = val, group = est, color = est, shape=est, linetype=est)) +
  geom_line(size=0.5) +
  geom_hline(yintercept = 0, size=0.3)+
  geom_point(size = 2)+
  facet_grid(cols=vars(n), scales = "free_y", labeller = labeller(n = c("500"= "n=500", "1000" = "n=1000", "5000" = "n=5000"))) +
  labs(title = " ",
       x = "Quantiles",
       y = "RMSE Reduction (%)",
       color= "Estimator",
       shape= "Estimator",
       linetype= "Estimator") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0.10,0.90,by=0.10)) +
  theme(panel.spacing = unit(1, "cm"))+
  scale_color_brewer(palette="Dark2")+
  scale_color_manual(values = c("OLS" = cb_colors[2], "Lasso" = cb_colors[3]),
                     labels = c("Linear adjustment", "ML adjustment")) +
  scale_shape_manual(values = c("OLS" = 17, "Lasso" = 15),
                     labels = c("Linear adjustment", "ML adjustment"))+
  scale_linetype_manual(values = c("OLS" = "dashed", "Lasso" = "solid"),
                        labels = c("Linear adjustment", "ML adjustment"))+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face="bold"),
        plot.title = element_text(size = 16, face="bold"),
        legend.text =element_text(size = 14, face="bold"),
        legend.title =element_text(size = 16, face="bold"),
        strip.text = element_text(size = 14, face = "bold"))+
  ylim(-4,52) 
ggsave("fig_qte_ratio.pdf", width = 11, height=3.5)  ## save 


##------------------------------------------------------------------------------
## Simulation for sequence of DGPs
##------------------------------------------------------------------------------
ratio_results= matrix(NA, 9, 10)

## Results from s=1,...,10
for (dgp_seq in 1:10){
  dgp_number = dgp_seq
  dgp = read.csv( paste0("./result/dgp_seq", dgp_number, "_ratio.csv")) %>%
    mutate(quantiles = seq(0.1, 0.9, 0.1)) %>%
    select(c( n, quantiles, "Lasso")) 
  
  ratio_results[,dgp_number] = round(dgp$Lasso,2)
}

index = 1:10
data = expand.grid(Quantile = seq(0.1,0.9, by=0.1), Index = factor(index))
data$Value = as.vector(ratio_results)

## Plot
data_summary = data %>% filter(Quantile==0.4) 

ggplot(data, aes(x = Quantile, y = Value, group=Index, color= Index)) +
  geom_line(linetype="dotdash") + 
  geom_point(size=3, shape=17) +
  geom_label(data = data_summary, aes(label = paste0("s=",Index)), hjust =-0.2, vjust = -0.2, size = 4, alpha = 0.7) +  # Add labels inside text boxes
  theme_minimal() +
  labs(title = " ",
       x = "Quantiles",
       y = "RMSE Reduction (%)") +
  scale_x_continuous(breaks = seq(0.10,0.90,by=0.10)) +
  guides(color = FALSE) +
  annotate("label", x = max(data$Quantile)-0.1, y = max(data$Value), label = "s=1: more relevant \ns=10: less relevant", color = "black", size = 4) +  # Add text at top left
  theme(
    axis.text.x = element_text(angle = 0),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face="bold"),
    plot.title = element_text(size = 16, face="bold"),
    legend.text =element_text(size = 14, face="bold"),
    legend.title =element_text(size = 16, face="bold"),
    strip.text = element_text(size = 14, face = "bold")) 

ggsave("ratio_relevance.pdf", width = 11, height=8)  ## save 


  
  