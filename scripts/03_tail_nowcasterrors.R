#------------------------------------------------------------------------#
# Tail Nowcast Errors.                                                   #
#                                                                        #
# Michael Pfarrhofer, Maximilian Boeck                                   #
#------------------------------------------------------------------------#
rm(list=ls())

# use realtime and/or rolling-window factors
rt = TRUE
rw = TRUE
if(!rt & rw) stop("not possible")

ne_auxarima <- TRUE

library(forecast)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(zoo)
library(reshape2)

# --------------------------------------------------------------------------------------
# setup
set.p <- c(0.05,0.1,0.16,0.5,0.84,0.9,0.95)

load(paste0("./results/GDP_quantiles",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),"_qr_dhslags4.rda"))
load("./results/SPF_mean_quantiles.rda")
load("./results/SPF_quantiles.rda")

rownames(quant_mcmc_post) = rownames(Y)
rownames(nc_quant)        = as.character(as.Date.ts(nc_quant))

# actual median forecast
nc_quant[,"mean"] <- window(nc_quant_mean[,4],start=start(nc_quant),end=end(nc_quant))

# plot the quantiles over time
quant_df   = data.frame("time"=rownames(Y),"GDP"=Y,quant_mcmc_post)
ncquant_df = data.frame("time"=as.character(as.Date.ts(nc_quant)),nc_quant)

colnames(quant_df)[-c(1,2)] = paste0("GDP_p",set.p*100)
colnames(quant_df)[2]       = "GDP_t"
colnames(ncquant_df)[-1]    = c(paste0("GDP_p",set.p*100),"GDP_t")

quant_m   = melt(quant_df)
ncquant_m = melt(ncquant_df)

colnames(quant_m)[3]   = "ACT"
colnames(ncquant_m)[3] = "SPF"

full_df <- left_join(quant_m,ncquant_m,by=c("time","variable")) %>%
  mutate("NCEraw" = ACT - SPF)  %>%
  subset(time %in% as.character(seq(as.Date("1970-01-01"),as.Date("2019-12-01"),by="month")))
nce_df <- full_df %>% dplyr::select(time, variable, NCEraw) %>%
  pivot_wider(names_from = variable, values_from = NCEraw)

# --------------------------------------------------------------------------------------
# purge predictable component from the nowcast errors
nce_clean_df <- nce_df
if(ne_auxarima){
  for(i in 2:ncol(nce_clean_df)){
    tmp <- auto.arima(nce_df[,i])
    nce_clean_df[,i] <- resid(tmp)
  } 
}
full_df <- nce_clean_df %>% melt() %>%
  rename("NCE"="value") %>%
  right_join(full_df)
levels(full_df$variable) <- c("Actual",paste0("p=",set.p))

pp1 <- full_df %>%
  subset(variable %in% c("Actual","p=0.5")) %>%
  ggplot() +
  geom_line(aes(x=as.Date(time),y=NCE,color=variable,linewidth=variable)) +
  scale_color_manual(name="",values=c("black","navy")) +
  scale_linewidth_manual(name="",values=c(0.8,0.8)) +
  geom_hline(yintercept = 0) +
  theme_cowplot() + ylab("Nowcast error") + xlab("") +
  coord_cartesian(expand = FALSE,ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), 
        legend.position = "bottom", legend.key.width=unit(1,"cm"))

pp2 <- full_df %>%
  subset(variable %in% c("p=0.1","p=0.9")) %>%
  ggplot() +
  geom_line(aes(x=as.Date(time),y=NCE,color=variable,linewidth=variable)) +
  scale_color_manual(name="",values=c("firebrick","navy")) +
  scale_linewidth_manual(name="",values=c(0.8,0.8)) +
  geom_hline(yintercept = 0) +
  theme_cowplot() + ylab("") + xlab("") +
  coord_cartesian(expand = FALSE,ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), 
        legend.position = "bottom", legend.key.width=unit(1,"cm"))

pdf(file=paste0("./nowcast_errors",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"),width=10,height=4)
plot_grid(pp1,pp2,ncol=2)
dev.off()

# --------------------------------------------------------------------------------------
# plot quantiles over time
full_tmp <- full_df
levels(full_tmp$variable) <- gsub("p=0.","p",levels(full_tmp$variable))

ACT_wide <- full_tmp %>% 
  dplyr::select(time,variable,ACT) %>%
  pivot_wider(names_from = variable, values_from = ACT)
SPF_wide <- full_tmp %>% 
  dplyr::select(time,variable,SPF) %>%
  pivot_wider(names_from = variable, values_from = SPF)

levels(full_tmp$variable) <- c("Actual","p5","p10","p16","TVP-QR","p84","p90","p95")
pp_gdp <- ggplot(ACT_wide,aes(x=as.Date(time))) +
  geom_ribbon(aes(ymin=p1,ymax=p9),alpha=0.3) + 
  
  geom_line(aes(y = ACT, x=as.Date(time), linetype = variable),linewidth = 0.8,
            data = subset(full_tmp, variable %in% c("Actual","TVP-QR"))) +
  scale_linetype_manual(name = "",values = c("solid","dashed")) +
  
  geom_hline(yintercept=0) +
  xlab("") + ylab("GDP") +
  theme_cowplot() + coord_cartesian(expand=FALSE) +
  theme(axis.line.x = element_blank(), 
        legend.position = "bottom",
        legend.key.width = unit(1,"cm"))

levels(full_tmp$variable) <- c("Median","p5","p10","p16","Ensemble","p84","p90","p95")
pp_spf <- ggplot(SPF_wide,aes(x=as.Date(time))) +
  geom_ribbon(aes(ymin=p1,ymax=p9),alpha=0.3) + 
  
  geom_line(aes(y = SPF, x=as.Date(time), linetype = variable),linewidth = 0.8,
            data = subset(full_tmp, variable %in% c("Median","Ensemble"))) +
  scale_linetype_manual(name = "",values = c("solid","dashed")) +
  
  geom_hline(yintercept=0) +
  xlab("") + ylab("SPF nowcast") +
  theme_cowplot() + coord_cartesian(expand=FALSE) +
  theme(axis.line.x = element_blank(), 
        legend.position = "bottom",
        legend.key.width = unit(1,"cm"))

pdf(file=paste0("./gdp_quantiles",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"),width=10,height=4)
plot_grid(pp_gdp,pp_spf,ncol=1)
dev.off()

# --------------------------------------------------------------------------------------
# scatterplots
library(GGally)
nce_density <- full_df %>% 
  subset(variable %in% c("Actual","p=0.1","p=0.5","p=0.9")) %>%
  dplyr::select(time,variable,NCE) %>%
  pivot_wider(names_from = variable, values_from = NCE)
colnames(nce_density) <- c("time","GDP_t","GDP_p10","GDP_p50","GDP_p90")
nce_density <- nce_density %>% dplyr::select(GDP_t, GDP_p10, GDP_p50, GDP_p90)

cor.xy <- round(cor(nce_density$GDP_p50,nce_density$GDP_t),digits=2)
pp1 <- ggplot(nce_density,aes(x = GDP_t, y = GDP_p50)) +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  geom_point(shape=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_cowplot() + coord_equal(expand=TRUE,xlim=c(-6,6),ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  ylab("p=0.5") +
  xlab(bquote("Actual,"~rho==.(cor.xy)))

cor.xy <- round(cor(nce_density$GDP_p50,nce_density$GDP_p10),digits=2)
pp2 <- ggplot(nce_density,aes(x = GDP_p50, y = GDP_p10)) +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  geom_point(shape=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_cowplot() + coord_equal(expand=TRUE,xlim=c(-6,6),ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  ylab("p=0.1") +
  xlab(bquote("p=0.5,"~rho==.(cor.xy)))

cor.xy <- round(cor(nce_density$GDP_p90,nce_density$GDP_p50),digits=2)
pp3 <- ggplot(nce_density,aes(x = GDP_p50, y = GDP_p90)) +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  geom_point(shape=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_cowplot() + coord_equal(expand=TRUE,xlim=c(-6,6),ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  ylab("p=0.9") +
  xlab(bquote("p=0.5,"~rho==.(cor.xy)))

cor.xy <- round(cor(nce_density$GDP_t,nce_density$GDP_p10),digits=2)
pp4 <- ggplot(nce_density,aes(x = GDP_t, y = GDP_p10)) +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  geom_point(shape=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_cowplot() + coord_equal(expand=TRUE,xlim=c(-6,6),ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  ylab("p=0.1") + 
  xlab(bquote("Actual,"~rho==.(cor.xy)))

cor.xy <- round(cor(nce_density$GDP_t,nce_density$GDP_p90),digits=2)
pp5 <- ggplot(nce_density,aes(x = GDP_t, y = GDP_p90)) +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  geom_point(shape=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_cowplot() + coord_equal(expand=TRUE,xlim=c(-6,6),ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  ylab("p=0.9") + 
  xlab(bquote("Actual,"~rho==.(cor.xy)))

cor.xy <- round(cor(nce_density$GDP_p10,nce_density$GDP_p90),digits=2)
pp6 <- ggplot(nce_density,aes(x = GDP_p10, y = GDP_p90)) +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  geom_point(shape=4) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_cowplot() + coord_equal(expand=TRUE,xlim=c(-6,6),ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  ylab("p=0.9") + 
  xlab(bquote("p=0.1,"~rho==.(cor.xy)))

pdf(file=paste0("./nowcast_scatter",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"),width=10,height=4)
plot_grid(pp1,pp2,pp3,ncol=3)
dev.off()

pdf(file=paste0("./nowcast_scatter_extended",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"),width=10,height=8)
plot_grid(pp4,pp1,pp5,
          pp2,pp6,pp3,ncol=3)
dev.off()

nce_density_m <- melt(nce_density)
levels(nce_density_m$variable) <- c("Actual","p=0.1","p=0.5","p=0.9")
pp_dens <- ggplot(nce_density_m) +
  geom_density(aes(y = value,fill=variable,color=variable,linewidth=variable,alpha = variable),
               position = "identity") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  scale_linewidth_manual(values = c(1,0.3,0.3,0.3),name="") +
  scale_alpha_manual(values = c(0,0.1,0.1,0.1),name="") +
  scale_fill_manual(values = c("white","firebrick3","black","navy"),name="") +
  scale_color_manual(values = c("black","firebrick3","black","navy"),name="") +
  theme_cowplot() + coord_cartesian(expand=TRUE,ylim=c(-6,6)) +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(0.1,"cm")) +
  ylab("NE") + xlab("Density")

pdf(file=paste0("./nowcast_scatter_v2",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".pdf"),width=10,height=3)
plot_grid(pp1,pp2,pp6,pp_dens,ncol=4,align="hv",axis="trbl",rel_widths = c(0.2,0.2,0.2,0.3))
dev.off()

save(full_df,file=paste0("./results/nc_errors_DF",ifelse(rt,"_rt",""),ifelse(rw,"_rw",""),".rda"))


