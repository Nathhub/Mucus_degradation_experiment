# Load library packages
library(tidyverse)
library(Hmisc)
library(moderndive)
library(gridExtra)

# Load data
#setwd("C:/Users/nathan/OneDrive - University of Southampton/Documents")
#setwd("C:/Users/natha/OneDrive - University of Southampton/Documents")
setwd("C:/Users/ndh1n17/OneDrive - University of Southampton/Documents")

dat <- as_tibble(read.csv("micro_data.csv", header=T))
dat$Treatment <- as.factor(dat$Treatment)
dat$Replicate <- paste(dat$Treatment, dat$Replicate, sep = "")
dat$Replicate <- as.factor(dat$Replicate)
dat$Subrep <- substr(dat$Replicate, 2,2)

legend <- c("Bacterial abundance (Cell ml-1)", "BCP (µg C L-1 h-1)", "DOC concentration (µg L-1)", "TDN concentration (µg L-1)")
scale <- list(scales::scientific, waiver(), waiver(), waiver())
list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(dat$Data)) {
  #Subsampling data
  dummy <- dat[dat$Data==i,]
  #Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=Conc, color=Treatment, group=Replicate)) +
                 geom_point(aes(shape=Subrep)) + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab(legend[p])) + scale_y_continuous(labels = scale[[p]])
 list_plots[[p]] <- plot
p <- p+1
}
grid.arrange(grobs=list_plots[1:4])

## BA at T0

mean(dat$Conc[dat$Data=="BA" & dat$Time==0])
sd(dat$Conc[dat$Data=="BA" & dat$Time==0])

#### Means and SD ###

agg <- aggregate(dat$Conc,by=list(dat$Data, dat$Time, dat$Treatment), mean, na.rm=T)
colnames(agg) <- c("Data","Time","Treatment","mean")
agg$SD <- aggregate(dat$Conc,by=list(dat$Data, dat$Time, dat$Treatment), sd, na.rm=T)$x

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(agg$Data)) {
  #Subsampling data
  dummy <- agg[agg$Data==i,]
  #Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=mean, color=Treatment)) +
                  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Treatment), alpha=0.3, colour=NA) + 
                  geom_point(size=2) + geom_line(aes(y=mean),size=1.2) + ggtitle(i) + xlab("Time (hour)") +
                  ylab(legend[p]) + theme_classic() + scale_y_continuous(labels = scale[[p]]))
  list_plots[[p]] <- plot
  p <- p+1
}
grid.arrange(grobs=list_plots[1:4])

# Calculate the difference from T0

dat$Diff <- dat$Conc
for (i in unique(dat$Data)) {
  for (j in unique(dat$Time)) {
    dat$Diff[dat$Time==j & dat$Data==i] <- dat$Conc[dat$Time==j & dat$Data==i] - dat$Conc[dat$Time==0 & dat$Data==i]
  }
}

agg <- aggregate(dat$Diff,by=list(dat$Data, dat$Time, dat$Treatment), mean, na.rm=T)
colnames(agg) <- c("Data","Time","Treatment","mean")
agg$SD <- aggregate(dat$Diff,by=list(dat$Data, dat$Time, dat$Treatment), sd, na.rm=T)$x

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(agg$Data)) {
  #Subsampling data
  dummy <- agg[agg$Data==i,]
  #Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=mean, color=Treatment)) +
                  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Treatment), alpha=0.3, colour=NA) + 
                  geom_point(aes(shape=Treatment), size=2.7) + geom_line(aes(y=mean),size=1) + ggtitle(i) + xlab("Time (hour)") +
                  ylab(legend[p]) + theme_classic() + scale_y_continuous(labels = scale[[p]]))
  list_plots[[p]] <- plot
  p <- p+1
}
grid.arrange(grobs=list_plots[1:4])

mean(dat$Diff[dat$Data=="DOC" & dat$Time==30] - dat$Diff[dat$Data=="DOC" & dat$Time==20])
sd(dat$Diff[dat$Data=="DOC" & dat$Time==30] - dat$Diff[dat$Data=="DOC" & dat$Time==20])

# C:N ratio of the dissolved matter

ratio <- dat[dat$Data=="DOC",]
ratio$TDN <- dat$Conc[dat$Data=="TDN"]
ratio$Ratio <- ratio$Conc/ratio$TDN

ggplot(ratio, aes(x = Time, y=Ratio, color=Treatment, group=Replicate)) + geom_point() + geom_line() + ggtitle("C:N ratio") +
  xlab("Time (hour)") + ylab("C:N ratio") + theme_classic() 

# Calculate the rate on each time point

dat$Time_diff <- dat$Time

for (j in 2:6) {
  dat$Time_diff[dat$Time==unique(dat$Time)[j]] <- dat$Time[dat$Time==unique(dat$Time)[j]] - dat$Time[dat$Time==unique(dat$Time)[j-1]]
}

dat$rate <- dat$Diff
for (i in unique(dat$Data)) {
  for (j in 2:6) {
    dat$rate[dat$Time==unique(dat$Time)[j] & dat$Data==i] <- (dat$Diff[dat$Time==unique(dat$Time)[j] & dat$Data==i] - 
      dat$Diff[dat$Time==unique(dat$Time)[j-1] & dat$Data==i])/unique(dat$Time_diff[dat$Time==unique(dat$Time)[j]])
  }
}

agg <- aggregate(dat$rate,by=list(dat$Data, dat$Time, dat$Treatment), mean, na.rm=T)
colnames(agg) <- c("Data","Time","Treatment","mean")
agg$SD <- aggregate(dat$rate,by=list(dat$Data, dat$Time, dat$Treatment), sd, na.rm=T)$x

legend <- c("Bacterial growth (Cells ml-1 h-1)", "Change in BCP (?g C L-1 h-2)", "Change in DOC concentration (?g C L-1 h-1)", "Change in TDN concentration (?g C L-1 h-1)")
list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(agg$Data)) {
  #Subsampling data
  dummy <- agg[agg$Data==i,]
  #Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=mean, color=Treatment)) + geom_bar(aes(fill=Treatment), stat="identity", position = "dodge") +
                  ggtitle(i) + xlab("Time (hour)") +
                  ylab(legend[p]) + theme_classic() + scale_y_continuous(labels = scale[[p]]))
  list_plots[[p]] <- plot
  p <- p+1
}
grid.arrange(grobs=list_plots[1:4])

### calculate the microbial growth rate constant

growth_rate <- dat[dat$Data=="BA",]
growth_rate$Log_N <- log10(growth_rate$Conc)
growth_rate$Time_point <- rep(1:6, each =9)
growth_rate$Gro_r <- NA

for (i in 2:6) {
  growth_rate$Gro_r[growth_rate$Time_point==i] <- ((growth_rate$Log_N[growth_rate$Time_point==i] - growth_rate$Log_N[growth_rate$Time_point==i-1])*2.303)/growth_rate$Time_diff[growth_rate$Time_point==i]
}

# average contents of carbon and nitrogen for coastal bacterial assemblages: 30.2 +/- 12.3 fg of C cell-1 and 5.8 +/- 1.5 fg of N cell-1 (Fukuda et al. 1998)
bacbio <- dat[dat$Data=="BA",]
bacbio$Biomass <- (bacbio$Conc * 30.2)/1000000 # ?g C L-1
bacbio$BM_diff <- (bacbio$Diff * 30.2)/1000000 # ?g C L-1
bacbio$DOC <- dat$Conc[dat$Data=="DOC"]
bacbio$DOC_diff <- dat$Diff[dat$Data=="DOC"]
bacbio$BCP <- dat$Conc[dat$Data=="BCP"]
bacbio$BR <- bacbio$DOC_diff-bacbio$BCP
bacbio$DOC_net <- bacbio$DOC - bacbio$Biomass
bacbio$DOC_rate <- dat$rate[dat$Data=="DOC"]
bacbio$TDN_rate <- dat$rate[dat$Data=="TDN"]
bacbio$TDN_diff <- dat$Diff[dat$Data=="TDN"]

ggplot(bacbio, aes(x = Time, y=DOC_diff, color=Treatment, group=Replicate)) + geom_point() + geom_line() + ggtitle("Bacterial Respiration") +
  xlab("Time (hour)") + ylab("Bacterial Respiration ?g C L-1 ") + theme_classic() 
# + scale_y_continuous(labels = scales::scientific)

mean(bacbio$Biomass[bacbio$Time==42 & bacbio$Treatment!="C"])
sd(bacbio$Biomass[bacbio$Time==42 & bacbio$Treatment!="C"])

