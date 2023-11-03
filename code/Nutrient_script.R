# Load library packages
library(tidyverse)
library(Hmisc)
library(moderndive)
library(gridExtra)


# Load data
#setwd("C:/Users/nathan/OneDrive - University of Southampton/Documents")
#setwd("C:/Users/natha/OneDrive - University of Southampton/Documents")
setwd("C:/Users/ndh1n17/OneDrive - University of Southampton/Documents")
#setwd("C:/Users/nathu/OneDrive - University of Southampton/Documents")

dat <- as_tibble(read.csv("mucus_data.csv", header=T))
dat <- dat[dat$Nutrient!="TRP",]
nutrient <- unique(dat$Nutrient)

### ---------------------------------------- Dialysis ----------------------------------------------------------

dia <- dat[dat$Exp=="DIA",]

zero <- NA
for (i in unique(dia$Nutrient[dia$Nutrient!="Nitrite" & dia$Nutrient!="Nitrate" & dia$Nutrient!="Phosphate" & dia$Nutrient!="Ammonium"])) {
  dummy <- dia[dia$Nutrient==i,]
  zero_dum <- dummy[c(1:4),]
  zero <- rbind(zero,zero_dum)
}
zero <- zero[-c(1,5,9),]
zero[c(3,6,8)] <- 0
zero$ID <- gsub(".$","0", zero$ID)
dia <- rbind(dia, zero)

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(dia$Nutrient)) {
  #Subsampling data
  dummy <- dia[dia$Nutrient==i,]
  #Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=Conc, color=Replicate, group=Replicate)) +
                  geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("Concentration (M)"))
  list_plots[[p]] <- plot
  p <- p+1
}
#grid.arrange(grobs=list_plots[1:9])
#grid.arrange(grobs=list_plots[10:18])
#grid.arrange(grobs=list_plots[19:25])

dia$Cor <- dia$Conc
for (i in unique(dia$Nutrient[dia$Nutrient!="DOC" & dia$Nutrient!="TDN"])) {
  for (j in 0:4) {
    dia$Cor[dia$Time_point==j & dia$Nutrient==i] <-  dia$Conc[dia$Time_point==j & dia$Nutrient==i] - dia$Conc[dia$Time_point==j & dia$Nutrient==i & dia$Replicate=="C"]
  }
}

dia$Cum <- dia$Cor/2 # to get the moles per L for 1g of mucus-FD
for (i in unique(dia$Nutrient)) {
  for (j in 1:4) {
    dia$Cum[dia$Time_point==j & dia$Nutrient==i] <-  dia$Cum[dia$Time_point==j-1 & dia$Nutrient==i] + dia$Cor[dia$Time_point==j & dia$Nutrient==i]/2
  }
}

# list_plots <- list()
# p <- 1
# #Creating a loop
# for (i in unique(dia$Nutrient)) {
#   # Subsampling data
#   dummy <- dia[dia$Nutrient==i,]
#   # Calling an empty plot
#   plot <- print(ggplot(dummy, aes(x = Time, y=Cum, color=Replicate, group=Replicate)) +
#                   geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("Cumulative nutrient release (nM)")+
#                   theme_classic())
#   list_plots[[p]] <- plot
#   p <- p+1
# }
# grid.arrange(grobs=list_plots[1:6])
# grid.arrange(grobs=list_plots[7:12])
# grid.arrange(grobs=list_plots[13:18])
# grid.arrange(grobs=list_plots[19:25])

# ---------------- Rates ------------------

table <- as.data.frame(matrix(data = NA, nrow = length(nutrient), ncol=2))
colnames(table) <- c("Mean", "SD")
rownames(table) <- nutrient
for (i in unique(dia$Nutrient)) {
  table$Mean[rownames(table)==i] <- mean(dia$Cum[dia$Nutrient==i & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)])
  table$SD[rownames(table)==i] <- sd(dia$Cum[dia$Nutrient==i & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)])
}

# DOC/TDN control correction

# From Tinta et al. 2020: DOC & TDN leached from Control= 178.67 ?M & 15.19 ?M-> moles= 89.335 & 7.595 ?moles
table$Mean[rownames(table)=="DOC"] <- mean(dia$Cum[dia$Nutrient=="DOC" & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)]-89.335)
table$SD[rownames(table)=="DOC"] <- sd(dia$Cum[dia$Nutrient=="DOC" & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)]-89.335)
table$Mean[rownames(table)=="TDN"] <- mean(dia$Cum[dia$Nutrient=="TDN" & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)]-7.595)
table$SD[rownames(table)=="TDN"] <- sd(dia$Cum[dia$Nutrient=="TDN" & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)]-7.595)

# Total AA
dia.aa <- dia[105:512,]
dia.tot.aa <- aggregate(dia.aa$Cum, by=list(dia.aa$Replicate, dia.aa$Time),sum)
colnames(dia.tot.aa) <- c("Replicate","Time", "Cum")

ggplot(data = dia.tot.aa, aes(x=Time, y=(Cum/1000), color=Replicate, group=Replicate)) + geom_point() + geom_line() + 
  ggtitle("Dialysis - Total amino-acids") + xlab("Time (hours)") + ylab("Cumulative AA release (µmoles g DW-1")+theme_classic()

AA.mean <- mean(dia.tot.aa$Cum[(dia.tot.aa$Replicate==1 | dia.tot.aa$Replicate==3) & dia.tot.aa$Time==24])
AA.sd <- sd(dia.tot.aa$Cum[(dia.tot.aa$Replicate==1 | dia.tot.aa$Replicate==3) & dia.tot.aa$Time==24])

table["AA.tot",] <- c(AA.mean,AA.sd)

### DIN calculation

dia.din <- dia[dia$Nutrient=="Nitrite" | dia$Nutrient=="Nitrate" | dia$Nutrient=="Ammonium",]
dia.din <- aggregate(dia.din$Cum[dia.din$Nutrient=="Nitrite" | dia.din$Nutrient=="Nitrate" | dia.din$Nutrient=="Ammonium"], by=list(dia.din$Replicate, dia.din$Time),sum)
colnames(dia.din) <- c("Replicate","Time", "Cum")

ggplot(data = dia.din, aes(x=Time, y=Cum, color=Replicate, group=Replicate)) + geom_point() + geom_line()

din.mean <- mean(dia.din$Cum[(dia.din$Replicate==1 | dia.din$Replicate==3) & dia.din$Time==24])
din.sd <- sd(dia.din$Cum[(dia.din$Replicate==1 | dia.din$Replicate==3) & dia.din$Time==24])

table["DIN",] <- c(din.mean,din.sd)

### DON calculation

dia.don <- dia[dia$Nutrient=="TDN",c(5,6,7,10)]
dia.din$Nutrient <- "DIN"
dia.nit <- rbind(dia.din[dia.din$Replicate!="C",],dia.don)
dia.don <- aggregate(dia.nit$Cum, by=list(dia.nit$Replicate, dia.nit$Time),diff)
colnames(dia.don) <- c("Replicate","Time", "Cum")

ggplot(data = dia.don, aes(x=Time, y=Cum, color=Replicate, group=Replicate)) + geom_point() + geom_line()

# From Tinta et al. 2020: DOC & TDN leached from Control= 178.67 ?M & 15.19 ?M-> moles= 89.335 & 7.595 ?moles

don.mean <- mean((dia$Cum[dia$Nutrient=="TDN" & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)]-7.595)-
                   dia.din$Cum[(dia.din$Replicate==1 | dia.din$Replicate==3) & dia.din$Time==24])
don.sd <- sd((dia$Cum[dia$Nutrient=="TDN" & dia$Time_point==4 & (dia$Replicate==1 | dia$Replicate==3)]-7.595)-
               dia.din$Cum[(dia.din$Replicate==1 | dia.din$Replicate==3) & dia.din$Time==24])

table["DON",] <- c(don.mean,don.sd)

### ---------------------------------------- Leaching ----------------------------------------------------------

leach <- dat[dat$Exp=="LE1" | dat$Exp=="LE2",]
leach <- leach[-(193:211),] # Removing ASW data for AA

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(leach$Nutrient)) {
  # Subsampling data
  dummy <- leach[leach$Nutrient==i & leach$Exp=="LE1",]
  # Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=Conc, color=Replicate, group=Replicate)) +
                  geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("Concentration (µM)")) +
    theme_classic()
  list_plots[[p]] <- plot
  p <- p+1
}
# grid.arrange(grobs=list_plots[1:6])
# grid.arrange(grobs=list_plots[7:12])
# grid.arrange(grobs=list_plots[13:18])
# grid.arrange(grobs=list_plots[19:24])
# grid.arrange(grobs=list_plots[25:27])

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(leach$Nutrient)) {
  # Subsampling data
  dummy <- leach[leach$Nutrient==i & leach$Exp=="LE2",]
  # Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=Conc, color=Replicate, group=Replicate)) +
                  geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("Concentration (nM)")+
                  theme_classic())
  list_plots[[p]] <- plot
  p <- p+1
}
# grid.arrange(grobs=list_plots[1:6])
# grid.arrange(grobs=list_plots[7:12])
# grid.arrange(grobs=list_plots[13:18])
# grid.arrange(grobs=list_plots[19:24])
# grid.arrange(grobs=list_plots[25:27])

leach$Diff <- leach$Conc
for (i in unique(leach$Nutrient)) {
  for (j in 0:3) {
    leach$Diff[leach$Time_point==j & leach$Nutrient==i] <-  leach$Conc[leach$Time_point==j & leach$Nutrient==i] - leach$Conc[leach$Time_point==0 & leach$Nutrient==i]
  }
}

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(leach$Nutrient)) {
  # Subsampling data
  dummy <- leach[leach$Nutrient==i & leach$Exp=="LE2",]
  # Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=Diff, color=Replicate, group=Replicate)) +
                  geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("change in nutrient (?M)") +
                  theme_classic())
  list_plots[[p]] <- plot
  p <- p+1
}
# grid.arrange(grobs=list_plots[1:6])
# grid.arrange(grobs=list_plots[7:12])
# grid.arrange(grobs=list_plots[13:18])
# grid.arrange(grobs=list_plots[19:24])
# grid.arrange(grobs=list_plots[25:27])

# ---------------- Rates ------------------

leach$Mol <- leach$Diff*4 # to get the moles per L and per grams of mucus-FD
final <- leach[leach$Time==24,]
initial <- leach[leach$Time==1 & (leach$Nutrient=="POC" | leach$Nutrient=="PON" | leach$Nutrient=="DOC" | leach$Nutrient=="TDN"),]


# -------  plots of DOM and POM

#### LE1

# Filter data
dummy <- leach[leach$Exp=="LE1" & (leach$Nutrient=="DOC" | leach$Nutrient=="POC"),]

# Create plot
ggplot(dummy, aes(x = Time, y=Conc, color=Nutrient, group=interaction(Nutrient, Replicate))) +
  geom_point(aes(shape=Replicate), size=2.5) +  # Increase the size of the points
  geom_line(linewidth=1) +  # Increase the size of the lines
  scale_shape_manual(values=c(15,16,17)) +
  scale_color_manual(values=c("#619CFF", "#00BA38")) +  # Use green and orange to represent the nutrients
  ggtitle("DOC and POC") + 
  xlab("Time (hour)") + 
  ylab("Concentration (µM)") +
  theme_classic()

# Filter data
dummy <- leach[leach$Exp=="LE1" & (leach$Nutrient=="TDN" | leach$Nutrient=="PON"),]

# Create plot
ggplot(dummy, aes(x = Time, y=Conc, color=Nutrient, group=interaction(Nutrient, Replicate))) +
  geom_point(aes(shape=Replicate), size=2.5) +  # Increase the size of the points
  geom_line(linewidth=1) +  # Increase the size of the lines
  scale_shape_manual(values=c(15,16,17)) +
  scale_color_manual(values=c("#619CFF", "#00BA38")) +  # Use green and orange to represent the nutrients
  ggtitle("TDN and DON") + 
  xlab("Time (hour)") + 
  ylab("Concentration (µM)") +
  theme_classic()

#### LE2

# Filter data
dummy <- leach[leach$Exp=="LE2" & (leach$Nutrient=="DOC" | leach$Nutrient=="POC"),]

# Create plot
ggplot(dummy, aes(x = Time, y=Conc, color=Nutrient, group=interaction(Nutrient, Replicate))) +
  geom_point(aes(shape=Replicate), size=2.5) +  # Increase the size of the points
  geom_line(linewidth=1) +  # Increase the size of the lines
  scale_shape_manual(values=c(15,16,17)) +
  scale_color_manual(values=c("#619CFF", "#00BA38")) +  # Use green and orange to represent the nutrients
  ggtitle("DOC and POC") + 
  xlab("Time (hour)") + 
  ylab("Concentration (µM)") +
  theme_classic()

# Filter data
dummy <- leach[leach$Exp=="LE2" & (leach$Nutrient=="TDN" | leach$Nutrient=="PON"),]

# Create plot
ggplot(dummy, aes(x = Time, y=Conc, color=Nutrient, group=interaction(Nutrient, Replicate))) +
  geom_point(aes(shape=Replicate), size=2.5) +  # Increase the size of the points
  geom_line(linewidth=1) +  # Increase the size of the lines
  scale_shape_manual(values=c(15,16,17)) +
  scale_color_manual(values=c("#619CFF", "#00BA38")) +  # Use green and orange to represent the nutrients
  ggtitle("TDN and DON") + 
  xlab("Time (hour)") + 
  ylab("Concentration (µM)") +
  theme_classic()

#### ------------------- Replicates combined ------------------------

#### Means and SD ###

agg <- aggregate(leach$Conc,by=list(leach$Exp, leach$Nutrient, leach$Time), mean, na.rm=T)
colnames(agg) <- c("Exp", "Nutrient","Time","mean")
agg$SD <- aggregate(leach$Conc,by=list(leach$Exp, leach$Nutrient, leach$Time), sd, na.rm=T)$x

# Subsampling data
dummy <- agg[agg$Exp=="LE1" & (agg$Nutrient=="TDN" | agg$Nutrient=="PON"),]

# Define color palette
color_palette <- c("#F8766D", "#619CFF")

# Calling an empty plot
plot <- print(ggplot(dummy, aes(x = Time, y=mean, color=Nutrient)) +
                geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Nutrient), alpha=0.3, colour=NA) +
                geom_point(aes(shape=Nutrient), size=4) + geom_line(aes(y=mean),size=1.2) +
                scale_color_manual(values=color_palette) + 
                scale_fill_manual(values=color_palette) +  # Add this line to change ribbon colors
                ggtitle("TDN and PON") + xlab("Time (hour)") +
                ylab("Concentration (µM)") + theme_classic())

# Subsampling data
dummy <- agg[agg$Exp=="LE2" & (agg$Nutrient=="TDN" | agg$Nutrient=="PON"),]
# Calling an empty plot
plot <- print(ggplot(dummy, aes(x = Time, y=mean, color=Nutrient)) +
                geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Nutrient), alpha=0.3, colour=NA) +
                geom_point(aes(shape=Nutrient), size=2.7) + geom_line(aes(y=mean),size=1.2) + ggtitle("TDN and PON") + xlab("Time (hour)") +
                ylab("Concentration (µM)") + theme_classic())

table.LE1 <- as.data.frame(matrix(data = NA, nrow = length(nutrient), ncol=2))
colnames(table.LE1) <- c("Mean", "SD")
rownames(table.LE1) <- nutrient

for (i in unique(leach$Nutrient)) {
  table.LE1$Mean[rownames(table.LE1)==i] <- mean(leach$Mol[leach$Exp=="LE1" & leach$Nutrient==i & leach$Time_point==3 & (leach$Replicate=="A" | leach$Replicate=="B")])
  table.LE1$SD[rownames(table.LE1)==i] <- sd(leach$Mol[leach$Exp=="LE1" & leach$Nutrient==i & leach$Time_point==3 & (leach$Replicate=="A" | leach$Replicate=="B")])
}

table.LE2 <- as.data.frame(matrix(data = NA, nrow = length(unique(leach$Nutrient)), ncol=2))
colnames(table.LE2) <- c("Mean", "SD")
rownames(table.LE2) <- unique(leach$Nutrient)
for (i in unique(leach$Nutrient)) {
  table.LE2$Mean[rownames(table.LE2)==i] <- mean(leach$Mol[leach$Exp=="LE2" & leach$Nutrient==i & leach$Time_point==3 & (leach$Replicate=="A" | leach$Replicate=="B")])
  table.LE2$SD[rownames(table.LE2)==i] <- sd(leach$Mol[leach$Exp=="LE2" & leach$Nutrient==i & leach$Time_point==3 & (leach$Replicate=="A" | leach$Replicate=="B")])
}

table.leach <- cbind(table.LE1,table.LE2)

### LE1 and LE2 combined

table.LE <- as.data.frame(matrix(data = NA, nrow = length(nutrient), ncol=2))
colnames(table.LE) <- c("Mean", "SD")
rownames(table.LE) <- nutrient

for (i in unique(leach$Nutrient)) {
  table.LE$Mean[rownames(table.LE)==i] <- mean(leach$Mol[leach$Nutrient==i & leach$Time_point==3 & (leach$Replicate=="A" | leach$Replicate=="B")])
  table.LE$SD[rownames(table.LE)==i] <- sd(leach$Mol[leach$Nutrient==i & leach$Time_point==3 & (leach$Replicate=="A" | leach$Replicate=="B")])
}

write.csv(table.LE, "Leaching_rates_combined.csv", row.names = T)

### ---------------------------------------- Mucus ----------------------------------------------------------

muc <- dat[dat$Exp=="MUC" & is.na(dat$Treatment)==F,]
muc$Subrep <- as.factor(substr(muc$Replicate,2,2))

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(muc$Nutrient)) {
  # Subsampling data
  dummy <- muc[muc$Nutrient==i,]
  # Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=Conc, color=Treatment, group=Replicate)) +
                  geom_point(aes(shape=Subrep)) + geom_line() + ggtitle(i) + xlab("Time (hour)") +
                  ylab("Concentration (nM)")+ theme_classic())
  list_plots[[p]] <- plot
  p <- p+1
}
grid.arrange(grobs=list_plots[1:6])
grid.arrange(grobs=list_plots[7:12])
grid.arrange(grobs=list_plots[13:18])
grid.arrange(grobs=list_plots[19:24])
grid.arrange(grobs=list_plots[25:27])

# DON and C:N ratio

muc.din <- muc[muc$Nutrient=="Nitrite" | muc$Nutrient=="Nitrate" | muc$Nutrient=="Ammonium",]
agg_din <- aggregate(muc.din$Conc, by=list(muc.din$Time,muc.din$Replicate, muc.din$Treatment), sum)
colnames(agg_din) <- c("Time", "Replicate","Treatment", "Conc")
ggplot(agg_din, aes(x = Time, y=Conc, color=Treatment, group=Replicate)) +
  geom_point() + geom_line() + ggtitle("DIN") + xlab("Time (hour)") +
  ylab("Concentration (nM)")+ theme_classic()

muc.TDN <- muc[muc$Nutrient=="TDN",c(-1,-2,-3,-7,-9,-10)]

muc.DON <- merge(agg_din, muc.TDN, by=c("Time", "Replicate", "Treatment"))
colnames(muc.DON) <- c("Time", "Replicate", "Treatment", "DIN","TDN","DON")
muc.DON$DON <- muc.DON$TDN - muc.DON$DIN

muc.DOC <- muc[muc$Nutrient=="DOC",c(-1,-2,-3,-7,-9,-10)]
muc.CN <- merge(muc.DON, muc.DOC, by=c("Time", "Replicate", "Treatment"))
colnames(muc.CN) <- c("Time", "Replicate", "Treatment", "DIN","TDN","DON","DOC")

df <- pivot_longer(muc.CN, 4:8, names_to = "Nutrient", values_to = "Conc" )
ggplot(data = df, aes(x=Time, y=Conc, col= Nutrient, shape=Replicate, group=interaction(Nutrient, Replicate)))+
  geom_point() + geom_line()

muc.CN$C_N <- muc.CN$DOC/muc.CN$DON

agg_N <- aggregate(df$Conc,by=list(df$Nutrient, df$Time, df$Treatment), mean, na.rm=T)
colnames(agg_N) <- c("Nutrient","Time","Treatment","mean")
agg_N$SD <- aggregate(df$Conc,by=list(df$Nutrient, df$Time, df$Treatment), sd, na.rm=T)$x

ggplot(agg_N[agg_N$Nutrient=="C_N",], aes(x = Time, y=mean, color=Treatment)) +
  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Treatment), alpha=0.3, colour=NA) +
  geom_point(aes(shape=Treatment), size=2.7) + geom_line(aes(y=mean),size=1.2) + ggtitle("C:N ratio") + xlab("Time (hour)") +
  ylab("Ratio") + theme_classic()

ggplot(agg_N[agg_N$Nutrient=="DOC",], aes(x = Time, y=mean, color=Treatment)) +
  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Treatment), alpha=0.3, colour=NA) +
  geom_point(aes(shape=Treatment), size=2.7) + geom_line(aes(y=mean),size=1.2) + ggtitle("DOC") + xlab("Time (hour)") +
  ylab("µM of DOC") + theme_classic()

ggplot(agg_N[agg_N$Nutrient=="DON",], aes(x = Time, y=mean, color=Treatment)) +
  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Treatment), alpha=0.3, colour=NA) +
  geom_point(aes(shape=Treatment), size=2.7) + geom_line(aes(y=mean),size=1.2) + ggtitle("DON") + xlab("Time (hour)") +
  ylab("µM of DON") + theme_classic()

df$Muc <- ifelse(df$Treatment %in% c("A","B"), "M", "C")
agg_N_M <- aggregate(df$Conc,by=list(df$Nutrient, df$Time, df$Muc), mean, na.rm=T)
colnames(agg_N_M) <- c("Nutrient","Time","Muc","mean")
agg_N_M$SD <- aggregate(df$Conc,by=list(df$Nutrient, df$Time, df$Muc), sd, na.rm=T)$x

ggplot(agg_N_M[agg_N_M$Nutrient=="C_N",], aes(x = Time, y=mean, color=Muc)) +
  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Muc), alpha=0.3, colour=NA) +
  geom_point(aes(shape=Muc), size=2.7) + geom_line(aes(y=mean),size=1.2) + ggtitle("C:N ratio") + xlab("Time (hour)") +
  ylab("Ratio") + theme_classic()


# rate of GLY consumption

mean(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="GLY"])/15
sd(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="GLY"])/15

t.test(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15, 
       muc$Conc[muc$Treatment=="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15)

wilcox.test(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15, 
            muc$Conc[muc$Treatment=="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15)

# rate of TAU consumption

mean(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="TAU"])
sd(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="TAU"])

mean(muc$Conc[muc$Treatment!="C" & muc$Time==42 & muc$Nutrient=="TAU"])
sd(muc$Conc[muc$Treatment!="C" & muc$Time==42 & muc$Nutrient=="TAU"])

mean((muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="TAU"]-
        muc$Conc[muc$Treatment!="C" & muc$Time==42 & muc$Nutrient=="TAU"])/37)

sd((muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="TAU"]-
      muc$Conc[muc$Treatment!="C" & muc$Time==42 & muc$Nutrient=="TAU"])/37)

t.test(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="TAU"]-
         muc$Conc[muc$Treatment!="C" & muc$Time==42 & muc$Nutrient=="TAU"], 
       muc$Conc[muc$Treatment=="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15)

wilcox.test(muc$Conc[muc$Treatment!="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15, 
            muc$Conc[muc$Treatment=="C" & muc$Time==5 & muc$Nutrient=="GLY"]/15)

aa <- unique(muc$Nutrient)[9:27]
muc.aa <- muc[muc$Nutrient==aa,]
leach.tot.aa <- aggregate(leach.aa$Mol, by=list(leach.aa$Exp, leach.aa$Replicate, leach.aa$Time),sum)
colnames(leach.tot.aa) <- c("Exp","Replicate","Time", "Mol")

ggplot(data = leach.tot.aa, aes(x=Time, y=Mol, color=Replicate, group=Replicate)) + geom_point() + geom_line() + facet_wrap(~Exp)

#### ------------------- Replicates combined ------------------------

#### Means and SD ###

agg <- aggregate(muc$Conc,by=list(muc$Nutrient, muc$Time, muc$Treatment), mean, na.rm=T)
colnames(agg) <- c("Nutrient","Time","Treatment","mean")
agg$SD <- aggregate(muc$Conc,by=list(muc$Nutrient, muc$Time, muc$Treatment), sd, na.rm=T)$x

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(muc$Nutrient)) {
  # Subsampling data
  dummy <- agg[agg$Nutrient==i,]
  # Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=mean, color=Treatment)) +
                  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD, fill=Treatment), alpha=0.3, colour=NA) +
                  geom_point(aes(shape=Treatment), size=2.7) + geom_line(aes(y=mean),size=1.2) + ggtitle(i) + xlab("Time (hour)") +
                  ylab("Concentration (nM)") + theme_classic())
  list_plots[[p]] <- plot
  p <- p+1
}
grid.arrange(grobs=list_plots[1:6])
grid.arrange(grobs=list_plots[7:12])
grid.arrange(grobs=list_plots[13:18])
grid.arrange(grobs=list_plots[19:24])
grid.arrange(grobs=list_plots[25:27])

list_plots[1:2]

muc$Diff <- muc$Conc
for (i in unique(muc$Nutrient)) {
  for (j in unique(muc$Time_point)) {
    muc$Diff[muc$Time_point==j & muc$Nutrient==i] <-  muc$Conc[muc$Time_point==j & muc$Nutrient==i] - muc$Conc[muc$Time_point==0 & muc$Nutrient==i]
  }
}
# 
# list_plots <- list()
# p <- 1
# #Creating a loop
# for (i in unique(muc$Nutrient)) {
#   # Subsampling data
#   dummy <- muc[muc$Nutrient==i,]
#   # Calling an empty plot
#   plot <- print(ggplot(dummy, aes(x = Time, y=Diff, color=Treatment, group=Replicate)) +
#                   geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("change in nutrient (?M)"))
#   list_plots[[p]] <- plot
#   p <- p+1
# }
# grid.arrange(grobs=list_plots[1:6])
# grid.arrange(grobs=list_plots[7:12])
# grid.arrange(grobs=list_plots[13:18])
# grid.arrange(grobs=list_plots[19:24])
# grid.arrange(grobs=list_plots[25:27])

### Mean concentration at beginning and end for phosphate and ammonium

mean(muc$Conc[muc$Nutrient=="Phosphate" & muc$Time==0])
sd(muc$Conc[muc$Nutrient=="Phosphate" & muc$Time==0])

mean(muc$Conc[muc$Nutrient=="Ammonium" & muc$Time==0])
sd(muc$Conc[muc$Nutrient=="Ammonium" & muc$Time==0])

mean(muc$Conc[muc$Nutrient=="Phosphate" & muc$Time==42])
sd(muc$Conc[muc$Nutrient=="Phosphate" & muc$Time==42])

mean(muc$Conc[muc$Nutrient=="Ammonium" & muc$Time==42])
sd(muc$Conc[muc$Nutrient=="Ammonium" & muc$Time==42])


# ---------------- Rates ------------------

muc$Mol <- muc$Diff

table.muc <- as.data.frame(matrix(data = NA, nrow = length(nutrient), ncol=2))
colnames(table.muc) <- c("Mean", "SD")
rownames(table.muc) <- nutrient
for (i in unique(muc$Nutrient)) {
  table.muc$Mean[rownames(table.muc)==i] <- mean(muc$Mol[muc$Nutrient==i & muc$Time_point==5 & muc$Treatment!="c"])
  table.muc$SD[rownames(table.muc)==i] <- sd(muc$Mol[muc$Nutrient==i & muc$Time_point==5 & muc$Treatment!="c"])
}

### Total AA calculation

aa <- unique(muc$Nutrient)[9:27]
muc.aa <- muc[433:1458,]
sum.tot.aa <- aggregate(muc.aa$Mol, by=list(muc.aa$Replicate, muc.aa$Time),sum, na.rm=T)
colnames(sum.tot.aa) <- c("Replicate","Time", "Mol")
sum.tot.aa$Treatment <- as.factor(substr(sum.tot.aa$Replicate,1,1))
sum.tot.aa$Subrep <- as.factor(substr(sum.tot.aa$Replicate,2,2))

ggplot(data = sum.tot.aa, aes(x=Time, y=Mol/100, color=Treatment, group=Replicate)) +
  geom_point(aes(shape=Subrep)) + geom_line()+ ggtitle("Total AA") + xlab("Time (hour)") + ylab("Concentration (µM)")+
  theme_classic()

muc.tot.aa <- aggregate(sum.tot.aa$Mol, by=list(sum.tot.aa$Treatment, sum.tot.aa$Time),mean, na.rm=T)
colnames(muc.tot.aa) <- c("Treatment","Time", "Mol")
muc.tot.aa$SD <- aggregate(sum.tot.aa$Mol, by=list(sum.tot.aa$Treatment, sum.tot.aa$Time),sd, na.rm=T)$x

ggplot(data = muc.tot.aa, aes(x=Time, y=Mol, color=Treatment, group=Treatment)) + 
  geom_ribbon(aes(ymin=Mol-SD, ymax=Mol+SD, fill=Treatment), alpha=0.3, colour=NA) +
  geom_point(size=2) + geom_line(size=1.2) + ggtitle("Total AA") + xlab("Time (hour)") + ylab("Concentration (?M)")

ggplot(data=muc.aa, aes(x=Time, y=Mol, colour=Nutrient, fill=Nutrient)) + geom_bar(stat = "identity", position = "stack") + facet_wrap(~Treatment)

AA.mean.LE2 <- mean(leach.tot.aa$Mol[leach.tot.aa$Exp=="LE2" & leach.tot.aa$Replicate!="C" & leach.tot.aa$Time==24])
AA.sd.LE2 <- sd(leach.tot.aa$Mol[leach.tot.aa$Exp=="LE2" & leach.tot.aa$Replicate!="C" & leach.tot.aa$Time==24])

table.leach["AA.tot",] <- c(AA.mean.LE1,AA.sd.LE1,AA.mean.LE2,AA.sd.LE2)

# Calculate the rate on each time point for mucus experiment

muc$Time_diff <- muc$Time

for (j in 2:6) {
  muc$Time_diff[muc$Time==unique(muc$Time)[j]] <- muc$Time[muc$Time==unique(muc$Time)[j]] - muc$Time[muc$Time==unique(muc$Time)[j-1]]
}

muc$rate <- muc$Mol
for (j in 2:6) {
  muc$rate[muc$Time==unique(muc$Time)[j]] <- (muc$Mol[muc$Time==unique(muc$Time)[j]] - 
                                                muc$Mol[muc$Time==unique(muc$Time)[j-1]])/unique(muc$Time_diff[muc$Time==unique(muc$Time)[j]])
}

list_plots <- list()
p <- 1
#Creating a loop
for (i in unique(muc$Nutrient)) {
  # Subsampling data
  dummy <- muc[muc$Nutrient==i,]
  # Calling an empty plot
  plot <- print(ggplot(dummy, aes(x = Time, y=rate, color=Treatment, group=Replicate)) +
                  geom_point() + geom_line() + ggtitle(i) + xlab("Time (hour)") + ylab("Concentration (M)"))
  list_plots[[p]] <- plot
  p <- p+1
}
grid.arrange(grobs=list_plots[1:6])
grid.arrange(grobs=list_plots[7:12])
grid.arrange(grobs=list_plots[13:18])
grid.arrange(grobs=list_plots[19:24])
grid.arrange(grobs=list_plots[25:27])

mean(muc$rate[muc$Treatment!="C" & muc$Time==30 & muc$Nutrient=="Phosphate"])
sd(muc$rate[muc$Treatment!="C" & muc$Time==30 & muc$Nutrient=="Phosphate"])

mean(muc$rate[muc$Treatment=="C" & muc$Time==30 & muc$Nutrient=="Phosphate"])
sd(muc$rate[muc$Treatment=="C" & muc$Time==30 & muc$Nutrient=="Phosphate"])

t.test(muc$rate[muc$Treatment!="C" & muc$Time==30 & muc$Nutrient=="Phosphate"], 
       muc$rate[muc$Treatment=="C" & muc$Time==30 & muc$Nutrient=="Phosphate"])
wilcox.test(muc$rate[muc$Treatment!="C" & muc$Time==30 & muc$Nutrient=="Phosphate"], 
            muc$rate[muc$Treatment=="C" & muc$Time==30 & muc$Nutrient=="Phosphate"])

mean(muc$rate[muc$Time==30 & muc$Nutrient=="Ammonium"])
sd(muc$rate[muc$Time==30 & muc$Nutrient=="Ammonium"])
