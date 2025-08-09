## Fry Size
# last run 2025 08 09
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# libraries
library(lme4) # version 1.1-37
library(lmerTest) # version 3.1-3
library(ggplot2) # version 3.5.2
library(rmcorr) # version 0.7.0


# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# file <- file.choose()
file <- "input files/Fry Size 2017.csv"

d <- read.csv(file)

d$Egg.Date <- as.numeric(format(as.Date(as.character(d$Egg.Date), "%m/%d/%y"), "%j"))

x1 <- lmer(Mean5.Weight ~ Trt + Fry.Count + Week + 
	Trt:Week + (1|Tank), data = d)
summary(x1)
anova(x1, type = "II")
ranova(x1)
hist(resid(x1))
shapiro.test(resid(x1))


x2 <- lmer(Mean5.Length ~ Trt + Fry.Count + Week + 
             Trt:Week + (1|Tank), data = d)
summary(x2)
anova(x2, type = "II")
ranova(x2)
hist(resid(x2))
shapiro.test(resid(x2))


colors <- c("#53B7E0", "#E0A473")
dg <- 0.2

fig1 <- ggplot(data = d[!(is.na(d$Mean5.Length)), ], aes(x = Week, y = Mean5.Length, 
	color = as.factor(Trt))) 
fig1 <- fig1 +
	geom_point(aes(x = Week, y = Mean5.Length, group = as.factor(Tank), 
		color = as.factor(Trt)), size = 2, shape = 19, alpha = 0.25, position = position_dodge(width = dg)) +
	geom_line(aes(x = Week, y = Mean5.Length, group = as.factor(Tank), 
		color = as.factor(Trt)), alpha = 0.25, position = position_dodge(width = dg)) +
  scale_color_discrete(type = colors, name = "Care Treatment") +
  stat_summary(fun.y = mean, fun.ymin = meanminus, fun.ymax = meanplus, 
               geom = "pointrange", size = 1.5, lwd = 1.5, position = position_dodge(width = dg)) +
  stat_summary(fun = mean, geom = "line", size = 1.25, 
               position = position_dodge(width = dg)) +
  ylab("Mean Standard Length of 5 Fish (mm)") + xlab("Week Since Fertilization") + theme_bw()
	
fig1


dg <- 0.2
fig2 <- ggplot(data = d[!(is.na(d$Mean5.Weight)), ], aes(x = Week, y = Mean5.Weight, 
                                                         color = as.factor(Trt))) 
fig2 <- fig2 +
  geom_point(aes(x = Week, y = Mean5.Weight, group = as.factor(Tank), 
                 color = as.factor(Trt)), size = 2, shape = 19, alpha = 0.25, position = position_dodge(width = dg)) +
  geom_line(aes(x = Week, y = Mean5.Weight, group = as.factor(Tank), 
                color = as.factor(Trt)), alpha = 0.25, position = position_dodge(width = dg)) +
  scale_color_discrete(type = colors, name = "Care Treatment") +
  stat_summary(fun.y = mean, fun.ymin = meanminus, fun.ymax = meanplus, 
               geom = "pointrange", size = 1.5, lwd = 1.5, position = position_dodge(width = dg)) +
  stat_summary(fun = mean, geom = "line", size = 1.25, 
               position = position_dodge(width = dg)) +
  ylab("Mean Weight of 5 Fish (g)") + xlab("Week Since Fertilization") + theme_bw()

fig2  

