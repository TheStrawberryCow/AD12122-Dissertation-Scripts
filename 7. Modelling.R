#-----------------------------------------------------
# Title: Modelling
# Date: 11/08/2023
# Author: Arif Ayman
# Description: Modelling and Prediction with PCs
#-----------------------------------------------------

#Generate data.frame
clr_pcs <- data.frame(
  "pc1" = ord_clr$CA$u[,1],
  "pc2" = ord_clr$CA$u[,2],
  "pc3" = ord_clr$CA$u[,3],
  "Group" = phyloseq::sample_data(ps_clr)$Group
)
clr_pcs$Group_num <- ifelse(clr_pcs$Group == "IE_High_Risk", 0, 1)
head(clr_pcs)


#Distributional characteristics 
dd <- datadist(clr_pcs)
options(datadist = "dd")

#Plot unconditional associations
a <- ggplot(clr_pcs, aes(x = pc1, y = Group_num)) +
  Hmisc::histSpikeg(Group_num ~ pc1, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC1", y = "Pr(Definite_IE)\n")
b <- ggplot(clr_pcs, aes(x = pc2, y = Group_num)) +
  Hmisc::histSpikeg(Group_num ~ pc2, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC2", y = "Pr(Definite_IE)\n")
c <- ggplot(clr_pcs, aes(x = pc3, y = Group_num)) +
  Hmisc::histSpikeg(Group_num ~ pc3, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC3", y = "Pr(Definite_IE)\n")
cowplot::plot_grid(a, b, c, nrow = 2, ncol = 2, scale = .9, labels = "AUTO")



#Spline fitting
m1 <- rms::lrm(Group_num ~ rcs(pc1, 3) + rcs(pc2, 3) + rcs(pc3, 3), data = clr_pcs, x = TRUE, y = TRUE)

#Penalty search
pentrace(m1, list(simple = c(0, 1, 2), nonlinear = c(0, 100, 200)))

pen_m1 <- update(m1, penalty = list(simple = 0, nonlinear = 200))
pen_m1

ggplot(Predict(pen_m1))


#Bootstrapping#
#Obtain optimism corrected estimates
(val <- rms::validate(pen_m1))

#Corrected c-statistic
(c_opt_corr <- 0.5 * (val[1, 5] + 1))

#Plot calibration
cal <- rms::calibrate(pen_m1, B = 200)
plot(cal)
summary (cal)



#Output predicted probability
probabilities <- pen_m1 %>%
  predict(clr_pcs, type = "fitted")
print (probabilities)



#For new data
#probabilities <- pen_m1 %>% 
#predict(newdata, type = "fitted")
#print (probabilities)
