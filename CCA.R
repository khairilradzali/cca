#Load the package
library(vegan)

#Read OTU table and run CCA
data <- read.csv("varespec.csv", row.names=1, header = TRUE)
my.ca <- cca(data)

#Visualize the proportion of variance explained by each CCA axis
barplot(my.ca$CA$eig/my.ca$tot.chi, names.arg = 1:my.ca$CA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="CA axis")
head(my.ca$CA$eig/my.ca$CA$tot.chi)

#Visualize the scatterplot of OTUs data
plot(my.ca)

#Take a look at different scaling for CCA
layout(matrix(1:4,2,2,byrow=TRUE))
plot(my.ca, scaling="none", main="scaling = 0, raw")
plot(my.ca, scaling="sites", main="scaling = 1, sites")
plot(my.ca, scaling="species", main="scaling = 2, species (default)")
plot(my.ca, scaling="symmetric", main="scaling = 3, both")

#2 panel to display the sitees and species on separate plots
layout(matrix(1,2))
# plot sites
site.cols = "blue"
plot(my.ca, display="sites", type="n", main="Sites", scaling="sites")
text(my.ca, display="sites", col=site.cols, scaling="sites")
# plot species
plot(my.ca, display="species", type="n", ylab="", main="Species", scaling="species")
text(my.ca, display="species", col="black", scaling="species")

#Read physicochemical input table
fixed <- read.csv("varechem.csv", row.names=1, header = TRUE)

#Combine abundance table and physicochemical and run CCA
my.cca <- cca(data~., data=fixed)

#Check for variance inflation factors (VIF) for multicollinearity, VIF above than 10 considered as redundancy predictor variables.
vif.cca(my.cca)

#We use drop1 function as indicator to remove redundant/problematic predictor variables
drop1(my.cca, test="perm")


fixed <- fixed[,-8]

# rerun the cca 
my.cca <- cca(data ~ ., data=fixed)

# VIFs
vif.cca(my.cca)

#Compare the change in AIC and drop, if necessary
drop1(my.cca, test="perm")

#Remove any problematic variables
fixed <- fixed[,-8]

# rerun the cca, check VIF, and inspect AIC if necessary
my.cca <- cca(data ~ ., data=fixed)
vif.cca(my.cca) #if VIF>10, inspect the AIC again
drop1(my.cca, test="perm")

#Remove any problematic variable
fixed <- fixed[,-7]

# rerun the cca, check VIF, and inspect AIC if necessary 
my.cca <- cca(data ~ ., data=fixed)
vif.cca(my.cca) #if VIF>10, inspect the AIC again
drop1(my.cca, test="perm")

#Remove any problematic variable
fixed <- fixed[,-1]

# rerun the cca, check VIF, and inspect AIC if necessary
my.cca <- cca(data ~ ., data=fixed)
vif.cca(my.cca)
drop1(my.cca, test="perm")

#Remove any problematic variable
fixed <- fixed[,-11]

# rerun the cca, check VIF, and inspect AIC if necessary
my.cca <- cca(data ~ ., data=fixed)
vif.cca(my.cca) #Now all VIF<10

#Visualize the proportion of variance explained by each CCA axis
barplot(my.cca$CA$eig/my.cca$tot.chi, names.arg = 1:my.cca$CA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="CCA axis")

#CCA plotting based on different scaling
# plot none
plot(my.cca, display=c("sites", "bp"), type="n", main="none", scaling="none")
text(my.cca, display="sites", col=site.cols, scaling="none", cex = 0.8)
text(my.cca, display="species", col="black", scaling="none", cex = 0.8)
text(my.cca, display="bp", col="red", cex = 0.8)

# plot sites
plot(my.cca, display=c("sites", "bp"), type="n", main="Sites", scaling="sites")
text(my.cca, display="sites", col=site.cols, scaling="sites")
text(my.cca, display="bp", col="red")

# plot species
plot(my.cca, display=c("species", "bp"), type="n", ylab="", main="Species", scaling="species")
text(my.cca, display="species", col="black", scaling="species")
text(my.cca, display="bp", col="red")

# plot symmetric
plot(my.cca, display=c("sites", "bp"), type="n", main="symmetric", scaling="symmetric")
text(my.cca, display="sites", col=site.cols, scaling="symmetric", cex = 0.8)
text(my.cca, display="species", col="black", scaling="symmetric", cex = 0.8)
text(my.cca, display="bp", col="red", cex = 0.8)

# Extract the percent of proportion of variance explained by predictor variables
my.cca$CCA$tot.chi/my.cca$tot.chi
#59.68% of variance in community composition

# Significance test for individual predictors (type 3 test)
anova(my.cca, by="margin") 
# Strongest signifcant predictor variable is "S"

# significance test for entire model
anova(my.cca)
#The full mode is statistically significant (p<0.01)