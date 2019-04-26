###
### March 20, 2018
###
### John Heisey 400151293
###
### STATS780 Final Project:
### Analysis of breast cancer data
###
###
### Dataset Source: https://www.openml.org/d/13
### University Medical Centre,
### Institute of Oncology, Ljubljana, Yugoslavia.
###
### (If you are viewing in RStudio, all code is folded)
###

#####
### Data Prep/Importing Libraries
#####
library("ggthemes")
library("GGally")
library("extracat")
library("hdrcde")
library("KernSmooth")
library("ggplot2")
library("gridExtra")
library("vcd")
library(faraway)


# Importing dataframe from csv source file
# Note: There are a small percentage of rows with
# missing attributes (~8), to ensure a fair comparison
# between all methods these were discarded.
# Additionally the raw data was in all quote for
# with the quotes being removed in excel prior to
# importing
bc <- read.csv('breast-cancer-noquote.csv', header = TRUE, stringsAsFactors = FALSE)


# Initial excel cleaning process accidentally formatted
# 5-9 to 43229, switching it back
bc[bc$tumor.size == "43229",]$tumor.size <- "5-9"


# creating subset where there are only 
# complete observations (from 286 -> 278)
bc <- bc[complete.cases(bc), ]
#missed this data point, removing it
bc <- bc[-241,]



# Converting output recurrence variable to 0/1
# for no recurrance/recurrance, respectively
bc[bc$Class == 'no-recurrence-events',]$Class <- 0
bc[bc$Class == 'recurrence-events',]$Class <- 1
bc$Class <- as.integer(bc$Class)


#converting age to numeric values (20-20 -> 2, 30-39 -> 3, etc.)
for (ageGrp in c("20-29","30-39","40-49","50-59","60-69","70-79")){
  bc[bc$age == ageGrp,]$age <- substr(ageGrp,1,1)
}
bc$age <- as.integer(bc$age)

bc[bc$node.caps == "yes",]$node.caps <- 1
bc[bc$node.caps == "no",]$node.caps <- 0
bc$node.caps <- as.integer(bc$node.caps)

#Converting to easier to work with linear format
#(for loop would have been effective here but 
#cost/benefit for effort/time proved it not useful)
bc[bc$inv.nodes == "0-2",]$inv.nodes <- 1
bc[bc$inv.nodes == "3-5",]$inv.nodes <- 4
bc[bc$inv.nodes == "6-8",]$inv.nodes <- 7
bc[bc$inv.nodes == "9-11",]$inv.nodes <- 10
bc[bc$inv.nodes == "12-14",]$inv.nodes <- 13
bc[bc$inv.nodes == "15-17",]$inv.nodes <- 16
bc[bc$inv.nodes == "18-20",]$inv.nodes <- 19
bc[bc$inv.nodes == "21-23",]$inv.nodes <- 22
bc[bc$inv.nodes == "24-26",]$inv.nodes <- 25
bc$inv.nodes <- as.integer(bc$inv.nodes)

bc[bc$tumor.size == "0-4",]$tumor.size <- "2"
bc[bc$tumor.size == "5-9",]$tumor.size <- "7"
bc[bc$tumor.size == "10-14",]$tumor.size <- "12"
bc[bc$tumor.size == "15-19",]$tumor.size <- "17"
bc[bc$tumor.size == "20-24",]$tumor.size <- "22"
bc[bc$tumor.size == "25-29",]$tumor.size <- "27"
bc[bc$tumor.size == "30-34",]$tumor.size <- "32"
bc[bc$tumor.size == "35-39",]$tumor.size <- "37"
bc[bc$tumor.size == "40-44",]$tumor.size <- "42"
bc[bc$tumor.size == "45-49",]$tumor.size <- "47"
bc[bc$tumor.size == "50-54",]$tumor.size <- "52"
bc$tumor.size <- as.integer(bc$tumor.size)

#density plot of age variables with normal curve overlay
hist(bc$age, breaks=21, freq = FALSE, col = "lightblue", xaxt = "n", xlab = "Age Group", main = "")
axis(1, at=1:7, labels=c("","20-29","30-39","40-49","50-59","60-69","70-79"))
curve(dnorm(x, mean=mean(bc$age), sd=sd(bc$age)), col="darkblue", lwd=2, add=TRUE)

#converting radiation treatment variable to binary
bc[bc$irradiat == "yes",]$irradiat <- 1
bc[bc$irradiat == "no",]$irradiat <- 0
bc$irradiat <- as.integer(bc$irradiat)

# Doubledecker plot including node capsule penetration
# versus recurrence classification
doubledecker(node.caps ~ Class, data = bc, gp = gpar(fill = c("lightblue", "indianred1")), margins = c(2,5,2,1))


#####
### Logstic Regression Method
#####

# applying GLM to all linear variables
bc_glm_all = glm(Class ~ age + inv.nodes + tumor.size + deg.malig, data = bc, family = binomial("logit"))
summary(bc_glm)
bc_glm_no_age  = glm(Class ~ inv.nodes + tumor.size + deg.malig, data = bc, family = binomial("logit"))
summary(bc_glm_no_age)

glm_age = glm(Class ~ age, data = bc, family = binomial("logit"))
summary(glm_age)
glm_tumor = glm(Class ~ tumor.size, data = bc, family = binomial("logit"))
summary(glm_tumor)
glm_nodes = glm(Class ~ inv.nodes, data = bc, family = binomial("logit"))
summary(glm_nodes)
glm_malig = glm(Class ~ deg.malig, data = bc, family = binomial("logit"))
summary(glm_malig)

#Reducing model (no age) and testing null hypothesis
full_glm_res <- 345.90
red_glm_res <-308.98
G <- full_glm_res - red_glm_res
#chi squared result
1 - pchisq(G,2,lower.tail=FALSE)
#p value result
pchisq(G,2,lower.tail=FALSE)


#Reducing model further (no tumour/no node)
#and comparing
bc_glm_no_age  = glm(Class ~ inv.nodes + tumor.size + deg.malig, data = bc, family = binomial("logit"))
summary(bc_glm_no_age)
bc_glm_no_age_tum = glm(Class ~ inv.nodes + deg.malig, data = bc, family = binomial("logit"))
summary(bc_glm_no_age_tum)

bc_glm_no_age_node = glm(Class ~ tumor.size + deg.malig, data = bc, family = binomial("logit"))
summary(bc_glm_no_age_node)

bc_glm_no_age_node = glm(Class ~ tumor.size, data = bc, family = binomial("logit"))
summary(bc_glm_no_age_node)


full_glm_res <- 308.98
red_glm_res <- 311.97
G <- -full_glm_res + red_glm_res
#chi squared result
1 - pchisq(G,2,lower.tail=FALSE)
#p value result
pchisq(G,2,lower.tail=FALSE)




#Odds ratio calculations
d_size <- 5 #change in tumor size (mm)
odds_tumor_size <- exp(0.03741*d_size) #corresp. O.R.

d_num <- 4 #change in num. of nodes
odds_num_node <- exp(0.16898*d_num) #corresp. O.R.

d_deg <- 1 # change in degree
odds_deg <- exp(0.9622*d_num)


#logit plots produced using the reduced GLM
plot(bc$Class ~ bc$inv.nodes, xlab="Involved Nodes", ylab="Prob. of Recurrence")
lines(sort(bc$inv.nodes),ilogit(-1.32066+0.16898*sort(bc$inv.nodes)), lwd=4, col="lightblue")


#Generating train/test splits to use for all methods
#setting the seed to make partition reproductible
set.seed(780) # 780 seemed quite fitting
library(caTools)
train_rows = sample.split(bc$Class, SplitRatio=0.75)
train = bc[ train_rows,]
test  = bc[!train_rows,]

train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]

#Initiate classification test with stratified samples
bc_glm_no_age  = glm(Class ~ inv.nodes + tumor.size + deg.malig, data = train, family = binomial("logit"))
summary(bc_glm_no_age)

# beginning prediction
pred_classes<-round(predict(bc_glm_no_age,newdata=test,type="response"))
tab0 <- table(test$Class,pred_classes)
classAgreement(tab0)$crand

#recording all predictions for logistic regression
glm_pred_accs <- c(0.7323944, 0.7323944, 0.7323944, 0.6901408, 0.7183099, 0.7464789, 0.7323944, 0.7323944, 0.6901408, 0.7042254)
mean(glm_pred_accs)
sd(glm_pred_accs)



######
### Classification Tree Analysis Data Prep.
######

#importing analysis specific libraries
library(rpart) # for the tree
library(RColorBrewer)
library(ggplot2)
library(randomForest)
library(MASS)
library(tree)
library(gbm)
library(e1071)
library("ggthemes")
library("GGally")
library("extracat")
library(hdrcde)
library(KernSmooth)
library("gridExtra")
library("vcd")
library(tree)

#factorizing variables for classification tree
bc$node.caps <- factor(bc$node.caps)
bc$breast <- factor(bc$breast)
bc$breast.quad <- factor(bc$breast.quad)
bc$menopause <- factor(bc$menopause)
bc$Class <- factor(bc$Class)
bc$irradiat <- factor(bc$irradiat)

#running gen. classification analysis
bc_tree <- rpart(Class ~ ., data=bc, method="class")
print(bc_tree)
plot(bc_tree,margin = 0.03)
text(bc_tree, cex=1, pretty = 1)

#variable for class plot colouring
recur <- bc$Class

#plotting two initial splitter variables
#coloured by class
ggplot(bc, aes(x=deg.malig, y=inv.nodes)) +
  geom_point(size=3, shape=23, aes(color = recur)) +
  geom_smooth(method = lm) +
  theme(panel.background = element_blank())

#partition plot (not used in report, not very informative)
partition_plot <- tree(Class ~ inv.nodes + deg.malig, data = bc)
partition.tree(seed_plot1)

#####
### General classification tree analysis
#####
bc_tree <- tree(Class ~ ., data=train, method="class")
pred_class <- predict(bc_tree, test, type = "class")

#Generate classification table and analyze
tab1 <- table(test[,10],pred_class)
1 - classAgreement(tab1)$diag
classAgreement(tab1)$diag
tab1
classAgreement(tab1)$crand

#recorded prediction accuracies for the 10 splits
gen_class_pred_acc <- c(0.7391304, 0.6714286, 0.6760563, 0.7183099, 0.7042254, 0.7042254, 0.7391304, 0.7183099, 0.6760563, 0.6901408)
mean(gen_class_pred_acc)
sd(gen_class_pred_acc)


#####
###Bagging (identical analysis process to general class.)
#####
train <- train[complete.cases(train), ]
bag.bc = randomForest(Class~.,data=train,mtry=9,importance=TRUE, type="class")
bag.bc
bc.pred.bag = predict(bag.bc,test,type ="class")
tab2<-table(test[,10],bc.pred.bag)
tab2
classAgreement(tab2)$diag
1-classAgreement(tab2)$diag
classAgreement(tab2)$crand

# Variable importance plots produced
importance(bag.bc)
Bag.Var.Importance <- bag.bc
varImpPlot(Bag.Var.Importance)

# Recorded prediciton accuracies for bagging method
bag_pred_acc = c(0.7323944, 0.7352941, 0.6764706, 0.6521739, 0.7183099, 0.7352941, 0.7323944, 0.6521739, 0.7323944, 0.7142857)
mean(bag_pred_acc)
sd(bag_pred_acc)


#####
### Random Forest Method (identical analysis process to bagging)
#####
train <- train[complete.cases(train), ]
rf.bc = randomForest(Class~.,data=train,mtry=3,importance=TRUE, type="class")
rf.bc

bc.pred.rf = predict(rf.bc,test,type ="class")
tab3<-table(test[,10],bc.pred.rf)
tab3
classAgreement(tab3)$diag
1-classAgreement(tab3)$diag
classAgreement(tab3)$crand

importance(rf.bc)
Bag.Var.Importance <- rf.bc
varImpPlot(Bag.Var.Importance)

#recorded RF accuracies
rf_pred <- c(0.7571429, 0.7571429, 0.7058824, 0.7571429, 0.7681159, 0.7794118, 0.7681159, 0.7536232, 0.7794118, 0.6857143)
mean(rf_pred)
sd(rf_pred)










#####
### Boosting method (identical analysis process to RF)
#####
train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]
boost.bc = gbm(Class~., data = train, distribution = "multinomial", n.trees = 500, interaction.depth = 3, shrinkage = 0.003)
par(las=1, cex=1, mar = c(5,6,1,2))
summary(boost.bc)
yhat.boost=predict(boost.bc,newdata=test,n.trees=500,distribution="multinomial",type="response")
class.pred<-rep(0,nrow(test))
for(i in 1:nrow(test)){
  which(yhat.boost[i,,1]==max(yhat.boost[i,,1]))->class.pred[i]
}
tab4<-table(test[,10],class.pred)
tab4
1-classAgreement(tab4)$diag
classAgreement(tab4)$diag
gbm.perf(boost.bc)

boost.pred1 = c(0.7571429, 0.826087, 0.7826087, 0.7826087, 0.7826087, 0.7352941, 0.7571429, 0.7352941,0.7391304, 0.7536232)
mean(boost.pred1)
sd(boost.pred1)


### Boosting with revised parameters
train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]
boost2.bc = gbm(Class~., data = train, distribution = "multinomial", n.trees = 237, interaction.depth = 3, shrinkage = 0.003)

summary(boost2.bc)
yhat.boost2=predict(boost2.bc,newdata=test,n.trees=273,distribution="multinomial",type="response")
class.pred2<-rep(0,nrow(test))
for(i in 1:nrow(test)){
  which(yhat.boost2[i,,1]==max(yhat.boost2[i,,1]))->class.pred2[i]
}
tab5<-table(test[,10],class.pred2)
tab5
1-classAgreement(tab5)$diag
classAgreement(tab5)$diag
classAgreement(tab5)$crand
gbm.perf(boost2.bc)

boost.pred2 = c(0.7391304, 0.7246377, 0.8428571, 0.8142857, 0.7352941, 0.7391304, 0.761194, 0.7285714, 0.7142857, 0.7714286)
mean(boost.pred2)
sd(boost.pred2)


#####
### Neural Network Classification method
#####

#sigmoid function used in NN plotted for explanatory purposes 
sigmoid = function(v) {
  1 / (1 + exp(-v))
}
v <- seq(-7,7,0.01)

plot(v, sigmoid(v), col='indianred1', xlab = "v", ylab = "sigma(v)", yaxt = "n")
axis(2, at=c(0,0.5,1))


library(nnet)
library(neuralnet)
train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]

train$age <- scale(train$age)
train$tumor.size <- scale(train$tumor.size)
train$inv.nodes <- scale(train$inv.nodes)
train$deg.malig <- scale(train$deg.malig)
test$age <- scale(test$age)
test$tumor.size <- scale(test$tumor.size)
test$inv.nodes <- scale(test$inv.nodes)
test$deg.malig <- scale(test$deg.malig)

nn_bc<-nnet(Class ~ deg.malig + tumor.size + inv.nodes + age, data = train, size=6, decay=0.25, maxit=200)

NN = neuralnet(Class ~ deg.malig + tumor.size + inv.nodes + age, data=train, hidden = 6)
plot(NN)

nn_pred<-predict(nn_bc, test, type="class")
tab6 <- table(test[,10],nn_pred)
tab6
1-classAgreement(tab6)$diag
classAgreement(tab6)$diag
classAgreement(tab6)$crand


nn_acc1 <- c(0.7464788732, 0.7826086957, 0.7714285714, 0.7246376812, 0.768115942, 0.7285714286,  0.8285714286, 0.7285714286, 0.776119403)
mean(nn_acc1)
sd(nn_acc1)


#####
### VSCC method
#####
library(vscc)
bc_vscc <- bc[complete.cases(bc),]
bc_vscc <- bc_vscc[,-c(2,7,8)]
bc_vscc$node.caps <- as.integer(bc_vscc$node.caps)
bc_vscc$irradiat <- as.integer(bc_vscc$irradiat)
bc_vscc$Class <- as.integer(bc_vscc$Class)

x <- scale(bc_vscc[,-7])

run<-vscc(x, G = 1:3)
run

head(run$topselected)

tab7<-table(bc_vscc[,7],run$bestmodel$classification)
1-classAgreement(tab7)$diag
classAgreement(tab7)$crand

plot(run)


