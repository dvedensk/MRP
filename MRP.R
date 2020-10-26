library(rstan)
library(tigris)
library(readr)
library(dplyr)
library(stringr)
library(sampling)
library(ggthemes)
library(mapproj)

minnesota<-read.csv("ss14pmn.csv")

### Pre-process
###
minnesota <- minnesota %>% select("POVPIP","RAC1P", "SEX", "AGEP", "PUMA","PWGTP")
minnesota <- minnesota[complete.cases(minnesota),]     #drop NAs
minnesota$inPov <- as.numeric(minnesota$POVPIP < 100) #make categorical 0/1 where in Poverty if ratio < 100
minnesota$Sex <- factor(minnesota$SEX,levels=c(1,2),labels=c("Male","Female"))
minnesota$Race <- factor(minnesota$RAC1P,levels=1:9,labels=c("White","Black","AmerInd","Alaska","AmerInd2",
                                                              "Asian","PI","Other","TwoOrMore"))
minnesota$Age <- minnesota$AGEP %>% cut(breaks=c(0,17,40,64,95), include.lowest=T) #break into 4 age categories
minnesota$countyID <- as.numeric(as.factor(minnesota$PUMA))

minnesota <- minnesota %>% select(inPov, Sex, Race, Age, PUMA, countyID, PWGTP)

PUMA_ID_map <- minnesota %>% select(PUMA,countyID) %>% unique %>% arrange(countyID)
###

### PPS sample with weights w_i+beta+I(poverty)
###
N <- 10000
B <- 200 #degree of informative sampling
minnesota$PPS_weights <- minnesota$PWGTP + B*minnesota$inPov
inclusion_probs <- inclusionprobabilities(minnesota$PPS_weights, N)
indices <- UPmidzuno(inclusion_probs)
minn_sample <- minnesota[indices==1,]
###

### Stan Fit
###
Xmat <- model.matrix(formula(~ (Sex+Race+Age) -1),
                      data=data.frame(minn_sample))

n_iter <- 4000
n_obs <- dim(Xmat)[1]
n_cov <- dim(Xmat)[2]
n_cty <- length(unique(minnesota$countyID))
cty <- minn_sample$countyID
y <- minn_sample$inPov
x <- Xmat
kappa <- 5
sigma_beta <- sqrt(10)

mrp_dat <- list(n_obs=n_obs, n_cov=n_cov, n_cty=n_cty, cty=cty,
                y=y, x=x, kappa=kappa, sigma_beta=sigma_beta)

mrp_fit0 <- stan("MRP.stan", data = mrp_dat, chains = 1, iter = 1)

mrp_fit <- stan(fit=mrp_fit0, data = mrp_dat, cores = 4, chains = 1,
                warmup = 2000, iter = n_iter, open_progress = FALSE)

save(mrp_fit, file="mrp_fit.Rdata")

betas_fit <- extract(mrp_fit)$beta
intercepts_fit <- extract(mrp_fit)$u
colnames(betas_fit) <- colnames(Xmat)
colnames(intercepts_fit) <- PUMA_ID_map[,1]
###

### Plotting
###
pumasMN <- pumas("MN", cb=T, class="sp")
pumasMN <- fortify(pumasMN, region='PUMACE10')
pumasMN$id<-gsub("(?<![0-9])0+", "", pumasMN$id, perl = TRUE) #remove leading 0's so we can join

#plot truth
#get total number of people in poverty per PUMA
PUMApov <- minnesota %>% group_by(PUMA,inPov) %>% filter(inPov==1) %>% 
                         add_tally(inPov, name="numInPov")  %>% 
                         select(PUMA,numInPov) %>% unique
#get total number of people per PUMA
PUMApop <- minnesota %>% group_by(PUMA) %>%  add_tally(name="PUMApop")  %>% 
                         select(PUMA,PUMApop) %>% unique

pov_vals <- PUMApov %>% left_join(PUMApop, by="PUMA") %>% arrange(PUMA)
pov_vals$PUMA <- as.character(pov_vals$PUMA)
pov_vals <- pov_vals %>% mutate(propPov=numInPov/PUMApop)

plotDF.True <- pumasMN %>% left_join(pov_vals, by=c("id"="PUMA"))
png("plots/truth.png")
ggplot(plotDF.True, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=0.3, aes(fill=propPov))+
  coord_map()+
  theme_map()+
  scale_fill_gradientn(colors=topo.colors(10))+
  ggtitle("True Poverty by County")
dev.off()

#plot Horvitz Thompson
minn_sample$PUMA <- as.character(minn_sample$PUMA)
PUMApop <- minn_sample %>% group_by(PUMA) %>% tally(PPS_weights,name="PUMApop")
PUMApov <- minn_sample %>% group_by(PUMA) %>% mutate(weightedPov=inPov*PPS_weights) %>%
                           tally(weightedPov, name="totWgtPov")

HT_vals <- PUMApov %>% left_join(PUMApop, by="PUMA") %>% arrange(PUMA)
HT_vals <- HT_vals %>% mutate(propPov=totWgtPov/PUMApop)

plotDF.HT <- pumasMN %>% left_join(HT_vals, by=c("id"="PUMA"))
png("plots/directHT.png")
ggplot(plotDF.HT, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=0.3, aes(fill=propPov))+
  coord_map()+
  theme_map()+
  scale_fill_gradientn(colors=topo.colors(10))+
  ggtitle("H-T estimate of Poverty")
dev.off()

### Post-stratify and plot
###
attach(minnesota)
cells <- expand.grid(levels(Sex), levels(Race),levels(Age),levels(as.factor(countyID)))
postStratCells <- minnesota %>% group_by(Race,Sex,Age,countyID,PUMA) %>% tally(name="cellCount")

Xpred <- model.matrix(formula(~ (Sex+Race+Age) -1),
                             data=data.frame(minnesota))

n_iter <- n_iter - 2000 #correct for warm up iters that weren't kept
logit.inv <- plogis 

#generate bernoulli predictions for every cell and every iteration
predictions <- matrix(NA,n_iter,nrow(minnesota))
for (iter in 1:n_iter) {
  predictions[iter,] <- logit.inv(Xpred %*% betas_fit[iter,] +  intercepts_fit[iter,minnesota$countyID])
  #Very slow!
  predictions[iter,] <- sapply(predictions[iter,], rbinom, n=1, size=1)
  #equivalent to sapply(rownames(alt_predictions), rbinom, n=1, size=1)? Or need an anonymous function?
}

#calculate proportions from the above bernoulli results
#should be vectorized to be quicker
predProp <- matrix(NA,n_iter, n_cty)
for (iter in 1:n_iter) {
  PUMA_pred <- data.frame(ID=minnesota$countyID, pred=predictions[iter,])
  total <- PUMA_pred %>% group_by(ID) %>% tally(name="PUMAtot")
  inPov <- PUMA_pred %>% group_by(ID) %>% filter(pred==1) %>% tally(name="PUMApov")
  tmp <- total %>% left_join(inPov, by="ID") %>% mutate(predProp=PUMApov/PUMAtot) %>% select(predProp) %>% as.list()
  predProp[iter,] <- tmp$predProp
}

#take the mean for every PUMA so we can plot
posteriorMeans <- data.frame(ID=1:n_cty, propPov=colMeans(predProp)) %>% 
                  left_join(PUMA_ID_map, by=c("ID"="countyID")) %>% select(-ID)
posteriorMeans$PUMA <- as.character(posteriorMeans$PUMA)

png("plots/poststrat.png")
plotDF.PS <- pumasMN %>% left_join(posteriorMeans, by=c("id"="PUMA"))
ggplot(plotDF.PS, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=0.3, aes(fill=propPov))+
  coord_map()+
  theme_map()+
  scale_fill_gradientn(colors=topo.colors(10))+
  ggtitle("Poststratification estimate of Poverty")
dev.off()

##MSE calculations
Y <- pov_vals$propPov

##H-T MSE
Y.hat <- HT_vals$propPov
MSE.HT <- sum((Y-Y.hat)^2)/n_cty #0.1168035

##Post-strat MSE
Y.hat <- posteriorMeans$propPov
MSE.PS <- sum((Y-Y.hat)^2)/n_cty #0.02240517
###

### Plot everything on a common gradient--truth ends up looking terrible
###
png("plots/common_gradient_truth.png")
ggplot(plotDF.True, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=0.3, aes(fill=propPov))+
  coord_map()+
  theme_map()+
  scale_fill_gradientn(colors=topo.colors(10), limits=c(0,.8))+
  ggtitle("True Poverty by County")
dev.off()

png("plots/common_gradient_HT.png")
ggplot(plotDF.HT, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=0.3, aes(fill=propPov))+
  coord_map()+
  theme_map()+
  scale_fill_gradientn(colors=topo.colors(10), limits=c(0,.8))+
  ggtitle("H-T estimate of Poverty")
dev.off()

png("plots/common_gradient_poststrat.png")
ggplot(plotDF.PS, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=0.3, aes(fill=propPov), limits=c(0,.8))+
  coord_map()+
  theme_map()+`
  scale_fill_gradientn(colors=topo.colors(10))+
  ggtitle("Poststratification estimate of Poverty")
dev.off()
