Nr_of_cores = 3 #total number of samples will be cores * n_samples
cl <- makeCluster(Nr_of_cores) 
registerDoParallel(cl) # register the cluster

popsize <- 100000	# size of population
n.max <- 30		# maximum number of participants per sample
n.min <- 20		# minimum sample size !!! HAS TO BE EQUAL TO VISIT LENGTH
items <- 5 #number of items per subject
n_samples <- 5  #Number of samples per item i.e number of bootstrapped trajectories in Sch?nbrodt terminology.

setwd("~/Corridor_of_stability")
covmat_data <- read.csv("autism_data.csv", sep = ";")
covmat_model <- lmer(CHI_MLU ~VISIT + Diagnosis + (1+VISIT|SUBJ),covmat_data)


#extracting covariance structure of the trained model
vc<-VarCorr(covmat_model) # variance component
sigmaSubject <- as.numeric(attr(vc[[1]],"stddev")[1]) # random intercept by subject
sigmaVisit <- as.numeric(attr(vc[[1]],"stddev")[2]) # random slope of visit over subject
sigmaResiduals <- as.numeric(attr(vc,"sc"))
sigmaCorrelation <- as.numeric(attr(vc[[1]],"correlation")[2])

#forimg covariance matrix
covmat<-matrix(c(sigmaSubject^2,
                 sigmaCorrelation*sigmaSubject*sigmaVisit,
                 sigmaCorrelation*sigmaSubject*sigmaVisit,
                 sigmaVisit^2),nrow=2)
rownames(covmat) <- colnames(covmat) <- c("intercept", "Visit")

#Sett
Intercept_val <- 2 #True mean of intercept
Visit_val <- 0.2       #True mean effect of Visit 
ASD_effect = 0     #True effect of diagnosis: ASD is 0, and TD is 0.5
TD_effect = 0.5
Distri = "gaussian"#Distribution of residuals
SD_res = 0.5    # SD of residuals for Gaussian distribution
item = 5


print(paste("\nComputing next population at ", Sys.time()))

### Generating population 

mlu_grid<- expand.grid(Child.ID = 1:popsize, Visit = 1:item)

mlu_grid$Child.ID <- factor(paste("Child.ID", mlu_grid$Child.ID, sep = ""))


mlu_grid <- arrange(mlu_grid, Child.ID)


mlu_grid$Diagnosis <- factor(rep(c(rep("ASD",item), rep("TD",item))),
                             levels= c("ASD", "TD"))


#Simulating response values
mlu_sim<- sim.glmm(design.data = mlu_grid,
                   fixed.eff = list(intercept = Intercept_val,
                                    Visit = Visit_val,
                                    Diagnosis = c(ASD = ASD_effect, TD = TD_effect)),
                   rand.V = list(Child.ID = covmat),
                   distribution = Distri,
                   SD = SD_res)


#DF to store coefficients of each sample.
sample_df <- c()



samp = 1

cl <- makeCluster(3)
registerDoParallel(cl)

res <- foreach(i = 1:3, #samples needs to be replaced here
              .combine = "cbind", 
              .packages = c("lmerTest","tidyverse","MuMIn")) %dopar% {
                
                max_items <- n.max * item		#the length of the sample data frame, subjects * items
                min_items <- n.min * item   #length of sample data frame at first calculation
                visit_nr = item    # number of items per subject
                
                # DF of subjects to sample from
                sample_people <- mlu_sim %>% 
                  filter(Visit == 1) 
                
                #Drawing random subjects from sample equal to max sample size n.max
                sam_people <- sample_people[sample(1:nrow(sample_people), size=n.max, replace=F), ]	
                
                #Adding all trials for each sampled subject
                sample_people <- mlu_sim %>% 
                  filter(mlu_sim$Child.ID %in% sam_people$Child.ID)
                
                
                
                #Estimates parameters for sample
                 result <- corEvol(sample_people, n.min = n.min, item_nr = item)
                
                
                #adding sample number to name of columns
                colnames(result) <- c(paste("item_nr",i, sep = ""), 
                                      paste("Intercept", i,sep = ""), 
                                      paste("Visit", i,sep = ""), 
                                      paste("DiagnosisTD", i, sep = ""),
                                      paste("R2m", i, sep = ""),
                                      paste("R2c", i, sep = ""))
                
                result
                
                
              }

stopCluster(cl)
              