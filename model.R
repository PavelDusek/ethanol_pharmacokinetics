#####################################################
# Two-compartment model of ethanol pharamcokinetics #
#####################################################

# two compartments: gastrointestinal tract (GIT) and blood
# GIT: 1st order kinetics with activity constant lambda
# blood: 0 order kinetics with activity constant activity

# for upgrades see https://pubmed.ncbi.nlm.nih.gov/3319346/

git <- function(N0 = 25, lambda = 0.0015, length = 10000) {
  # model 1st order kinetics of adsorption from GIT to blood with constant lambda
  # numerical model, time in seconds
  git_ethanol <- vector()
  git_ethanol[1] <- N0
  for (i in 2:length) {
    git_ethanol[i] <- git_ethanol[i-1] - git_ethanol[i-1] * lambda
    if (git_ethanol[i] < 0) {
      git_ethanol[i] <- 0
    }
  }
  git_ethanol
}

blood <- function(git, activity = 0.0015, N0 = 0 ) {
  # model 0-order kinetics of methabolism in blood constant lambda
  # numerical model, time in seconds
  blood_ethanol <- vector()
  blood_ethanol[1] <- N0 + git[2] - git[1]
  for (i in 2:(length(git)-1)) {
    blood_ethanol[i] <- blood_ethanol[i-1] + max(git[i] - git[i+1], 0) - activity
    if (blood_ethanol[i] < 0) {
      blood_ethanol[i] <- 0
    }
  }
  blood_ethanol
}

intake <- function( dataframe, lambda ) {
  #parse dataframe with data of alcohol intake, count time differences between data in seconds
  #and model alcohol content in GIT with constant lambda
  #dataframe is a dataframe with columns "time" and "intake"

  no_of_rows <- dim(dataframe)[1]
  dataframe$time <- as.POSIXct(dataframe$time)

  #first alcohol intake
  timediff <- difftime( dataframe$time[2], dataframe$time[1], units="secs" )
  git_ethanol <- git( N0 = dataframe$intake[1], lambda = lambda, length = timediff )

  #subsequent intakes
  for (i in 2:(no_of_rows - 1) ) {
    time_btwn_intakes <- difftime( dataframe$time[i+1], dataframe$time[i], units="secs" )
    last_value <- git_ethanol[ length(git_ethanol) ]
    actual <- git( N0 = ( dataframe$intake[i] + last_value ), lambda = lambda, length = time_btwn_intakes )
    git_ethanol <- c(git_ethanol, actual)
  }
  #last intake
  last_value <- git_ethanol[ length(git_ethanol) ]
  actual <- git( N0 = ( dataframe$intake[no_of_rows] + last_value ), lambda = lambda )
  git_ethanol <- c(git_ethanol, actual)
  git_ethanol
}

ethanol_model <- function( intake_df, git_constant, blood_constant, body_weight  ) {
  #body weight in kilograms
  git_ethanol <- intake(dataframe = intake_df, lambda = git_constant)
  blood_ethanol <- blood(git_ethanol, activity = blood_constant)
  git_ethanol <- git_ethanol[1:(length(git_ethanol)-1)] #remove last element so that the vector sizes match
  model_df <- data.frame(
    time     = 1:length(git_ethanol)/3600, #time in hours
    git      = git_ethanol,
    blood    = blood_ethanol,
    permille = blood_ethanol / body_weight #per mille of ethanol weight concentration in blood
  )
  model_df
}

#Ethanol intake data recorded
intake_df <- read.csv("intake.csv", header=TRUE)

#Ethanol concentration data measured:
tester <- read.csv("tester.csv", header=TRUE)
tester$time <- as.POSIXct(tester$time)
start <- tester$time[1]
tester$seconds <- difftime( tester$time, start, units="secs" )
tester$hours <- difftime( tester$time, start, units="hours" )
real_data  <- tester$breath

#Some good constants estimates:
#model <- ethanol_model( intake_df = intake_df, git_constant = 0.0015, blood_constant = 0.0015, body_weight = 77)

#Find model constants that best fits the real data (minimum eucledian distance):
git_constant_tested   <- vector()
blood_constant_tested <- vector()
distances             <- vector()
for (git_constant in seq(from = 0.00300, to = 0.00400, by = 0.00001 ) ) {
  for (blood_constant in seq(from = 0.00100, to = 0.00200, by = 0.00001 ) ) {
    git_constant_tested   <- c(  git_constant_tested,   git_constant)
    blood_constant_tested <- c(blood_constant_tested, blood_constant)

    model <- ethanol_model( intake_df = intake_df, git_constant = git_constant, blood_constant = blood_constant, body_weight = 77)
    model_data <- c( 0, model$permille[tester$seconds])
    data_matrix <- matrix( c( model_data, real_data ), ncol=2 )
    distance <- (dist( t( data_matrix) ))
    distances <- c(distances, distance)
    #print( c( git_constant, blood_constant, distance ) )
  }
}
distance_df <- data.frame( git = git_constant_tested, blood = blood_constant_tested, distance = distances)
write.csv(distance_df, "distance.csv", row.names=FALSE)

min_dist <- which.min(distance_df$distance)
print("Best model fit:")
print(distance_df[min_dist,])

#Show distance of real data and model data based on different constants:
library(plot3D)
png('distance.png')
scatter3D( distance_df$git, distance_df$blood, distance_df$distance, main="Model error", xlab="absorption rate [1/s]", ylab="metabolism rate [g/s]", zlab="error [1]" )
dev.off()
git_const <- distance_df[min_dist,1]
blood_const <- distance_df[min_dist,2]
#For best blood constant
distance_best_blood <- distance_df[distance_df$blood == blood_const, ]
png('distance_best_blood.png')
plot( distance_best_blood$git, distance_best_blood$distance, main="Model error", xlab="absorption rate [1/s]", ylab="error [1]" )
dev.off()
distance_best_git <- distance_df[distance_df$git == git_const,]
png('distance_best_git.png')
plot( distance_best_git$blood, distance_best_git$distance, main="Model error", xlab="metabolism rate [g/s]", ylab="error [1]" )
dev.off()

model <- ethanol_model( intake_df = intake_df, git_constant = git_const, blood_constant = blood_const, body_weight = 77)
model_data <- c( 0, model$permille[tester$seconds])

png('plot_git_blood.png')
par(mfrow=c(2,1))
plot(model$time, model$git, type='l', main="Ethanol content in gastrointestinal tract (model)", xlab="time [s]", ylab="ethanol content [g]" )
plot(model$time, model$blood, type='l', main="Ethanol content in blood (model)", xlab="time [s]", ylab="ethanol content [g]" )
dev.off()

png('plot.png')
plot(tester$hours, tester$breath, main="Ethanol content in blood (model) and breath (measured)", xlab="time [h]", ylab="ethanol content [permille]")
lines(model$time, model$permille)
dev.off()

#plot(model_data, col=2)
#points(real_data, col=3)
