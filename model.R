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

intake_df <- read.csv("intake.csv", header=TRUE)
git_ethanol <- intake(dataframe = intake_df, lambda = 0.0015)
blood_ethanol <- blood(git_ethanol, activity = 0.0015)

png('plot_git_blood.png')
par(mfrow=c(2,1))
time <- 1:length(git_ethanol)/3600
plot(time, git_ethanol, type='l')
time <- 1:length(blood_ethanol)/3600
plot(time, blood_ethanol, type='l')
dev.off()

tester <- read.csv("tester.csv", header=TRUE)
start <- intake_df$time[1]
tester$seconds <- difftime( tester$time, start, units="secs" )
tester$hours <- difftime( tester$time, start, units="hours" )

vaha <- 77
prom <- blood_ethanol / vaha

#predicted data by the model:
model_data <- c( 0, prom[tester$seconds])
real_data  <- tester$breath

png('plot.png')
plot(tester$hours, tester$breath)
lines(time, prom)
dev.off()

plot(model_data, col=2)
points(real_data, col=3)
data_matrix <- matrix( c( model_data, real_data ), ncol=2 )

print("Distance between real and model data:")
print(dist( t( data_matrix) ))
