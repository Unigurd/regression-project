library(riskCommunicator)
data(framingham, package="riskCommunicator")

# This file is just a script to split the data into a training and a validation set.
set.seed(8088)
n <- nrow(framingham)
validation_size <- floor(n/10)
i <- sample(n, validation_size)
validation_data <- framingham[i,]
training_data <- framingham[-i,]
write.csv(validation_data, "validation_data.csv")
write.csv(training_data, "training_data.csv")
