filename = "Interactions Average.csv"
interactions = read.csv(filename, header = TRUE)

library(ggplot2)
ggplot(interactions, aes(x = Yeast)) + geom_bar()
       