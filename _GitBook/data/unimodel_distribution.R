library(tidyverse)

distr <- data.frame(temp = rnorm(1000, mean = 12, sd = 3))
ggplot(data = distr, aes(x = temp)) +
  geom_density(size = 0.8, colour = "red3") +
  xlab("Temperature (°C)") + ylab("Density")
ggsave(filename = "data/unimodel_distribution.jpg", width = 6, height = 3)
