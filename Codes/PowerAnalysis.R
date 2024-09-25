#install.packages('pwr')
library(pwr)
pwr.r.test(r = 0.519, sig.level = 0.05, power = 0.8)

# Calculate power for a correlation test with a given sample size
#result <- pwr.r.test(n = 18, r = 0.86, sig.level = 0.001)

# Print the result
#print(result)
pwr.r.test(n = 18, r = 0.86, sig.level = 0.001)
