####### testing the 'entropart' pacakage

# q = 0 corresponds to the weighted harmonic mean (reciprocal of the arithmetic mean)
# q = 1 to the weighted geometric mean
# q = 2 to the weighted arithmetic mean

library(entropart)
df <- data.frame(C1 = c(10, 10, 10, 10), C2 = c(0, 20, 35, 5), C3 = c(25, 15, 0, 2), row.names = c("sp1", "sp2", "sp3", "sp4"))
df


w <- c(1, 2, 1)
MC <- MetaCommunity(Abundances = df, Weights = w)
class(MC)

## The "P-sub-i" for Shannon entropy which quantifies the uncertainty in 
##predicting the species identity of an individual that is taken at random from the dataset
MC$Ps

###  =
data("Paracou618", package = "entropart")
summary(Paracou618.MC)

# It can also do phylogenetic and functional diversity

Tsallis(Ps = Paracou618.MC$Ps, q = 1)  # Order 1 is the weighted geometric mean (see above)
Shannon(Ps = Paracou618.MC$Ps)  # Notice that is the same!

# Sample coverage - probability for a species of the community to be observed in the actual sample.
Coverage(Ns = Paracou618.MC$Ns)

# Bias-corrected estimators - improve the estimation of entropy despite unobserved species
bcTsallis(Ns = Paracou618.MC$Ns, q = 1)

# Bias-corrected Shannon - This one is the one we should use!
bcShannon(Ns = Paracou618.MC$Ns)

# Effective numbers of species - Entropy converted into “true diversity” via Hill's numbers
expq(Simpson(Ps = Paracou618.MC$Ps), q = 2)
Diversity(Ps = Paracou618.MC$Ps, q = 2)

# Effective species number with bias correction
bcDiversity(Ns = Paracou618.MC$Ns, q = 2)

# Meta-community functions
# Alpha entropy
e <- AlphaEntropy(Paracou618.MC, q = 1)
summary(e)


# Diversity partitioning among alpha, beta and gamma diversity
p <- DivPart(q = 1, MC = Paracou618.MC, Biased = FALSE)
summary(p)

#  Diversity estimates for each component
de <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", Simulations = 1000)
summary(de)
plot(de)
#  Diversity profile
dp <- DivProfile(seq(0, 2, 0.2), Paracou618.MC, Biased = FALSE)
summary(dp)
plot(dp)
