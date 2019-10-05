library(CRFutil)
library(rbenchmark)
library(microbenchmark)

# Make up a random graph
g <- erdos.renyi.game(16, 1, typ="gnp")
dev.off()
plot(g)

# Get its adjacency matrix and genrate an MRF sample
adj <- as.matrix(as_adj(g))

f0       <- function(y){ as.numeric(c((y==1),(y==2)))}
rmod     <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

# "true" theta
rmod$par <- runif(rmod$n.par,-1.5,1.5)
rmod$par

# Make true pots from true theta
out.pot <- make.pots(parms = rmod$par,  crf = rmod,  rescaleQ = T, replaceQ = T)
rmod$edges
rmod$node.pot
rmod$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 500
samps <- sample.exact(rmod, num.samps)
colnames(samps) <- 1:ncol(samps)
mrf.sample.plot(samps)

configs.and.counts <- as.data.frame(ftable(data.frame(samps)))
head(configs.and.counts)
configs <- configs.and.counts[,1:ncol(samps)]
configs <- as.matrix(sapply(configs, as.numeric))
head(configs)


theta.pars   <- fix_node_and_edge_par(node_par = rmod$node.par, edge_par = rmod$edge.par)
jridx <- 1
jridx <- sample(1:nrow(configs), size = 1)
configs[jridx,]

test.phi.vec <- phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par)
test.phi.vec
#test.phi.vec <- NULL



test.phi.vec2 <- phi.features(config = configs[jridx,], 
             edges.mat = rmod$edges, 
             node.par = rmod$node.par, 
             edge.par = rmod$edge.par, 
             ff=f0)
length(test.phi.vec2) == sum(test.phi.vec == test.phi.vec2)

# Check to make sure we get the same phi vectors between the R and C code: 
chk.phis <- rep(FALSE, nrow(configs))
for(i in 1:nrow(configs)) {
  jridx <- i

  test.phi.vec  <- phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par)
  test.phi.vec2 <- phi.features(config = configs[jridx,], 
                                edges.mat = rmod$edges, 
                                node.par = rmod$node.par, 
                                edge.par = rmod$edge.par, 
                                ff=f0)
  
  chk.phis[i] <- length(test.phi.vec2) == sum(test.phi.vec == test.phi.vec2)
  
}
# Number of mistakes:
sum(!chk.phis)


# Speed check 1:
for(i in 1:nrow(configs)) {
  jridx <- i
  
  test.phi.vec  <- phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par)

  print(i)
  
}

# Speed check 2:
benchmark(
  test.phi.vec  <- phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par),
  test.phi.vec2 <- phi.features(config = configs[jridx,], 
                                edges.mat = rmod$edges, 
                                node.par = rmod$node.par, 
                                edge.par = rmod$edge.par, 
                                ff=f0),replications = 1000
)


microbenchmark(
  phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par),
  phi.features(config = configs[jridx,], 
                                edges.mat = rmod$edges, 
                                node.par = rmod$node.par, 
                                edge.par = rmod$edge.par, 
                                ff=f0)
)

526.45265/30.64928

#fix_node_and_edge_par2(node_par = rmod$node.par, edge_par = rmod$edge.par)

# C loop: (much faster, but num_par pass in causing a bug....)
junk <- compute_model_matrix(configs, rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par)
junk <- NULL
dim(junk)
head(junk)


# R loop:
rmod$n.par
junk <- t(sapply(1:nrow(configs), function(xx){
  #phi.features(config = configs[xx,], edges.mat, node.par, edge.par, ff)
  phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par)  
}))


# Compare: 
benchmark(
  junk  <- compute_model_matrix(configs, rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par),
  junk2 <- t(sapply(1:nrow(configs), function(xx){
    phi_features_C(as.matrix(configs[jridx,]), rmod$edges, node_par = theta.pars$node_par, edge_par = theta.pars$edge_par)  
  }))
)

# RcppParallel ???
