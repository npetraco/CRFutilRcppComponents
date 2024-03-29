library(CRFutil)
library(CRFutilRcppComponents)
library(rbenchmark)
library(microbenchmark)

# Make up a random graph
num.nodes <- 12
g <- erdos.renyi.game(num.nodes, 1, typ="gnp")
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

plot(1:nrow(configs.and.counts),configs.and.counts[,ncol(configs.and.counts)],typ="h", ylab="config freqs",xlab="config #")

# Here we remove the extra index put in by CRF using the C function
theta.pars   <- fix_node_and_edge_par(node_par = rmod$node.par, edge_par = rmod$edge.par)

#delta_alpha(as.matrix(samps), rmod$node.par, rmod$edge.par, num_params_in=0)
Da.matC <- delta_alpha(
  samples       = as.matrix(samps), 
  node_par      = theta.pars$node_par, 
  edge_par      = theta.pars$edge_par,
  edge_mat      = rmod$edges,
  adj_nodes     = rmod$adj.nodes, 
  num_params_in = 0)

Delta.alpha.info <- delta.alpha(crf = rmod, samples = samps, printQ = F)
Delta.alpha <- Delta.alpha.info$Delta.alpha


dim(Da.matC)
dim(Delta.alpha)

# Check C vs R output
sum(Da.matC[,1] != Delta.alpha[,1])
sapply(1:ncol(Da.matC), function(xx){sum(Da.matC[,xx] != Delta.alpha[,xx])})

rmod$node.par
class(samps)

theta.pars$node_par
testslicedetect(theta.pars$node_par)
testslicedetect(rmod$node.par)
