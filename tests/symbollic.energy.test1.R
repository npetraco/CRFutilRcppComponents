library(CRFutil)
library(CRFutilRcppComponents)
library(rbenchmark)
library(microbenchmark)

# Make up a random graph
g <- erdos.renyi.game(5, 1, typ="gnp")
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

# Tests:
tesf <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

#jridx <- 1
jridx <- sample(1:nrow(configs), size = 1)

# R version:
symbolic.conditional.energy(config                   = as.matrix(configs[jridx,]), 
                            condition.element.number = 1, 
                            crf                      = tesf,
                            ff                       = f0,
                            format                   = "conditional.phi",
                            printQ                   = T)
# C version:
tesf$adj.nodes
symbolic_conditional_energy(config                   = as.matrix(configs[jridx,]), 
                            condition_element_number = 2, 
                            node_par                 = theta.pars$node_par, 
                            edge_par                 = theta.pars$edge_par, 
                            edge_mat                 = rmod$edges,
                            adj_nodes                = rmod$adj.nodes
                            )
as.matrix(configs[jridx,])
configs[jridx,]
t(as.matrix(configs[jridx,]))
class(tesf$node_par)
tesf$node.par
theta.pars
rmod$adj.nodes[[1]]
rmod$adj.nodes
