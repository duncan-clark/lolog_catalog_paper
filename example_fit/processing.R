# This script processed the .gexf file for the network,
# and adds covariates from the covariates from the covariate.csv file
# additional objects are created for the purpose of modelling edge attributes.

#Required Packages
library("rgexf")
library("statnet")
library("Matrix")

factor_convert <- function(factor) {
  levels(factor)[match(factor, levels(factor))]
}

dat <- read.gexf("data.gexf")

att <- read.csv("data_covariates.csv")
att <- att[, c(1, 2, 4, 5, 6, 7)]
att <- cbind(att, (att[, 5] + att[, 6]))

#switch order so that 53 is highest
att$Lokationposition1highest53lowest <-
  54 - att$Lokationposition1highest53lowest
names(att)[match("Lokationposition1highest53lowest", names(att))] <-
  "rank"

att[[2]] <-
  factor_convert(att[[2]])

#list names of sweets givers or repeater
att$Label[as.logical(att[, 7])]
#tmp <- sna::degree(net,cmode = "indegree")
#View(cbind(tmp,att$Label))

#Need to make adjacency matrix since the using the edge list results in some missing nodes:
dat$adj <-
  matrix(rep(0, length(dat$nodes[, 1]) ** 2), nrow = length(dat$nodes[, 1]))
for (i in 1:length(dat$nodes[, 1])) {
  for (j in 1:length(dat$nodes[, 1][-i])) {
    if (sum((dat$edges$source == i) * (dat$edges$target == j)) == 1) {
      dat$adj[i, j] <- 1
    }
  }
}

net <-
  network(
    dat$adj,
    vertex.attr = cbind(att, att, att) ,
    vertex.attrnames = list(
      "id",
      "names",
      "rank",
      "handicapped",
      "repeater",
      "sweetgiver",
      "repeater_or_sweets",
      "id_in",
      "names_in",
      "rank_in",
      "handicapped_in",
      "repeater_in",
      "sweetgiver_in",
      "repeater_or_sweets_in",
      "id_out",
      "names_out",
      "rank_out",
      "handicapped_out",
      "repeater_out",
      "sweetgiver_out",
      "repeater_or_sweets_out"
    )
  )

#add in duplicated nodal covariates so can cope with in and out degree terms
edges <-
  expand.grid(seq(1, length(dat$nodes[, 1])), seq(1, length(dat$nodes[, 1])))
#var1 is source, var2 is destination
edges$rank1 <- att$rank[match(edges$Var1, att$Id)]
edges$rank2 <- att$rank[match(edges$Var2, att$Id)]
edges$uprank <- edges$rank2 > edges$rank1
edges$rankdiff <- abs(edges$rank2 - edges$rank1)

edges_adj_uprank <-
  matrix(rep(0, length(dat$nodes[, 1]) ** 2), nrow = length(dat$nodes[, 1]))
for (i in 1:length(dat$nodes[, 1])) {
  for (j in (1:length(dat$nodes[, 1]))[-i]) {
    tmp <- which(edges$Var1 == i & edges$Var2 == j)
    if (edges$uprank[tmp] == TRUE) {
      edges_adj_uprank[i, j] <- 1
    }
  }
}

edges_adj_rankdiff <-
  matrix(rep(0, length(dat$nodes[, 1]) ** 2), nrow = length(dat$nodes[, 1]))
for (i in 1:length(dat$nodes[, 1])) {
  for (j in (1:length(dat$nodes[, 1]))[-i]) {
    tmp <- which(edges$Var1 == i & edges$Var2 == j)
    edges_adj_rankdiff[i, j] <- edges$rankdiff[tmp]
  }
}

set.edge.value(net, attrname = "uprank", value = edges_adj_uprank)
set.edge.value(net, attrname = "rankdiff", value = edges_adj_uprank)

plot_func_repeater <- function(net, arrows, title) {
  colours = get.vertex.attribute(net, "repeater")
  colours[colours == 0] <- "Deep Sky Blue"
  colours[colours == 1] <- "Dark Orange"
  gplot(
    net,
    vertex.col = colours,
    vertex.cex = (sna::degree(net, cmode = "indegree") / 10 + 0.5),
    arrowhead.cex = 100 / network.edgecount(net),
    usearrows = arrows,
    mode = "kamadakawai",
    main = title,
    jitter = TRUE,
    displayisolates = TRUE
  )
  pdf(file = paste(title, ".pdf", sep = ""))
  gplot(
    net,
    vertex.col = colours,
    vertex.cex = (sna::degree(net, cmode = "indegree") / 10 + 0.5),
    arrowhead.cex = 100 / network.edgecount(net),
    usearrows = arrows,
    mode = "kamadakawai",
    main = title,
    jitter = TRUE,
    displayisolates = TRUE
  )
  dev.off()
}

plot_func_repeater_or_sweets <- function(net, arrows, title) {
  colours = get.vertex.attribute(net, "repeater_or_sweets")
  colours[colours == 0] <- "Deep Sky Blue"
  colours[colours == 1] <- "Dark Orange"
  
  gplot(
    net,
    vertex.col = colours,
    vertex.cex = (sna::degree(net, cmode = "indegree") / 10 + 0.5),
    edge.lwd =  50 / network.edgecount(net),
    arrowhead.cex = 50 / network.edgecount(net),
    usearrows = arrows,
    mode = "kamadakawai",
    main = title,
    jitter = TRUE,
    displayisolates = TRUE
  )
  pdf(file = paste(title, ".pdf", sep = ""))
  gplot(
    net,
    vertex.col = colours,
    vertex.cex = (sna::degree(net, cmode = "indegree") / 10 + 0.5),
    edge.lwd =  50 / network.edgecount(net),
    arrowhead.cex = 50 / network.edgecount(net),
    usearrows = arrows,
    mode = "kamadakawai",
    main = title,
    jitter = TRUE,
    displayisolates = TRUE
  )
  dev.off()
}

plot_func_handicap <- function(net, arrows, title) {
  colours = get.vertex.attribute(net, "handicapped")
  colours[colours == 0] <- "Deep Sky Blue"
  colours[colours == 1] <- "Dark Orange"
  gplot(
    net,
    vertex.col = colours,
    vertex.cex = (sna::degree(net, cmode = "indegree") / 10 + 0.5),
    arrowhead.cex = 100 / network.edgecount(net),
    usearrows = arrows,
    mode = "kamadakawai",
    main = title,
    jitter = TRUE,
    displayisolates = TRUE
  )
  pdf(file = paste(title, ".pdf", sep = ""))
  gplot(
    net,
    vertex.col = colours,
    vertex.cex = (sna::degree(net, cmode = "indegree") / 10 + 0.5),
    arrowhead.cex = 100 / network.edgecount(net),
    usearrows = arrows,
    main = title,
    jitter = TRUE,
    displayisolates = TRUE
  )
  dev.off()
}

#plotting looks  reasonable
plot_func_repeater(net,
                   TRUE,
                   "Network of 19th Century German Schoolboys with repeaters orange")
plot_func_repeater_or_sweets(
  net,
  TRUE,
  "Network of 19th Century German Schoolboys with repeaters_or_sweetgivers orange"
)
plot_func_handicap(net,
                   TRUE,
                   "Network of 19th Century German Schoolboys with handicapped orange")

save(net, edges_adj_uprank, edges_adj_rankdiff, file = "processed_data.RData")
