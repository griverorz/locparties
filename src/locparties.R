## Description: Location of Galician political parties
## Author: @griverorz
## Date: Sat Mar 22 14:17:58 PDT 2014

setwd("~/Documents/datablog/locparties/")

library(rjags)
library(plyr)
library(reshape2)
library(mvtnorm)
library(ggplot2)

set.seed(201301411)

## Read data in
read_cis <- function(syntaxf, dataf) {
    ## @syntaxf is the location of the syntax file (starting with ES)
    ## @dataf is the location of the data file (starting with DA)
    require(foreign)
    
    ## Add EOL to the last line
    system(paste("sed -i -e '$a\\'", syntaxf, sep = " "))
    
    ## Convert to utf8
    new_syntaxf <- paste("utf8_", basename(syntaxf), sep = "")
    out_syntaxf <- file.path(dirname(syntaxf), new_syntaxf)
    system(paste("iconv -f ISO-8859-1 -t UTF-8", 
                 syntaxf, ">", out_syntaxf, sep = " "))
    
    ## Create output file
    out_dataf <- paste(basename(dataf), ".sav", sep = "")
    exportline <- paste("SAVE OUTFILE='", basename(dataf), ".sav'.", sep = "")
    system(paste("echo", exportline, ">>", out_syntaxf, sep = " "))
    
    ## Replace commas with dots for decimals in datafile
    system(paste("sed -i -e 's/,/./g'", dataf, sep = " "))
    
    ## Run, damnit, run!
    cwd <- getwd()
    setwd(dirname(dataf))
    system(paste("pspp", basename(out_syntaxf), sep = " "), wait = TRUE)
    setwd(cwd)
    
    ## Read SAV file 
    outdf <- read.spss(file.path(dirname(dataf), out_dataf), 
                       to.data.frame = TRUE, use.value.labels = FALSE)
    return(outdf)
}


pgal <- read_cis("./dta/MD2963/ES2963", "./dta/MD2963/DA2963")
codena <- function(x, maxcode = 90) ifelse(x > maxcode, NA, x)

ideo <- codena(pgal$P44)
nacl <- codena(pgal$P35)

ideo.pp <- codena(pgal$P3601)
ideo.psoe <- codena(pgal$P3602)
ideo.age <- codena(pgal$P3603)
ideo.bng <- codena(pgal$P3604)

nacl.pp <- codena(pgal$P4501)
nacl.psoe <- codena(pgal$P4502)
nacl.age <- codena(pgal$P4503)
nacl.bng <- codena(pgal$P4504)

Mloc_nacl <- c(cbind(nacl.pp, nacl.psoe, nacl.age, nacl.bng))
Mloc_ideo <- c(cbind(ideo.pp, ideo.psoe, ideo.age, ideo.bng))
Mloc <- cbind(Mloc_nacl, Mloc_ideo)
party <- c(mapply(rep, c(1,2,3,4), nrow(Mloc)/4))

parties <- 4
dims <- 2

notcomplete <- apply(Mloc, 1, function(x) any(is.na(x)) & !all(is.na(x)))
pos <- vector("list", sum(notcomplete))

t <- 1
for (i in 1:nrow(Mloc)) {
    if (notcomplete[i]) {
        pos[[t]] <- c(i, which(is.na(Mloc[i,]), arr.ind = TRUE))
        t <- t + 1
    }
}
noobs <- do.call(rbind, pos)

dataList <- list(y = Mloc[!notcomplete,],
                 ymiss = Mloc[notcomplete,],
                 noobs = noobs,
                 N = nrow(Mloc[!notcomplete,]),
                 Nnoobs = nrow(noobs),
                 P = party[!notcomplete],
                 Pnoobs = party[notcomplete],
                 Ncluster = parties,
                 muCluster = rep(0, dims),
                 tauCluster = 0.001*diag(dims),
                 lambda = parties,
                 Omega = 0.001*diag(dims))

jmodel <- jags.model("src/locparties.jags", data = dataList,
                     n.chains = 1, n.adapt = 5000)

jsamples <- coda.samples(jmodel, variable.names = c("mu", "Tau"), 
                         n.iter = 5E3, thin = 10)

simulation <- melt(jsamples[[1]])

param_party <- function(simulation, party) {
    mugrep <- paste("mu*.{,3}", party, "\\.$", sep = "")
    Taugrep <- paste("Tau*.{,6}", party, "\\.$", sep = "")
    mu <- grep(mugrep, names(simulation))
    Tau <- grep(Taugrep, names(simulation))
    return(list("mu" = as.matrix(simulation[, mu]), 
                "Tau" = as.matrix(simulation[, Tau])))
}

simulate_party <- function(simulation, party) {
    params <- param_party(simulation, party)
    nn <- nrow(params$mu)
    outlist <- vector("list", 5*nn)
    for (i in 1:(5*nn)) {
        muv <- params$mu[((i - 1) %% nn) + 1, ]
        Tauv <- matrix(params$Tau[((i - 1) %% nn) + 1, ], nrow = 2, ncol = 2)        
        outlist[[i]] <- rmvnorm(1, muv, Tauv)
    }
    return(cbind(do.call(rbind, outlist), party))
}

simvalues <- vector("list", parties)
for (i in 1:parties) {
    simvalues[[i]] <- simulate_party(simulation, i)    
}

simvalues <- do.call(rbind, simvalues)
simvalues <- as.data.frame(simvalues)
names(simvalues) <- c("ideo", "nacl", "party")

ggplot(simvalues, aes(x = ideo, y = nacl, group = party)) + 
    stat_density2d() + 
    scale_x_continuous(limits = c(1, 10)) + 
    scale_y_continuous(limits = c(1, 10)) +
    xlab("Ideology") + 
    ylab("Nationalism")
