#'@export
ahistory2genealogy <- function(ahistory){
  nsteps <- length(ahistory)
  nparticles <- length(ahistory[[1]])
  genealogy.df <- data.frame()
  for (istep in 2:nsteps){
    genealogy.df <- rbind(genealogy.df, data.frame(step = rep(istep-1, nparticles),
                                                   currentgen = ahistory[[istep]], nextgen = 1:nparticles))
  }
  bs <- matrix(nrow = nparticles, ncol = nsteps)
  bs[,nsteps] <- 1:nparticles
  for (iparticle in 1:nparticles){
    for (istep in (nsteps-1):1){
      bs[iparticle,istep] <- ahistory[[istep+1]][bs[iparticle,istep+1]]
    }
  }
  bs.df <- reshape2::melt(bs) %>% rename(particle=Var1, time=Var2, b=value)
  dendro.df <- data.frame()
  for (istep in 2:nsteps){
    from = bs.df %>% filter(time == istep-1) %>% pull(b)
    to = bs.df %>% filter(time == istep) %>% pull(b)
    dendro.df <- rbind(dendro.df,
                       data.frame(from = paste0(from, "t", istep-1), to = paste0(to, "t", istep),
                                  istep = istep-1))
  }
  return(list(genealogy = genealogy.df, lineage = bs.df, dendro = dendro.df))
}
