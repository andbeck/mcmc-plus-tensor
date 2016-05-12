load('~/Dropbox/####workingprojects/APB_Reger_MVPlast_Qst/Gmatrix_2013/Workspaces/G_IndLHFit_5x5.RData')

source('~/Documents/Repos/mcmc-plus-tensor/derived.stats2.R')
derived.stats2(indFM1, indM1, no.traits = 5)$CoreStats

source('~/Documents/Repos/mcmc-plus-tensor/TestBed/derived.stats2_MS.R')
derived.stats2_MS(indFM1, indM1, no.traits = 5)$CoreStats
