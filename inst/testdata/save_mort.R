setwd(dirname( rstudioapi::getSourceEditorContext()$path ))

load('~/Library/CloudStorage/Box-Box/macbook/Warping/mort/1983-2013Female.RData')
save(yearCom, qf, prob, optns, country,
     file = 'mort_female_1983_2013.Rda')
