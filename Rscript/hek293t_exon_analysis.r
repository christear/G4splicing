#### code for regenerating figures based on HEK293T cassete exons sequence analysis, and Hela-HEK293T shared temperature-related exons
#### realted figures: Figure 1A and 1D, Figure S1B-S1E, Figure 2A and 2B
### loading the required library
library(ggpubr)
library(gplots)
library(data.table)

### setting default thresholds
mins = 30
lencut = 20
pcut = 0.05
dcut = -0.1

###
cat('HEK293T splicing data analysis...\n')
# loading the data
#data_path='data'
#save(list = c('hekavepsi','hekgcontentl','subg4l'),file = 'hela.exon.sequence.features.Rdata')
load(paste0(data_path,'/hek293t.exon.sequence.features.Rdata'))
hektypes = c('dnsexon_ss3','exon','ss3','ss5','upsexon_ss5')
# G4-hunter score: id splitted by ::
hekss5g4hs = as.data.frame(fread(paste0(data_path,'/hek293t.ss5.g4h.score.tab'),head = T))
hekss3g4hs = as.data.frame(fread(paste0(data_path,'/hek293t.ss3.g4h.score.tab'),head = T))
hekss5g4hso = hekss5g4hs[order(gsub('::.*','',hekss5g4hs[,1])),]
hekss3g4hso = hekss3g4hs[order(gsub('::.*','',hekss3g4hs[,1])),]
#
if (!(file.exists('raw_pdf'))){
    system('mkdir -p raw_pdf')
}

# cassete exon
cat('processing cassete exons ...\n')
tmp1 = apply(hekavepsi[,9:10],1,function(x) min(x[!(is.na(x))]))
tmp2 = apply(hekavepsi[,9:10],1,function(x) max(x[!(is.na(x))]))
hekexonsid = paste(hekavepsi[,1],hekavepsi[,2],sep = ':')
hekceid = hekexonsid[tmp1 < 0.95 & tmp2 > 0.05]

# splicing changes at different temperatures
oa33_dpsi = hekavepsi$aveOA_33 - hekavepsi$aveDMSO_33
oa39_dpsi = hekavepsi$aveOA_39 - hekavepsi$aveDMSO_39_2
tg35_dpsi = hekavepsi$aveTG003_35 - hekavepsi$aveDMSO_35
tg39_dpsi = hekavepsi$aveTG003_39 - hekavepsi$aveDMSO_39
cls2 = get_palette('npg',2)
lwd = 2
# related to Figure 1A and Figure S1B
pdf('raw_pdf/HEK293T.OA.TG03.dpsi.pdf')
par(mfrow = c(2,2))
plot(ecdf(abs(tg35_dpsi)),col = cls2[2],xlim = c(0,0.7),lwd = lwd,main = 'TG003 vs DMSO',xlab = '|Delta PSI|',ylab = 'Cumulative fraction')
lines(ecdf(abs(tg39_dpsi)),col = cls2[1],lwd = lwd)
legend('bottomright',legend = c('35','39'),col = rev(cls2),lty = 1)
abline(h = 0.5,lty = 2)
abline(v = 0,lty = 2)
plot(ecdf(abs(oa33_dpsi)),col = cls2[2],xlim = c(0,0.5),lwd = lwd,main = 'OA vs DMSO',xlab = '|Delta PSI|',ylab = 'Cumulative fraction')
lines(ecdf(abs(oa39_dpsi)),col = cls2[1],lwd = lwd)
legend('bottomright',legend = c('33','39'),col = rev(cls2),lty = 1)
abline(h = 0.5,lty = 2)
abline(v = 0,lty = 2)
dev.off()

# related to Figure 1D and S1C
dmso_creid = hekexonsid[hekavepsi$Pval_DMSO < pcut & hekavepsi$aveDMSO_39 - hekavepsi$aveDMSO_35 >= abs(dcut)]
dmso_cieid = hekexonsid[hekavepsi$Pval_DMSO < pcut & hekavepsi$aveDMSO_39 - hekavepsi$aveDMSO_35 <= dcut]
dmso_bgeid = hekexonsid[hekavepsi$Pval_DMSO > 0.1 & abs(hekavepsi$aveDMSO_39 - hekavepsi$aveDMSO_35) < abs(dcut)]
dmso_bgeid = intersect(dmso_bgeid,hekceid)
tg03_creid = hekexonsid[hekavepsi$Pval_TG003 < pcut & hekavepsi$aveTG003_39 - hekavepsi$aveTG003_35 >= abs(dcut)]
tg03_cieid = hekexonsid[hekavepsi$Pval_TG003 < pcut & hekavepsi$aveTG003_39 - hekavepsi$aveTG003_35 <= dcut]
tg03_bgeid = hekexonsid[hekavepsi$Pval_TG003 > 0.1 & abs(hekavepsi$aveTG003_39 - hekavepsi$aveTG003_35) < abs(dcut)]
tg03_bgeid = intersect(tg03_bgeid,hekceid)

### G content
cat('G content analysis ...\n')
hekgcontentl2 = lapply(hekgcontentl,function(x) cbind(as.data.frame(gsub("::.*","",x[,1])),x[,2:3]))
# depend on CLK
depecregcontentl = lapply(hekgcontentl2,function(x) x[x[,1] %in% setdiff(dmso_creid,tg03_creid),])
depeciegcontentl = lapply(hekgcontentl2,function(x) x[x[,1] %in% setdiff(dmso_cieid,tg03_cieid),])
# independent of CLK
indecregcontentl = lapply(hekgcontentl2,function(x) x[x[,1] %in% intersect(tg03_creid,dmso_creid),])
indeciegcontentl = lapply(hekgcontentl2,function(x) x[x[,1] %in% intersect(tg03_cieid,dmso_cieid),])
### using the DMSO bgeid as NS for CLK-dep and CLK-indep
#combbgegcontentl = lapply(hekgcontentl2,function(x) x[x[,1] %in% unique(dmso_bgeid),])
combbgegcontentl1 = lapply(hekgcontentl2,function(x) x[x[,1] %in% setdiff(dmso_bgeid,tg03_bgeid),])
combbgegcontentl2 = lapply(hekgcontentl2,function(x) x[x[,1] %in% intersect(dmso_bgeid,tg03_bgeid),])
#
depegcda = data.frame(content = c(unlist(lapply(depecregcontentl,function(x) x[,2])),unlist(lapply(combbgegcontentl1,function(x) x[,2])),unlist(lapply(depeciegcontentl,function(x) x[,2]))),region = factor(c(rep(hektypes,unlist(lapply(depecregcontentl,nrow))),rep(hektypes,unlist(lapply(combbgegcontentl1,nrow))),rep(hektypes,unlist(lapply(depeciegcontentl,nrow)))),levels = c('upsexon_ss5','ss3','exon','ss5','dnsexon_ss3')),condition = rep(c('heat_induced_exon','bg','cold_induced_exon'),c(length(unlist(depecregcontentl))/3,length(unlist(combbgegcontentl1))/3,length(unlist(depeciegcontentl))/3)))
depegcda$condition = factor(depegcda$condition,levels = c('cold_induced_exon','bg','heat_induced_exon'))
indegcda = data.frame(content = c(unlist(lapply(indecregcontentl,function(x) x[,2])),unlist(lapply(combbgegcontentl2,function(x) x[,2])),unlist(lapply(indeciegcontentl,function(x) x[,2]))),region = factor(c(rep(hektypes,unlist(lapply(indecregcontentl,nrow))),rep(hektypes,unlist(lapply(combbgegcontentl2,nrow))),rep(hektypes,unlist(lapply(indeciegcontentl,nrow)))),levels = c('upsexon_ss5','ss3','exon','ss5','dnsexon_ss3')),condition = rep(c('heat_induced_exon','bg','cold_induced_exon'),c(length(unlist(indecregcontentl))/3,length(unlist(combbgegcontentl2))/3,length(unlist(indeciegcontentl))/3)))
indegcda$condition = factor(indegcda$condition,levels = c('cold_induced_exon','bg','heat_induced_exon'))

# related to Figure S1D and S1E
pdf('raw_pdf/HEK293T.cassete.exons.G.content.pdf')
ggboxplot(depegcda[depegcda$region != 'exon',],x = 'region',y = 'content',fill = 'condition',palette = 'npg',outlier.shape = NA)
ggboxplot(indegcda[indegcda$region != 'exon',],x = 'region',y = 'content',fill = 'condition',palette = 'npg',outlier.shape = NA)
dev.off()

#
hekexonsid = paste(hekavepsi[,1],hekavepsi[,2],sep = ':')
hekdpi = cbind(hekavepsi$aveDMSO_39 - hekavepsi$aveDMSO_35,hekavepsi$aveTG003_39 - hekavepsi$aveTG003_35,hekavepsi$aveDMSO_39_2 - hekavepsi$aveDMSO_33,hekavepsi$aveOA_39 - hekavepsi$aveOA_33)
rownames(hekdpi) = hekexonsid
colnames(hekdpi) = c('DMSO_39-35','TG003_39-35','DMSO_39-33','OA_39-33')
hekdpi = hekdpi[order(rownames(hekdpi)),]
# loading Hela data
load(paste0(data_path,"/hela.exon.sequence.features.Rdata"))
# id removing strand inf \([+-])
ss5g4hs = as.data.frame(fread(paste0(data_path,'/hela.ss5.g4h.score.tab'),head = T))
ss3g4hs = as.data.frame(fread(paste0(data_path,'/hela.ss3.g4h.score.tab'),head = T))
ss5g4hso = ss5g4hs[order(gsub('\\(.*','',ss5g4hs[,1])),]
ss3g4hso = ss3g4hs[order(gsub('\\(.*','',ss3g4hs[,1])),]
#
subexons = helaexons[helaexons[,7] - helaexons[,6] > lencut,]
subepsiave = as.data.frame(cbind(apply(subexons[,14:17],1,mean),apply(subexons[,18:21],1,mean),apply(subexons[,24:27],1,mean),subexons[,c(12,22)]))
colnames(subepsiave)[1:3] = c('avePSI_32','avePSI_37','avePSI_40')
subepsiave$exon_coord = paste('chr',subexons[,4],':',subexons[,6],':',subexons[,7],sep = '')
subepsiave$id = paste(subexons[,1],subexons[,3],sep = ':')
subepsiave$DPSI_WT_32vs37 = subexons$DPSI_WT_32vs37
subepsiave$DPSI_WT_37vs40 = subexons$DPSI_WT_37vs40
#
heldpi = cbind(subepsiave$avePSI_37 - subepsiave$avePSI_32,subepsiave$avePSI_40 - subepsiave$avePSI_37)
rownames(heldpi) = subepsiave$id
colnames(heldpi) = c('37-32','40-37')
heldpi = heldpi[order(rownames(heldpi)),]
#
subhekss5g4hso = hekss5g4hso[gsub('::.*','',hekss5g4hso[,1]) %in% rownames(hekdpi),]
subhekss3g4hso = hekss3g4hso[gsub('::.*','',hekss3g4hso[,1]) %in% rownames(hekdpi),]
subhekss5g4hsom = sapply(2:(ncol(subhekss5g4hso) - 25),function(i) apply(subhekss5g4hso[,i:(i+25)],1,max))
subhekss3g4hsom = sapply(2:(ncol(subhekss3g4hso) - 25),function(i) apply(subhekss3g4hso[,i:(i+25)],1,max))
subss5g4hso = ss5g4hso[gsub('\\(.*','',ss5g4hso[,1]) %in% rownames(heldpi),]
subss3g4hso = ss3g4hso[gsub('\\(.*','',ss3g4hso[,1]) %in% rownames(heldpi),]
subss5g4hsom = sapply(2:(ncol(subss5g4hso) - 25),function(i) apply(subss5g4hso[,i:(i+25)],1,max))
subss3g4hsom = sapply(2:(ncol(subss3g4hso) - 25),function(i) apply(subss3g4hso[,i:(i+25)],1,max))

# correlaltion between G4 score and delta PSI
cat('correlaltion between G4 score and delta PSI...\n')
#gcut = 0.5
gcut = 0
hekc3 = sapply(1:ncol(hekdpi), function(j) sapply(2:ncol(subhekss3g4hso),function(i) cor(subhekss3g4hso[subhekss3g4hso[,i] > gcut,i],-hekdpi[subhekss3g4hso[,i]>gcut,j],use = 'p')))
hekc5 = sapply(1:ncol(hekdpi), function(j) sapply(2:ncol(subhekss5g4hso),function(i) cor(subhekss5g4hso[subhekss5g4hso[,i] > gcut,i],-hekdpi[subhekss5g4hso[,i]>gcut,j],use = 'p')))
helc3 = sapply(1:ncol(heldpi), function(j) sapply(2:ncol(subss3g4hso),function(i) cor(subss3g4hso[subss3g4hso[,i] > gcut,i],-heldpi[subss3g4hso[,i] > gcut,j],use = 'p')))
helc5 = sapply(1:ncol(heldpi), function(j) sapply(2:ncol(subss5g4hso),function(i) cor(subss5g4hso[subss5g4hso[,i] > gcut,i],-heldpi[subss5g4hso[,i] > gcut,j],use = 'p')))
# correlation with max G4 score within 25nt window
hekcm3 = sapply(1:ncol(hekdpi), function(j) sapply(2:ncol(subhekss3g4hsom),function(i) cor(subhekss3g4hsom[,i],-hekdpi[,j],use = 'p')))
hekcm5 = sapply(1:ncol(hekdpi), function(j) sapply(2:ncol(subhekss5g4hsom),function(i) cor(subhekss5g4hsom[,i],-hekdpi[,j],use = 'p')))
helcm3 = sapply(1:ncol(heldpi), function(j) sapply(2:ncol(subss3g4hsom),function(i) cor(subss3g4hsom[,i],-heldpi[,j],use = 'p')))
helcm5 = sapply(1:ncol(heldpi), function(j) sapply(2:ncol(subss5g4hsom),function(i) cor(subss5g4hsom[,i],-heldpi[,j],use = 'p')))
#
cls4 = get_palette('npg',7)
lwd = 2
xseq = seq(-200,175)
# related to Figure 2A
pdf('raw_pdf/G4s.dpsi.cor.coe.pdf')
plot(xseq,hekc5[,1],type = 'l',ylim = c(-0.1,0.1),ylab = 'Correlation with delta PSI(SS5)',xlab = 'Relative position to SS5',col = cls4[1],lwd = lwd)
lines(xseq,helc5[,1],col = cls4[3],lwd = lwd)
lines(xseq,-helc5[,2],col = cls4[4],lwd = lwd)
abline(h = 0,lty = 2)
abline(v = 0,lty = 2)
legend('topleft',legend = c('HEK_35-39','Hela_32-37','Hela_40-37'),col = cls4[-2],lty = 1)
dev.off()

###
cat('overlapping exons in HEK293T and Hela ...\n')
hekavepsi$exon_coord = paste(hekavepsi$Chr,hekavepsi$ExonStart,hekavepsi$ExonEnd,sep = ':')
hekavepsi$DPSI_DMSO_35vs39 = hekavepsi$aveDMSO_35 - hekavepsi$aveDMSO_39
oid = intersect(hekavepsi$exon_coord,subepsiave$exon_coord)
# if one exon have multiple rows, select the one with smallest pvalue
get_me2 = function(mat,cn){
    mat = mat[order(mat[,cn]),]
    me = mat[1,]
    me
}
ouniqhek = as.data.frame(t(sapply(1:length(oid),function(i) get_me2(hekavepsi[hekavepsi$exon_coord == oid[i],],7))))
ouniqhela = as.data.frame(t(sapply(1:length(oid),function(i) get_me2(subepsiave[subepsiave$exon_coord == oid[i],],4))))
combuniqres = cbind(ouniqhek[,c(1:6,7,20)],ouniqhela[,6:7],ouniqhela[,c(4,8)])
#write.table(apply(combuniqres,2,as.character),file = 'combined.HEK293.Hela.coldshock.txt',sep = '\t',quote = F)
#combuniqres = read.table('combined.HEK293.Hela.coldshock.txt',head = T)

sigt1 = sigt2 = rep('NS',nrow(combuniqres))
sigt1[combuniqres[,7] < pcut & combuniqres[,8] > 0] = 'UP'
sigt1[combuniqres[,7] < pcut & combuniqres[,8] < 0] = 'DN'
sigt2[combuniqres[,11] < pcut & combuniqres[,12] > 0] = 'UP'
sigt2[combuniqres[,11] < pcut & combuniqres[,12] < 0] = 'DN'
signmat = matrix(table(paste(sigt1,sigt2,sep= ':')),3)
cols7 = get_palette('npg',7)
rownames(signmat) = colnames(signmat) = c('CRE','NS','CIE')
# related to Figure 2B
pdf('./raw_pdf/coldshock.exons.hek.hela.overlap.pdf')
barplot(apply(signmat,2,function(x) x/sum(x)),beside = F,border = F,ylim = c(0,1),col = cols7[3:5])
legend('topleft',legend = c('CRE','NS','CIE'),fill = cols7[3:5])
dev.off()

# add G4 score
hekss3g4hs_max = apply(hekss3g4hs[,176:225],1,max)
hekss5g4hs_max = apply(hekss5g4hs[,176:225],1,max)
names(hekss3g4hs_max) = names(hekss5g4hs_max) = gsub('::.*','',hekss3g4hs[,1])
combid = paste(combuniqres[,1],combuniqres[,2],sep = ':')
subhekss3g4hs_max = hekss3g4hs_max[names(hekss3g4hs_max) %in% combid]
subhekss5g4hs_max = hekss5g4hs_max[names(hekss5g4hs_max) %in% combid]
combuniqres2 = cbind(combuniqres[order(combid),],subhekss3g4hs_max[order(names(subhekss3g4hs_max))],subhekss5g4hs_max[order(names(subhekss5g4hs_max))])
colnames(combuniqres2)[(ncol(combuniqres2) - 1):ncol(combuniqres2)] = c('SS3','SS5')
#write.table(combuniqres2,file = 'combined.HEK293.Hela.coldshock.addG4s.txt',sep = '\t',quote = F)
siguniqres2 = combuniqres2[combuniqres2$PvalueWT_32vs37 < pcut & combuniqres2$DPSI_WT_32vs37 < 0 & combuniqres2$Pval_DMSO < pcut & combuniqres2$DPSI_DMSO_35vs39 < 0,]
#siguniqres2_cut0.1 = combuniqres2[combuniqres2$PvalueWT_32vs37 < pcut & combuniqres2$DPSI_WT_32vs37 < -0.1 & combuniqres2$Pval_DMSO < pcut & combuniqres2$DPSI_DMSO_35vs39 < -0.1,]
sigug4snum = cbind(c(sum(siguniqres2$SS3 < 0.6),sum(siguniqres2$SS3 >= 0.6 & siguniqres2$SS3 < 0.8),sum(siguniqres2$SS3 >= 0.8 & siguniqres2$SS3 < 1),sum(siguniqres2$SS3 >= 1)),c(sum(siguniqres2$SS5 < 0.6),sum(siguniqres2$SS5 >= 0.6 & siguniqres2$SS5 < 0.8),sum(siguniqres2$SS5 >= 0.8 & siguniqres2$SS5 < 1),sum(siguniqres2$SS5 >= 1)),c(sum(apply(siguniqres2[,13:14],1,max) < 0.6),sum(apply(siguniqres2[,13:14],1,max) >= 0.6 & apply(siguniqres2[,13:14],1,max) < 0.8),sum(apply(siguniqres2[,13:14],1,max) >= 0.8 & apply(siguniqres2[,13:14],1,max) < 1),sum(apply(siguniqres2[,13:14],1,max) >= 1)),c(sum(apply(siguniqres2[,13:14],1,min) < 0.6),sum(apply(siguniqres2[,13:14],1,min) >= 0.6 & apply(siguniqres2[,13:14],1,min) < 0.8),sum(apply(siguniqres2[,13:14],1,min) >= 0.8 & apply(siguniqres2[,13:14],1,min) < 1),sum(apply(siguniqres2[,13:14],1,min) >= 1)))
rownames(sigug4snum) = c('<0.6','[0.6,0.8)','[0.8,1)','>=1')
colnames(sigug4snum) = c('SS3','SS5','SS3|SS5','SS3&SS5')

cls5_2 = c('#91003f','#e7298a','#df65b0','#d4b9da')
# related to Figure 2B
pdf('raw_pdf/CRE.ss3.ss5.g4s.pdf')
barplot(apply(sigug4snum,2,function(x) x/sum(x))[,1:3],horiz = T,border = F,col = rev(cls5_2))
legend('top',legend = c('<0.6','[0.6,0.8)','[0.8,1.0)','>1.0'),fill = rev(cls5_2))
dev.off()













