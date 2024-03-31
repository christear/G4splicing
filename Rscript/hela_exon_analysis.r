#### code for regenerating figures based on Hela cassete exons sequence analysis, including G content, G4 motif frequency, G4-hunter score analysis
#### realted figures: Figure 1E-1H, Figure S1G
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
cat('Hela splicing data analysis...\n')
# loading the data
#data_path='data_for_uploading'
#save(list = c('helaexons','gcontentl','subg4l'),file = 'hela.exon.sequence.features.Rdata')
load(paste0(data_path,'/hela.exon.sequence.features.Rdata'))
helatypes = c('exon','dnsexon_ss3','ss3','ss5','upsexon_ss5')
# G4-hunter score
ss5g4hs = as.data.frame(fread(paste0(data_path,'/hela.ss5.g4h.score.tab'),head = T))
ss3g4hs = as.data.frame(fread(paste0(data_path,'/hela.ss3.g4h.score.tab'),head = T))
ss5g4hso = ss5g4hs[order(gsub('\\(.*','',ss5g4hs[,1])),]
ss3g4hso = ss3g4hs[order(gsub('\\(.*','',ss3g4hs[,1])),]
# 
if (!(file.exists('raw_pdf'))){
    system('mkdir -p raw_pdf')
}

###
cat('processing cassete exons ...\n')
helaexons$DPSI_32vs40 = apply(helaexons[,grep('PSI_WT32_',colnames(helaexons))],1,mean) - apply(helaexons[,grep('PSI_WT40_',colnames(helaexons))],1,mean)
subexons = helaexons[helaexons[,7] - helaexons[,6] > lencut,]
subepsiave = as.data.frame(cbind(apply(subexons[,14:17],1,mean),apply(subexons[,18:21],1,mean),apply(subexons[,24:27],1,mean),subexons[,c(12,22)]))
colnames(subepsiave)[1:3] = c('avePSI_32','avePSI_37','avePSI_40')
subepsiave$exon_coord = paste('chr',subexons[,4],':',subexons[,6],':',subexons[,7],sep = '')
subepsiave$id = paste(subexons[,1],subexons[,3],sep = ':')
subepsiave$DPSI_WT_32vs37 = subexons$DPSI_WT_32vs37
subepsiave$DPSI_WT_37vs40 = subexons$DPSI_WT_37vs40
# get cassete exon id
ces = subexons[(subepsiave[,1] > 0.05 & subepsiave[,1] < 0.95) | (subepsiave[,2] > 0.05 & subepsiave[,2] < 0.95) | (subepsiave[,3] > 0.05 & subepsiave[,3] < 0.95),]
ceid = paste(ces[,1],ces[,3],sep = ':')
# classify cassete exons into three categories
# related to Figure 1E
cre1 = subexons[subexons$PvalueWT_32vs37 < pcut & subexons$DPSI_WT_32vs37 < dcut,]
cre2 = subexons[subexons$PvalueWT_37vs40 < pcut & subexons$DPSI_WT_37vs40 < dcut,]
cie1 = subexons[subexons$PvalueWT_32vs37 < pcut & subexons$DPSI_WT_32vs37 > -dcut,]
cie2 = subexons[subexons$PvalueWT_37vs40 < pcut & subexons$DPSI_WT_37vs40 > -dcut,]
cie1id = paste(cie1[,1],cie1[,3],sep = ':')
cie2id = paste(cie2[,1],cie2[,3],sep = ':')
cre1id = paste(cre1[,1],cre1[,3],sep = ':')
cre2id = paste(cre2[,1],cre2[,3],sep = ':')
bge1 = subexons[subexons$PvalueWT_32vs37 > 0.1 & abs(subexons$DPSI_WT_32vs37) < abs(dcut),]
bge2 = subexons[subexons$PvalueWT_37vs40 > 0.1 & abs(subexons$DPSI_WT_37vs40) < abs(dcut),]
bge1id = intersect(paste(bge1[,1],bge1[,3],sep = ':'),ceid)
bge2id = intersect(paste(bge2[,1],bge2[,3],sep = ':'),ceid)
bgeid = intersect(bge1id,bge2id)
# G content analysis
cat('G content analysis ...\n')
subgcontentl = lapply(gcontentl,function(x) cbind(as.data.frame(gsub("\\(.*","",x[,1])),x[,2:3]))
cre1gcontentl = lapply(subgcontentl,function(x) x[x[,1] %in% cre1id,])
cre2gcontentl = lapply(subgcontentl,function(x) x[x[,1] %in% cre2id,])
cie1gcontentl = lapply(subgcontentl,function(x) x[x[,1] %in% cie1id,])
cie2gcontentl = lapply(subgcontentl,function(x) x[x[,1] %in% cie2id,])
bgegcontentl = lapply(subgcontentl,function(x) x[x[,1] %in% bgeid,])
#
hela1da = data.frame(content = c(unlist(lapply(cre1gcontentl,function(x) x[,2])),unlist(lapply(bgegcontentl,function(x) x[,2])),unlist(lapply(cie1gcontentl,function(x) x[,2]))),region = factor(c(rep(helatypes,unlist(lapply(cre1gcontentl,nrow))),rep(helatypes,unlist(lapply(bgegcontentl,nrow))),rep(helatypes,unlist(lapply(cie1gcontentl,nrow)))),levels = c('upsexon_ss5','ss3','exon','ss5','dnsexon_ss3')),condition = rep(c('CRE','NS','CIE'),c(length(unlist(cre1gcontentl))/3,length(unlist(bgegcontentl))/3,length(unlist(cie1gcontentl))/3)))
hela1da$condition = factor(hela1da$condition,levels = c('CIE','NS','CRE'))

hela2da = data.frame(content = c(unlist(lapply(cre2gcontentl,function(x) x[,2])),unlist(lapply(bgegcontentl,function(x) x[,2])),unlist(lapply(cie2gcontentl,function(x) x[,2]))),region = factor(c(rep(helatypes,unlist(lapply(cre2gcontentl,nrow))),rep(helatypes,unlist(lapply(bgegcontentl,nrow))),rep(helatypes,unlist(lapply(cie2gcontentl,nrow)))),levels = c('upsexon_ss5','ss3','exon','ss5','dnsexon_ss3')),condition = rep(c('CRE','NS','CIE'),c(length(unlist(cre2gcontentl))/3,length(unlist(bgegcontentl))/3,length(unlist(cie2gcontentl))/3)))
hela2da$condition = factor(hela2da$condition,levels = c('CIE','NS','CRE'))
# related to Figure 1F and Figure S1G
pdf('./raw_pdf/Hela.cassete.exons.G.content.pdf')
ggboxplot(hela1da[hela1da$region != 'exon',],x = 'region',y = 'content',fill = 'condition',outlier.shape = NA,palette = 'npg') -> fp1
print(fp1)
ggboxplot(hela2da[hela2da$region != 'exon',],x = 'region',y = 'content',fill = 'condition',outlier.shape = NA,palette = 'npg') -> fp2
print(fp2)
dev.off()

# G4 motif analysis
cat('G4 motif analysis ...\n')
cre1g4l = lapply(subg4l,function(x) x[x[,1] %in% cre1id,])
bgeg4l = lapply(subg4l,function(x) x[x[,1] %in% bgeid,])
cie1g4l = lapply(subg4l,function(x) x[x[,1] %in% cie1id,])
g4prop = rbind(unlist(lapply(cre1g4l,nrow))/length(cre1id),unlist(lapply(bgeg4l,nrow))/length(bgeid),unlist(lapply(cie1g4l,nrow))/length(cie1id)) * 100
rownames(g4prop) = c('CRE','NS','CIE')
colnames(g4prop) = helatypes
g4propda = data.frame(proportion = c(g4prop),region = rep(colnames(g4prop),each = nrow(g4prop)),condition = factor(rep(rownames(g4prop),ncol(g4prop)),levels = c('CIE','NS','CRE')))
g4propda$region = factor(g4propda$region,levels = c('upsexon_ss5','ss3','exon','ss5','dnsexon_ss3'))
# related to Figure 1G
pdf('./raw_pdf/Hela.cassete.exon.G4prop.pdf')
ggbarplot(data = g4propda[g4propda$region != 'exon',],x = 'region',y = 'proportion',fill = 'condition',position = position_dodge(0.8),palette = 'npg') + ylim(c(0,20)) -> fp3
print(fp3)
dev.off()

# G4 hunter analysis
cat('G4hunter score analysis ...\n')
bgess3hc = ss3g4hs[gsub('\\(.*','',ss3g4hs[,1]) %in% bgeid,-1]
bgess5hc = ss5g4hs[gsub('\\(.*','',ss5g4hs[,1]) %in% bgeid,-1]
cre1ss3hc = ss3g4hs[gsub('\\(.*','',ss3g4hs[,1]) %in% cre1id,-1]
cre1ss5hc = ss5g4hs[gsub('\\(.*','',ss5g4hs[,1]) %in% cre1id,-1]
cie1ss3hc = ss3g4hs[gsub('\\(.*','',ss3g4hs[,1]) %in% cie1id,-1]
cie1ss5hc = ss5g4hs[gsub('\\(.*','',ss5g4hs[,1]) %in% cie1id,-1]
xseq = seq(-200,175)
cols = c('#00FFFF','#FF0000','#127D09')
lwd = 1.5
# related to Figure 1G
pdf('./raw_pdf/Hela.sss.g4hscore.pdf')
plot(xseq,apply(bgess3hc,2,mean),type = 'l',lwd = lwd,ylim = c(-0.4,0.4),ylab = 'G4Hunter score',xlab = 'Relative position to SS3',col = cols[1])
lines(xseq,apply(cre1ss3hc,2,mean),col = cols[2],lwd = lwd)
lines(xseq,apply(cie1ss3hc,2,mean),col = cols[3],lwd = lwd)
abline(h = 0,lty = 2)
abline(v = 0,lty = 2)
legend('topright',legend = c('NS','CRE','CIE'),col = cols,lty = 1)
plot(xseq,apply(bgess5hc,2,mean),type = 'l',lwd = lwd,ylim = c(-0.4,0.4),ylab = 'G4Hunter score',xlab = 'Relative position to SS5',col = cols[1])
lines(xseq,apply(cre1ss5hc,2,mean),col = cols[2],lwd = lwd)
lines(xseq,apply(cie1ss5hc,2,mean),col = cols[3],lwd = lwd)
abline(h = 0,lty = 2)
abline(v = 0,lty = 2)
legend('topright',legend = c('NS','CRE','CIE'),col = cols,lty = 1)
dev.off()





