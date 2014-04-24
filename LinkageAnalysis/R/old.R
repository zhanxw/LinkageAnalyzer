# # G1 female is REF, G1 male is HET # 1) Without G2 data # G2 female (offspring)
# is REF/HET with 50% probability # G2 male is HET (the same as G1 male) #
# REF:HET:VAR=3:4:1 # 2) if G2 mother is REF, REF:HET:VAR=1:1:0 # 3) if G2 mother
# is HET, REF:HET:VAR=1:2:1 # this function forms the matrix for exact test
# generate_exact<-function(mothers,children,observed) { children_trunc=children
# children_trunc[children_trunc>observed]=observed
# exact=matrix(data=0,ncol=length(mothers),nrow=prod(1+children_trunc))
# colnames(exact)=mothers for (r in 1:dim(exact)[2]) # praprogate the added 1 if
# a digit is full { if (r==1) {exact[,r]=rep(0:children_trunc[1])}else
# {exact[,r]=rep(rep(0:children_trunc[r]),each=prod(children_trunc[1:(r-1)]+1))}
# } exact=exact[apply(exact,1,sum)<=observed,] return(exact) } # this function
# finds the probability of producing G3 mice with VAR genotype given maternal
# genotype single_prob<-function(mother_gt) { if (is.null(mother_gt) ||
# is.na(mother_gt)) {mother_gt='unknown'} # corresponding G2 data does not exist
# if (mother_gt=='REF') {prob=0} if (mother_gt=='HET') {prob=1/4} if (!mother_gt
# %in% c('REF','HET')) {prob=1/8} # unknown, FAILED, FALSE return(prob) } # this
# function tests whether each gene is a homozygous lethal with G2 data # the H0
# is a gene is not homozygous lethal
# single_lethal<-function(genotype,genes,phenotype,G2) {
# lethal=rep(1,dim(genotype)[1]) # vector of lethality p values for (i in
# 1:dim(genotype)[1]) { # skip genes with too many invalid values or on chrX
# gt=unlist(genotype[i,]) gt=convert_gt(gt,'recessive') if (sum(!is.na(gt))<3 ||
# genes$chr[i]=='X') {next} observed=sum(gt[!is.na(gt)]==2) # observed number of
# VAR pval=0 children=table(phenotype$mother[!is.na(gt)]) # number of children
# per mother mothers=names(children) # mother names
# exact=generate_exact(mothers,children,observed) for (k in 1:dim(exact)[1]) { if
# (sum(exact[k,])>observed) {next} # larger than the observed # of VAR pvalet=1
# for (mother in mothers) {
# mother_gt=G2[paste(genes$Gene[i],genes$Coord[i]),mother]
# prob=single_prob(mother_gt)
# pvalet=pvalet*dbinom(exact[k,mother],children[mother],prob) } pval=pval+pvalet
# } lethal[i]=pval } return(lethal) } # this function finds the probability to
# produce VAR,HET; HET,VAR and VAR, VAR for a given pair of maternal genotype
# double_prob<-function(mother_gt1,mother_gt2) { if (is.null(mother_gt1) ||
# is.na(mother_gt1) || mother_gt1 %in% c('FALSE','FAILED'))
# {mother_gt1='unknown'} if (is.null(mother_gt2) || is.na(mother_gt2) ||
# mother_gt2 %in% c('FALSE','FAILED')) {mother_gt2='unknown'} if
# (mother_gt1=='REF' && mother_gt2=='REF') {prob=0} # cannot give the genotype we
# desired for sure if (mother_gt1=='HET' && mother_gt2=='HET') {prob=5/16} if
# ((mother_gt1=='HET' && mother_gt2=='REF') || (mother_gt1=='REF' &&
# mother_gt2=='HET')) {prob=1/8} if (mother_gt1=='unknown' &&
# mother_gt2=='unknown') {prob=(1/8+1/8+5/16)/4} if ((mother_gt1=='unknown' &&
# mother_gt2=='REF') || (mother_gt1=='REF' && mother_gt2=='unknown'))
# {prob=(1/8)/2} if ((mother_gt1=='unknown' && mother_gt2=='HET') ||
# (mother_gt1=='HET' && mother_gt2=='unknown')) {prob=(1/8+15/16)/2} return(prob)
# } # this function tests for synthetic lethality of two genes (VAR,HET; HET,VAR;
# and VAR,VAR) double_lethal<-function(data,input,i,j) { if
# (any(input$genes[c(i,j),'chr']=='X')) {return(1)} # don't handle gene on chrX
# for the moment children=table(data$mother) # number of children per mother
# mothers=names(children) # mother names observed=sum((data$gt1+data$gt2)>=3) #
# observed # of VAR,HET; HET,VAR; VAR,VAR pval=0
# exact=generate_exact(mothers,children,observed) # generate the matrix for exact
# tests for (k in 1:dim(exact)[1]) { pvalet=1 for (mother in mothers) {
# mother_gt1=input$G2[paste(input$genes$Gene[i],input$genes$Coordination[i]),mother]
# mother_gt2=input$G2[paste(input$genes$Gene[j],input$genes$Coordination[j]),mother]
# prob=double_prob(mother_gt1,mother_gt2)
# pvalet=pvalet*dbinom(exact[k,mother],children[mother],prob) } pval=pval+pvalet
# } return(pval) } 
