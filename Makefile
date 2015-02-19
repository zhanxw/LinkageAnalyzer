all:

VER = $(shell grep -m1 Version LinkageAnalyzer/DESCRIPTION |cut -f2- -d' ')

doc:
	Rscript genDoc.R

dist:
	R CMD build LinkageAnalyzer


inst: dist
	R CMD INSTALL LinkageAnalyzer_$(VER).tar.gz

purcell: dist
	rsync LinkageAnalyzer_$(VER).tar.gz zhanxw@purcell.swmed.edu:~/forDavid
