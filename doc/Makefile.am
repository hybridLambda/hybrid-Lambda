#docfiles = manual.tex 
if HAVE_PDFLATEX
doc_DATA = manual.pdf
MANNAME = manual
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf
DOCDIR=$(abs_top_builddir)/doc
CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANAUX) 
#$(MANNAME).bbl $(MANNAME).blg 

$(MANPDF): $(srcdir)/$(MANTEXSRC)
	$(PDFLATEX) $(MANTEXSRC)
#	$(BIBTEX) $(MANAUX)
#	$(PDFLATEX) $(MANTEXSRC)
	$(PDFLATEX) $(MANTEXSRC)
endif
#dist_doc_DATA = $(docfiles)

EXTRA_DIST = manual.pdf 

#SUBDIRS =

net1fig: 
	rm net1-figure*
	latex --shell-escape net1.tex
	convert -density 300 net1-figure0.eps -resize 1024x1024 net1-figure0.png
	convert -density 300 net1-figure1.eps -resize 1024x1024 net1-figure1.png
	convert -density 300 net1-figure2.eps -resize 1024x1024 net1-figure2.png

#if [ -a net1-figure0.md5 ] ; then rm net1-figure0.md5 ; fi
