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
	if [ -a net1-figure0.md5 ] ; then rm net1-figure0.md5 ; fi
	latex --shell-escape net1.tex
