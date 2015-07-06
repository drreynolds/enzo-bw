#===========================================================================
# Make dependency file for PDFLaTeX/BibTeX documents.
#===========================================================================
# To customize this makefile for use with another document, you will
# need to edit the macros in the first section, and to customize this
# makefile for use with other compilation programs, you will need to 
# edit the macros in the second section.
#
# Target   : the primary document name.
# TexFiles : lists actual names of files referenced in the primary document,
#            e.g. with LaTeX's \include command.
# FigFiles : lists names of figure files 
# BibFile  : lists actual names of .bib files referenced in your
#	     document with LaTeX's \bibliography command.
#
# Macros which exceed one line in length are joined with a backslash (\)
# at the end of the line. Macros with more than one value, e.g. multiple
# file names, must be surrounded by parentheses, either () or {}.
#===========================================================================

 Target   = dualfld_UG
 TexFiles = 
 FigFiles = 
 BibFile  = sources.bib
 StyFiles = 

#===========================================================================
# LaTeX and BibTeX commands
#===========================================================================

 LaTeX      = pdflatex
 BibTeX     = bibtex

#===========================================================================
# Targets, options, and dependencies
#===========================================================================

all:	${Target}.pdf

${Target}.pdf:	${TexFiles} ${FigFiles} ${StyFiles} ${Target}.tex ${BibFile}
	${LaTeX} ${Target}
	${LaTeX} ${Target}
	${BibTeX} ${Target}
	${LaTeX} ${Target}
	${LaTeX} ${Target}

clean:
	rm -f "\#"*"\#" *.aux *.log *.bbl *.blg *.toc *.lof *.lot

realclean: clean
	rm -f *~ *.pdf

#===========================================================================
# End of Makefile
#===========================================================================
