#===========================================================================
# Make dependency file for LaTeX/BibTeX documents.
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

 Target   = paper_I
 TexFiles = 
 FigFiles = 
 BibFile  = sources.bib
 StyFiles = 
#===========================================================================
# LaTeX and BibTeX commands
#===========================================================================

 LaTeX      = pdflatex
 BibTeX     = bibtex
 DVItoPS    = dvips
 DVItoPDF   = dvipdfm

#===========================================================================
# Targets, options, and dependencies
#===========================================================================

all:	${Target}.pdf

# ${Target}.pdf:	${Target}.dvi
# 	${DVItoPDF} ${Target}.dvi
# 	make clean

${Target}.pdf: ${Target}.tex ${Target}.bbl
	${LaTeX} ${Target}
	${LaTeX} ${Target}
	${LaTeX} ${Target}
	make clean

# ${Target}.ps:	${Target}.dvi
# 	${DVItoPS} ${Target}.dvi -o ${Target}.ps
# 	make clean

# ${Target}.dvi:	${TexFiles} ${FigFiles} ${StyFiles} ${Target}.tex ${Target}.bbl 
# 	${LaTeX} ${Target}
# 	${LaTeX} ${Target}

${Target}.bbl:	$(BibFile)
	${LaTeX} ${Target}
	${LaTeX} ${Target}
	${BibTeX} ${Target}


clean:
	rm -f "\#"*"\#" *.dvi *.aux *.log *.bbl *.blg *.toc *.lof *.lot

realclean: clean
	rm -f *~ *.pdf *.ps

#===========================================================================
# Viewing and editing options
#===========================================================================

edit: ${Target}.tex
	${Editor} ${Target}.tex &

view: ${Target}.ps
	${PSview} ${Target}.ps &

read: ${Target}.pdf
	${PDFview} ${Target}.pdf &

xview: ${Target}.dvi
	${DVIview} ${Target}.dvi &

#===========================================================================
# Archive file options
#===========================================================================

# Create tar archive of latex sources (text files only)
archive: 
	make all
	make clean
	rm -Rf $(ArkDir)
	mkdir $(ArkDir)
	cp -R $(ArkFiles) $(ArkDir)
	tar -czf $(ArkDir).tgz $(ArkDir)
	rm -Rf $(ArkDir)

# Create tar archive of all sources
fullarchive: 
	make all
	make clean
	rm -Rf $(ArkDirFull)
	mkdir $(ArkDirFull)
	cp -R $(ArkFiles) $(FigFiles) $(ArkDirFull)
	tar -czf $(ArkDirFull).tgz $(ArkDirFull)
	rm -Rf $(ArkDirFull)


#===========================================================================
# End of Makefile
#===========================================================================
