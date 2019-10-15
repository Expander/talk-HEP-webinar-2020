# http://www.gnu.org/software/make/manual/make.html

OUTFILENAME := talk.pdf
PLOTS       := \
		plots/SOFTSUSY/Mh_MS_TB-20_Xt--sqrt6.pdf \
		plots/SOFTSUSY/SS_TB-20_Xt--sqrt6.pdf

TEXDIRS     := $(PLOTSDIR)
BIBTEX      := bibtex

.PHONY: all clean

all: $(OUTFILENAME)

plots/SOFTSUSY/SS_TB-20_Xt--sqrt6.pdf: plots/SOFTSUSY/plot_SS.sh plots/SOFTSUSY/*.dat plots/SOFTSUSY/HSSUSY-3L/*.dat
	cd plots/SOFTSUSY/ && ./plot_SS.sh

plots/SOFTSUSY/Mh_MS_TB-20_Xt--sqrt6.pdf: plots/SOFTSUSY/plot_SS.sh plots/SOFTSUSY/*.dat
	cd plots/SOFTSUSY/ && ./plot_SS.sh

%.pdf: %.tex $(PLOTS)
	pdflatex $<
	cd Feynman && ./makeall.sh
	pdflatex $<

clean:
	rm -f *~ *.bak *.aux *.log *.toc *.bbl *.blg *.nav *.out *.snm *.backup
	rm -f plots/*.aux plots/*.log
	rm -f $(PLOTS)

distclean: clean
	rm -f $(OUTFILENAME)
