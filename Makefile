# http://www.gnu.org/software/make/manual/make.html

PDF := talk.pdf
TEX := talk.tex
PLOTS := \
		plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt-0_lam-0.1_kap-0.1.pdf \
		plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt-0_lam-0.3_kap-0.3.pdf \
		plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt--2_lam-0.1_kap-0.1.pdf \
		plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt--2_lam-0.3_kap-0.3.pdf \
		plots/SOFTSUSY/Mh_MS_TB-20_Xt--sqrt6.pdf \
		plots/SOFTSUSY/SS_TB-20_Xt--sqrt6.pdf \
		plots/SOFTSUSY/SS_TB-20_Xt--sqrt6_individual.pdf \
		plots/SOFTSUSY/HSSUSY_TB-20_Xt--sqrt6_individual.pdf \
		plots/SOFTSUSY/Mh_2L_vs_3L_MS_TB-20_Xt--sqrt6.pdf \
		plots/Mh3L/scan_Mh_MS_TB-20_Xt--sqrt6_uncertainty_Qpole.pdf \
		plots/FlexibleEFTHiggs-2L/Mh_MS_TB-20_Xt-0.pdf \
		plots/FlexibleEFTHiggs-2L/Mh_MS_TB-20_Xt--sqrt6.pdf \
		plots/HSSUSY-3L/scan_Mh_MS_TB-10_Xt--sqrt6.pdf \
		plots/HSSUSY-3L/scan_Mh_MS_TB-10_Xt--sqrt6_diff.pdf

.PHONY: all clean

all: $(PDF)

plots/SOFTSUSY/SS_TB-20_Xt--sqrt6.pdf: plots/SOFTSUSY/plot_SS.sh plots/SOFTSUSY/*.dat plots/SOFTSUSY/HSSUSY-3L/*.dat
	cd plots/SOFTSUSY/ && ./plot_SS.sh

plots/SOFTSUSY/SS_TB-20_Xt--sqrt6_individual.pdf: plots/SOFTSUSY/plot_SS_individual.sh plots/SOFTSUSY/*.dat plots/SOFTSUSY/HSSUSY-3L/*.dat
	cd plots/SOFTSUSY/ && ./plot_SS_individual.sh

plots/SOFTSUSY/Mh_MS_TB-20_Xt--sqrt6.pdf: plots/SOFTSUSY/plot_SS.sh plots/SOFTSUSY/*.dat
	cd plots/SOFTSUSY/ && ./plot_SS.sh

plots/SOFTSUSY/HSSUSY_TB-20_Xt--sqrt6_individual.pdf: plots/SOFTSUSY/plot_HSSUSY_individual.sh plots/SOFTSUSY/*.dat
	cd plots/SOFTSUSY/ && ./plot_HSSUSY_individual.sh

plots/SOFTSUSY/Mh_2L_vs_3L_MS_TB-20_Xt--sqrt6.pdf: plots/SOFTSUSY/plot_SS.sh plots/SOFTSUSY/*.dat
	cd plots/SOFTSUSY/ && ./plots/SOFTSUSY/plot_SS.sh

plots/Mh3L/scan_Mh_MS_TB-20_Xt--sqrt6_uncertainty_Qpole.pdf: plots/Mh3L/plot_scale_uncertainty.sh plots/Mh3L/*.dat
	cd plots/Mh3L/ && ./plot_scale_uncertainty.sh

plots/FlexibleEFTHiggs-2L/Mh_MS_TB-20_Xt-0.pdf: plots/FlexibleEFTHiggs-2L/plot.sh plots/FlexibleEFTHiggs-2L/*.dat
	cd plots/FlexibleEFTHiggs-2L && ./plot.sh

plots/FlexibleEFTHiggs-2L/Mh_MS_TB-20_Xt--sqrt6.pdf: plots/FlexibleEFTHiggs-2L/plot.sh plots/FlexibleEFTHiggs-2L/*.dat
	cd plots/FlexibleEFTHiggs-2L && ./plot.sh

plots/HSSUSY-3L/scan_Mh_MS_TB-10_Xt--sqrt6.pdf: plots/HSSUSY-3L/plot_Mh_MS_Xt--sqrt6.sh plots/HSSUSY-3L/*.dat
	cd plots/HSSUSY-3L/ && ./plot_Mh_MS_Xt--sqrt6.sh

plots/HSSUSY-3L/scan_Mh_MS_TB-10_Xt--sqrt6_diff.pdf: plots/HSSUSY-3L/plot_Mh_MS_Xt--sqrt6_diff.sh plots/HSSUSY-3L/*.dat
	cd plots/HSSUSY-3L/ && ./plot_Mh_MS_Xt--sqrt6_diff.sh

plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt-0_lam-0.1_kap-0.1.pdf \
plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt-0_lam-0.3_kap-0.3.pdf \
plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt--2_lam-0.1_kap-0.1.pdf \
plots/NMSSMEFTHiggs/DMh_MS_TB-5_Xt--2_lam-0.3_kap-0.3.pdf: \
plots/NMSSMEFTHiggs/plot_DMh_MS.sh \
plots/NMSSMEFTHiggs/*.dat
	$<

$(PDF): $(TEX) $(PLOTS)
	pdflatex $<
	cd Feynman && ./makeall.sh
	pdflatex $<

clean:
	rm -f *~ *.bak *.aux *.log *.toc *.bbl *.blg *.nav *.out *.snm *.backup
	rm -f plots/*.aux plots/*.log
	rm -f $(PLOTS)

distclean: clean
	rm -f $(PDF)
