SRCS  := $(addprefix lb_walkthroughs/,update_beta.qmd update_L.qmd update_F.qmd update_tau.qmd fit_modular.qmd) findings/2026-03-31-pdac-data-audit/PDAC_data_audit.qmd
PDFS  := $(SRCS:.qmd=.pdf)
HTMLS := $(SRCS:.qmd=.html)

all: $(PDFS) $(HTMLS)

%.pdf: %.qmd
	quarto render $< --to pdf

%.html: %.qmd
	quarto render $< --to html

clean:
	rm -f $(PDFS) $(HTMLS)
	rm -rf lb_walkthroughs/*_files/ findings/2026-03-31-pdac-data-audit/*_files/

.PHONY: all clean
