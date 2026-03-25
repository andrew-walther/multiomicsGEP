SRCS  := update_beta.qmd update_L.qmd update_F.qmd update_tau.qmd
PDFS  := $(SRCS:.qmd=.pdf)
HTMLS := $(SRCS:.qmd=.html)

all: $(PDFS) $(HTMLS)

%.pdf: %.qmd
	quarto render $< --to pdf

%.html: %.qmd
	quarto render $< --to html

clean:
	rm -f $(PDFS) $(HTMLS)
	rm -rf *_files/

.PHONY: all clean
