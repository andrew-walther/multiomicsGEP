SRCS  := update_beta.md update_L.md update_F.md update_tau.md
PDFS  := $(SRCS:.md=.pdf)
HTMLS := $(SRCS:.md=.html)

all: $(PDFS) $(HTMLS)

%.pdf: %.md
	pandoc $< --pdf-engine=xelatex -V geometry:margin=1in \
	  -V mainfont="Helvetica" -V monofont="Menlo" \
	  --toc --syntax-highlighting=tango -o $@

%.html: %.md
	pandoc $< --standalone --toc --syntax-highlighting=tango \
	  --metadata title="$*" -o $@

clean:
	rm -f $(PDFS) $(HTMLS)

.PHONY: all clean
