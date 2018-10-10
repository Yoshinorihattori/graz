TARGET = abst.pdf
BASENAME = $(basename $(TARGET))
all: $(TARGET)

.SUFFIXES: .pdf .tex

$(TARGET): $(BASENAME).tex
	platex $(BASENAME)
	dvipdfmx -p a4 -o $(TARGET) $(BASENAME)
	rm -rf $(BASENAME).aux $(BASENAME).log $(BASENAME).dvi *~

clean:
	rm -rf $(BASENAME).aux $(BASENAME).log $(BASENAME).dvi *~
	rm -rf $(TARGET)
