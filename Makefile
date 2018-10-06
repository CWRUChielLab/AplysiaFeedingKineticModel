EXEC=model

$(EXEC): Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ -o $@ "$<"

.PHONY: debug
debug: Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ -D debug -o $(EXEC) "$<"


.PHONY: clean
clean:
	rm -rf $(EXEC) SlugOutput2.txt SlugOutput2-csvTranslated.csv Figures plot.pdf Izhikevich.txt Izhikevich-csvTranslated.csv animationinfo.txt animationinfo-csvTranslated.csv animation.mp4 rasterplotinfo.txt rasterplotinfo-csvTranslated.csv RasterPlot.pdf output-mechanics.csv *~ .*~


.PHONY: check
check: $(EXEC)
	./$< SwallowPerturbed
	@echo "SwallowPerturbed:" `diff SlugOutput2.txt Check/Output-SwallowPerturbed.txt | wc -l` "deviations"
	@echo
	./$< Bite
	@echo "Bite:            " `diff SlugOutput2.txt Check/Output-Bite.txt             | wc -l` "deviations"
	@echo
	./$< RejectionA
	@echo "RejectionA:      " `diff SlugOutput2.txt Check/Output-RejectionA.txt       | wc -l` "deviations"
	@echo
	./$< RejectionB
	@echo "RejectionB:      " `diff SlugOutput2.txt Check/Output-RejectionB.txt       | wc -l` "deviations"
	@echo
	./$< SwallowA
	@echo "SwallowA:        " `diff SlugOutput2.txt Check/Output-SwallowA.txt         | wc -l` "deviations"
	@echo
	./$< SwallowB
	@echo "SwallowB:        " `diff SlugOutput2.txt Check/Output-SwallowB.txt         | wc -l` "deviations"


.PHONY: blessBite
blessBite: $(EXEC)
	./$< Bite
	mv SlugOutput2.txt Check/Output-Bite.txt

.PHONY: blessRejectionA
blessRejectionA: $(EXEC)
	./$< RejectionA
	mv SlugOutput2.txt Check/Output-RejectionA.txt

.PHONY: blessRejectionB
blessRejectionB: $(EXEC)
	./$< RejectionB
	mv SlugOutput2.txt Check/Output-RejectionB.txt

.PHONY: blessSwallowA
blessSwallowA: $(EXEC)
	./$< SwallowA
	mv SlugOutput2.txt Check/Output-SwallowA.txt

.PHONY: blessSwallowB
blessSwallowB: $(EXEC)
	./$< SwallowB
	mv SlugOutput2.txt Check/Output-SwallowB.txt

.PHONY: blessSwallowPerturbed
blessSwallowPerturbed: $(EXEC)
	./$< SwallowPerturbed
	mv SlugOutput2.txt Check/Output-SwallowPerturbed.txt


FIGDIR=Figures

.PHONY: figures
figures:
	make $(FIGDIR)/Plot-Bite.pdf
	make $(FIGDIR)/Plot-RejectionA.pdf
	make $(FIGDIR)/Plot-RejectionB.pdf
	make $(FIGDIR)/Plot-SwallowA.pdf
	make $(FIGDIR)/Plot-SwallowB.pdf
	make $(FIGDIR)/Plot-SwallowPerturbed.pdf

$(FIGDIR)/Plot-Bite.pdf: $(EXEC)
	./$< Bite
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-RejectionA.pdf: $(EXEC)
	./$< RejectionA
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-RejectionB.pdf: $(EXEC)
	./$< RejectionB
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-SwallowA.pdf: $(EXEC)
	./$< SwallowA
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-SwallowB.pdf: $(EXEC)
	./$< SwallowB
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-SwallowPerturbed.pdf: $(EXEC)
	./$< SwallowPerturbed
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@
