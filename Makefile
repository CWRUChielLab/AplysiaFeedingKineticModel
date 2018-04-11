EXEC=model

$(EXEC): Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ -o $@ "$<"


.PHONY: clean
clean:
	rm -rf $(EXEC) SlugOutput2.txt SlugOutput2-csvTranslated.csv Figures plot.pdf Izhikevich.txt Izhikevich-csvTranslated.csv animationinfo.txt animationinfo-csvTranslated.csv animation.mp4 rasterplotinfo.txt rasterplotinfo-csvTranslated.csv RasterPlot.pdf *~ .*~


.PHONY: check
check: $(EXEC)
	./$(EXEC) SwallowPerturbed
	@echo "SwallowPerturbed:" `diff SlugOutput2.txt Check/Output-SwallowPerturbed.txt | wc -l` "deviations"
	@echo
	./$(EXEC) Bite
	@echo "Bite:            " `diff SlugOutput2.txt Check/Output-Bite.txt             | wc -l` "deviations"
	@echo
	./$(EXEC) RejectionA
	@echo "RejectionA:      " `diff SlugOutput2.txt Check/Output-RejectionA.txt       | wc -l` "deviations"
	@echo
	./$(EXEC) RejectionB
	@echo "RejectionB:      " `diff SlugOutput2.txt Check/Output-RejectionB.txt       | wc -l` "deviations"
	@echo
	./$(EXEC) SwallowA
	@echo "SwallowA:        " `diff SlugOutput2.txt Check/Output-SwallowA.txt         | wc -l` "deviations"
	@echo
	./$(EXEC) SwallowB
	@echo "SwallowB:        " `diff SlugOutput2.txt Check/Output-SwallowB.txt         | wc -l` "deviations"


.PHONY: blessBite
blessBite: $(EXEC)
	./$(EXEC) Bite
	mv SlugOutput2.txt Check/Output-Bite.txt

.PHONY: blessRejectionA
blessRejectionA: $(EXEC)
	./$(EXEC) RejectionA
	mv SlugOutput2.txt Check/Output-RejectionA.txt

.PHONY: blessRejectionB
blessRejectionB: $(EXEC)
	./$(EXEC) RejectionB
	mv SlugOutput2.txt Check/Output-RejectionB.txt

.PHONY: blessSwallowA
blessSwallowA: $(EXEC)
	./$(EXEC) SwallowA
	mv SlugOutput2.txt Check/Output-SwallowA.txt

.PHONY: blessSwallowB
blessSwallowB: $(EXEC)
	./$(EXEC) SwallowB
	mv SlugOutput2.txt Check/Output-SwallowB.txt

.PHONY: blessSwallowPerturbed
blessSwallowPerturbed: $(EXEC)
	./$(EXEC) SwallowPerturbed
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
	./$(EXEC) Bite
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-RejectionA.pdf: $(EXEC)
	./$(EXEC) RejectionA
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-RejectionB.pdf: $(EXEC)
	./$(EXEC) RejectionB
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-SwallowA.pdf: $(EXEC)
	./$(EXEC) SwallowA
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-SwallowB.pdf: $(EXEC)
	./$(EXEC) SwallowB
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@

$(FIGDIR)/Plot-SwallowPerturbed.pdf: $(EXEC)
	./$(EXEC) SwallowPerturbed
	python PlotVariablesSlugOutput.py
	mkdir -p $(FIGDIR)
	mv plot.pdf $@