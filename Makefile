EXEC = model

# When a target is not specified, the default executable is built
.PHONY: default
default: $(EXEC)

model: DEFINES =
debug: DEFINES = -D debug

$(EXEC) debug: Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ $(DEFINES) -o $(EXEC) "$<"


.PHONY: clean
clean:
	rm -rf $(EXEC) SlugOutput2.txt SlugOutput2-csvTranslated.csv Figures plot.pdf Izhikevich.txt Izhikevich-csvTranslated.csv animationinfo.txt animationinfo-csvTranslated.csv animation.mp4 rasterplotinfo.txt rasterplotinfo-csvTranslated.csv RasterPlot.pdf output-mechanics.csv *~ .*~


BEHAVIORS = \
	Bite \
	RejectionA \
	RejectionB \
	SwallowA \
	SwallowB \
	SwallowPerturbed
CHECKBEHAVIORS = $(addprefix CHECK-, $(BEHAVIORS))

.PHONY: check
check: $(CHECKBEHAVIORS)

.PHONY: $(CHECKBEHAVIORS)
$(CHECKBEHAVIORS): CHECK-%: $(EXEC)
	@echo -n "Checking $*...\t"
	@./$< $* # run the model for this behavior
	@echo `diff SlugOutput2.txt Check/Output-$*.txt | wc -l` "deviations"


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
