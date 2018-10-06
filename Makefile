EXEC = model
FIGDIR = Figures


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
BLESSBEHAVIORS = $(addprefix BLESS-, $(BEHAVIORS))
FIGBEHAVIORS   = $(addprefix $(FIGDIR)/Plot-, $(addsuffix .pdf, $(BEHAVIORS)))


.PHONY: check
check: $(CHECKBEHAVIORS)

.PHONY: $(CHECKBEHAVIORS)
$(CHECKBEHAVIORS): CHECK-%: $(EXEC)
	@echo -n "Checking $*...\t"
	@./$< $* # run the model for this behavior
	@echo `diff SlugOutput2.txt Check/Output-$*.txt | wc -l` "deviations"


.PHONY: bless
bless: $(BLESSBEHAVIORS)

.PHONY: $(BLESSBEHAVIORS)
$(BLESSBEHAVIORS): BLESS-%: $(EXEC)
	@echo -n "Blessing $*...\t"
	@./$< $* # run the model for this behavior
	@mv SlugOutput2.txt Check/Output-$*.txt
	@echo "done"


.PHONY: figures
figures: $(FIGBEHAVIORS)

$(FIGBEHAVIORS): $(FIGDIR)/Plot-%.pdf: $(EXEC)
	@echo -n "Plotting $*...\t"
	@./$< $* # run the model for this behavior
	@python PlotVariablesSlugOutput.py
	@mkdir -p $(FIGDIR)
	@mv plot.pdf $@
	@echo "done"
