EXEC = model
FIGDIR = Figures
ANIDIR = Animations
OUTPUTS = \
	animation.mp4 \
	animationinfo.csv \
	Izhikevich.csv \
	output-mechanics.csv \
	plot.pdf \
	RasterPlot.pdf \
	rasterplotinfo.csv \
	SlugOutput2.csv


# When a target is not specified, the default executable is built
.PHONY: default
default: $(EXEC)

$(EXEC): DEFINES =
debug:   DEFINES = -D debug

.PHONY: debug
$(EXEC) debug: Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ $(DEFINES) -o $(EXEC) "$<"


.PHONY: clean
clean:
	rm -rf $(EXEC) $(FIGDIR) $(ANIDIR) $(OUTPUTS) *~ .*~


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
ANIBEHAVIORS   = $(addprefix $(ANIDIR)/Animation-, $(addsuffix .mp4, $(BEHAVIORS)))


.PHONY: check
check: $(CHECKBEHAVIORS)

.PHONY: $(CHECKBEHAVIORS)
$(CHECKBEHAVIORS): CHECK-%: $(EXEC)
	@echo -n "Checking $*...\t"
	@./$< $* # run the model for this behavior
	@echo `diff --strip-trailing-cr SlugOutput2.csv Check/Output-$*.csv | wc -l` "deviations"


.PHONY: bless
bless: $(BLESSBEHAVIORS)

.PHONY: $(BLESSBEHAVIORS)
$(BLESSBEHAVIORS): BLESS-%: $(EXEC)
	@echo -n "Blessing $*...\t"
	@./$< $* # run the model for this behavior
	@mv SlugOutput2.csv Check/Output-$*.csv
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


.PHONY: animations
animations: $(ANIBEHAVIORS)

$(ANIBEHAVIORS): $(ANIDIR)/Animation-%.mp4: $(EXEC)
	@echo "Animating $*...\t"
	@./$< $* # run the model for this behavior
	@python newanimation.py
	@mkdir -p $(ANIDIR)
	@mv animation.mp4 $@
