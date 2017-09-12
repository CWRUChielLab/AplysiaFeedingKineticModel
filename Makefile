EXEC=model

$(EXEC): Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ -o $@ "$<"

.PHONY: clean
clean:
	rm -rf $(EXEC) SlugOutput2.txt SlugOutput2-csvTranslated.csv plot.pdf *~ .*~

.PHONY: check
check: $(EXEC)
	./$(EXEC) SwallowPerturbed
	@echo “SwallowPerturbed:   ” `diff SlugOutput2.txt Check/Output-SwallowPerturbed.txt | wc -l` “deviations”
	./$(EXEC) Bite
	@echo “Bite:   ” `diff SlugOutput2.txt Check/Output-Bite.txt | wc -l` “deviations”
	./$(EXEC) RejectionA
	@echo “RejectionA:   ” `diff SlugOutput2.txt Check/Output-RejectionA.txt | wc -l` “deviations”
	./$(EXEC) RejectionB
	@echo “RejectionB:   ” `diff SlugOutput2.txt Check/Output-RejectionB.txt | wc -l` “deviations”
	./$(EXEC) SwallowA
	@echo “SwallowA:   ” `diff SlugOutput2.txt Check/Output-SwallowA.txt | wc -l` “deviations”
	./$(EXEC) SwallowB
	@echo “SwallowB:   ” `diff SlugOutput2.txt Check/Output-SwallowB.txt | wc -l` “deviations”

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