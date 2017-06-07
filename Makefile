EXEC=model

$(EXEC): Aplysia\ Feeding\ Kinetic\ Model.cpp
	g++ -o $@ "$<"

.PHONY: clean
clean:
	rm -rf $(EXEC) *~ .*~
