COMPILER = g++
CFLAGS = -std=c++0x

pointset:	main.cpp point.cpp
	$(COMPILER) $(CFLAGS) -o main main.cpp point.cpp
	@echo "Done compiling pointset code."