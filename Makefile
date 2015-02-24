COMPILER = g++
CFLAGS = -std=c++0x

pointset:	main.cpp Point.cpp
	$(COMPILER) $(CFLAGS) -o pointset main.cpp Point.cpp
	@echo "Done compiling pointset code."