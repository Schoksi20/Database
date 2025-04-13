CC=g++ -flto -O3
CFLAGS=-c -Iinclude -std=c++11

# Source and include directories
SRC_DIR=src
INCLUDE_DIR=include

# Create include directory if it doesn't exist
$(shell if not exist $(INCLUDE_DIR) mkdir $(INCLUDE_DIR))

all: prepare-dirs MC-BRB

prepare-dirs:
	@if not exist build mkdir build
	@if not exist $(INCLUDE_DIR) mkdir $(INCLUDE_DIR)
	@copy $(SRC_DIR)\*.h $(INCLUDE_DIR)\ >nul 2>&1 || (echo No header files to copy)

MC-BRB: build/main.o build/Graph.o build/HelperAlgorithms.o build/MCBRBSolver.o
	$(CC) build/main.o build/Graph.o build/HelperAlgorithms.o build/MCBRBSolver.o -o MC-BRB
	@if exist build\*.o del build\*.o  # Remove .o files in build directory

build/main.o: $(SRC_DIR)/main.cpp $(INCLUDE_DIR)/Graph.h $(INCLUDE_DIR)/MCBRBSolver.h
	@if not exist build mkdir build
	$(CC) $(CFLAGS) -o build/main.o $(SRC_DIR)/main.cpp

build/Graph.o: $(SRC_DIR)/Graph.cpp $(INCLUDE_DIR)/Graph.h
	@if not exist build mkdir build
	$(CC) $(CFLAGS) -o build/Graph.o $(SRC_DIR)/Graph.cpp

build/HelperAlgorithms.o: $(SRC_DIR)/HelperAlgorithms.cpp $(INCLUDE_DIR)/HelperAlgorithms.h $(INCLUDE_DIR)/Graph.h
	@if not exist build mkdir build
	$(CC) $(CFLAGS) -o build/HelperAlgorithms.o $(SRC_DIR)/HelperAlgorithms.cpp

build/MCBRBSolver.o: $(SRC_DIR)/MCBRBSolver.cpp $(INCLUDE_DIR)/MCBRBSolver.h $(INCLUDE_DIR)/Graph.h $(INCLUDE_DIR)/HelperAlgorithms.h
	@if not exist build mkdir build
	$(CC) $(CFLAGS) -o build/MCBRBSolver.o $(SRC_DIR)/MCBRBSolver.cpp

clean:
	@if exist build\*.o del build\*.o  # Remove .o files
	@if exist build rmdir /s /q build  # Remove the build directory
	@if exist MC-BRB.exe del MC-BRB.exe  # Remove the executable
	@if not exist build mkdir build  # Recreate the build directory after cleaning