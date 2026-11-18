CXX ?= g++
MPICXX ?= mpicxx
CXXFLAGS ?= -O3 -std=c++17 -fopenmp -Wall -Wextra
SRC_DIR := src
BUILD_DIR := build

SOURCES := common.cpp graph.cpp io.cpp sam_dfs.cpp hamsearch.cpp tsp.cpp verify.cpp visualize.cpp mpi_driver.cpp main.cpp
OMP_OBJECTS := $(addprefix $(BUILD_DIR)/omp/,$(SOURCES:.cpp=.o))
MPI_OBJECTS := $(addprefix $(BUILD_DIR)/mpi/,$(SOURCES:.cpp=.o))
HEADERS := $(wildcard $(SRC_DIR)/*.h)

.PHONY: all omp mpi tools clean test bench-quick

all: omp tools

tools: run_bench gen_synthetic

run_bench: bench/run_bench.cpp
	$(CXX) -O2 -std=c++17 -Wall -Wextra -o $@ $<

gen_synthetic: bench/gen_synthetic.cpp
	$(CXX) -O2 -std=c++17 -Wall -Wextra -o $@ $<

omp: hampath

mpi: hampath_mpi

hampath: $(OMP_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

hampath_mpi: $(MPI_OBJECTS)
	$(MPICXX) $(CXXFLAGS) -DUSE_MPI -o $@ $^

$(BUILD_DIR)/omp/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/mpi/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(MPICXX) $(CXXFLAGS) -DUSE_MPI -c $< -o $@

test: hampath run_bench
	./run_bench --suite smoke --threads 1,2 --out results/smoke.csv

bench-quick: hampath run_bench
	./run_bench --suite smoke,graphdb,tsplib --threads 1,2,4 --time-limit 60 --out results/quick.csv

clean:
	rm -rf $(BUILD_DIR) hampath hampath_mpi run_bench gen_synthetic
