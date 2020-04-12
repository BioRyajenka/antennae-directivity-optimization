PROJECT=sphereml
VERSION=1.0
CXXFLAGS=-MD -DNDEBUG -O3 -Wall -std=c++11 -fPIC
# LDFLAGS=-lpybind11
SRC_DIR := ./
OUT_DIR := build
OBJ_DIR := $(OUT_DIR)
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
SRC_MPI := $(SRC_DIR)/joptimize.cpp  $(SRC_DIR)/jade.cpp 
SRC_PY := $(SRC_DIR)/pybind_sphereml.cpp
SRC_CC := $(filter-out $(SRC_MPI) $(SRC_PY), $(SRC_FILES))

OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
OBJ_MPI := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_MPI))
OBJ_PY := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_PY))
OBJ_CC := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CC))
OBJ_MAINS := $(OBJ_DIR)/main.o $(OBJ_DIR)/pybind_sphereml.o $(OBJ_DIR)/joptimize.o

DEPS=$(OBJ_FILES:$(OBJ_DIR)/%.o=$(OBJ_DIR)/%.d)

-include $(DEPS)


all: directivity lib joptimize

.PHONY : clean

directivity: $(OBJ_DIR)/main.o $(filter-out $(OBJ_MAINS) $(OBJ_MPI), $(OBJ_FILES))
	c++ $(LDFLAGS) -o $@ $^ -std=c++11 

joptimize: $(OBJ_DIR)/joptimize.o $(filter-out $(OBJ_MAINS) $(OBJ_PY), $(OBJ_FILES))
	mpic++ $(LDFLAGS) -o $@ $^ -std=c++11

lib: $(OBJ_DIR)/pybind_sphereml.o $(filter-out $(OBJ_MAINS)  $(OBJ_MPI), $(OBJ_FILES))
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` $^ -o sphereml`python3-config --extension-suffix`

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	@echo -n Comiling $< ...
	@if [ "" != "$(findstring $<,$(SRC_MPI))" ]; then   \
		mpic++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<; \
		echo MPI compiled; \
	elif [ "" != "$(findstring $<,$(SRC_PY))" ]; then  \
		c++ $(CPPFLAGS) -shared $(CXXFLAGS) `python3 -m pybind11 --includes` -c -o $@ $<; \
		echo Python compiled; \
	else                                                \
		c++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<;    \
		echo CC compiled; \
	fi

clean:
	rm -rf directivity
	rm -rf $(OUT_DIR)
	find . -name '*.pyc' -delete
	find . -name '*.o' -delete
	find . -name '*.d' -delete
	find . -name '*.so' -delete
	find . -name '*~' -delete
