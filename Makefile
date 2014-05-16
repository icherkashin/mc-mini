EXECUTABLE   := mc-mini

MODULES      = parser matrixForms geometry problem output

# Build directories and files
BUILD_DIR    := build/
TARGET       := $(addprefix $(BUILD_DIR),$(EXECUTABLE))
# Source directories and files
SRC_DIR      := source/
SRC_DIRS     := $(foreach mdir,$(MODULES),$(addprefix $(SRC_DIR),$(mdir)/)) $(SRC_DIR)
SOURCES      := $(foreach sdir,$(SRC_DIRS),$(wildcard $(sdir)*.cpp))
# include directories and files
INCLUDE_DIR  := include/
INCLUDE_DIRS := $(foreach mdir,$(MODULES),$(addprefix $(INCLUDE_DIR),$(mdir)/)) $(INCLUDE_DIR)
INCLUDES     := $(foreach idir,$(INCLUDE_DIRS),$(wildcard $(idir)/*.h))
# Object directories and files
OBJ_DIR      := $(addprefix $(BUILD_DIR),obj/)
OBJ_DIRS     := $(foreach mdir,$(MODULES),$(addprefix $(OBJ_DIR),$(mdir)/)) $(OBJ_DIR)
OBJECTS      := $(foreach src,$(SOURCES),$(patsubst $(SRC_DIR)%.cpp,$(OBJ_DIR)%.o,$(src)))

# C/C++ compiler
CC           := h5c++
# C/C++ compiler flags
CFLAGS       += -Wall -c -Iinclude -fopenmp -std=c++11 -fopenmp 

# C/C++ linker
LD           := h5c++
# C/C++ linker flags
LDFLAGS      += -fopenmp 


vpath %.cpp $(patsubst ' ',':',$(SRC_DIRS))
vpath %.o $(patsubst ' ',':',$(OBJ_DIRS))

define make-goal
$1/%.o %.cpp
  $(CC) $(CFLAGS) -c $$< -o $$@
endef

.PHONY: all checkdirs clean style test

all : checkdirs $(TARGET)

$(OBJ_DIR)%.o : %.cpp
	@echo =====\($(EXECUTABLE)\)====================== Building $<
	@$(CC) $(CFLAGS) -c $< -o $@

$(TARGET) : $(OBJECTS)
	@echo =====\($(EXECUTABLE)\)====================== Linking $@
	@$(LD) $(LDFLAGS) $^ -o $@

checkdirs : $(BUILD_DIR) $(OBJ_DIRS)

$(BUILD_DIR) $(OBJ_DIRS):
	@mkdir -p $@

clean : 
	@rm -rf $(BUILD_DIR)

style :

test :
	@echo $(SRC_DIR) $(SRC_DIRS) $(SOURCES)
# DO NOT DELETE

source/geometry/geometry.o: include/geometry/geometry.h
source/geometry/geometry.o: include/parser/parser.h
source/geometry/geometry.o: include/parser/paramTree.h
