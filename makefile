CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -fopenmp
INCLUDES = -Isrc -I/usr/include/gsl
LDFLAGS = -lgsl -lgslcblas -lm

SRC_DIR = src
CORE_DIR = $(SRC_DIR)/core
UTILS_DIR = $(SRC_DIR)/utils
PRICING_DIR = $(SRC_DIR)/pricing

SRCS = \
	$(CORE_DIR)/UpOutCallOption.cpp \
	$(UTILS_DIR)/MathUtils.cpp \
	$(PRICING_DIR)/CustomerTridiagonalSolver.cpp \
	$(PRICING_DIR)/GslTridiagonalSolver.cpp \
	main.cpp

OBJS = $(SRCS:.cpp=.o)
TARGET = option_pricing

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
