# Compiler and tools
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11 -fopenmp -g
LDFLAGS = -lgsl -lgslcblas -ldl

# Source files and objects
SRCS = TestUoc.cpp UpOutCallOption.cpp optimizer.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = TestUoc.out

# Default target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET)

.PHONY: all clean
