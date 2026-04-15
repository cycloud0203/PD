CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall
SRCS = src/main.cpp src/partitioner.cpp
TARGET = bin/fm

all: $(TARGET)

$(TARGET): $(SRCS) src/partitioner.h
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)
