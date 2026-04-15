CXX = g++
CXXFLAGS = -std=c++17 -O3 -fno-strict-aliasing -Wall -pthread
SRCDIR = src
BINDIR = bin
TARGET = $(BINDIR)/fp

SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ -pthread

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -f $(SRCDIR)/*.o $(TARGET)

.PHONY: all clean
