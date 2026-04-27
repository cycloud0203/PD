CXX = g++
CXXFLAGS = -std=c++17 -O3 -fno-strict-aliasing -Wall -pthread
SRCDIR = src
BINDIR = bin
TARGET = $(BINDIR)/fp

SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)

# Plain "touch" can leave mtimes from the file server clock; GNU make uses this host's
# clock and then warns "modification time in the future". Force mtimes to this host.
MTIME := $(shell date +%s)
$(shell touch -d "@$(MTIME)" Makefile $(SRCS) 2>/dev/null; \
	for o in $(OBJS); do test -f "$$o" && touch -d "@$(MTIME)" "$$o"; done; \
	test -d $(BINDIR) && touch -d "@$(MTIME)" $(BINDIR) || true; \
	test -f $(TARGET) && touch -d "@$(MTIME)" $(TARGET) || true)

all: $(TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ -pthread
	touch -d "@$(shell date +%s)" $@ $(BINDIR)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	touch -d "@$(shell date +%s)" $@

$(BINDIR):
	mkdir -p $(BINDIR)
	touch -d "@$(shell date +%s)" $(BINDIR)

clean:
	rm -f $(SRCDIR)/*.o $(TARGET)

.PHONY: all clean
