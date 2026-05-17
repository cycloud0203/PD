CXX       := g++
CXXFLAGS  := -std=c++17 -O3 -flto -funroll-loops -fomit-frame-pointer \
              -fopenmp -Wall -D_GLIBCXX_USE_CXX11_ABI=0 \
              -fno-pie

LDFLAGS   := -fopenmp -flto -no-pie
LIBS      := -Llib -lDetailPlace -lGlobalPlace -lLegalizer -lPlacement -lParser -lPlaceCommon

SOURCES   := src/ObjectiveFunction.cpp src/Optimizer.cpp src/GlobalPlacer.cpp src/main.cpp
OBJECTS   := $(patsubst src/%.cpp,src/%.o,$(SOURCES))
EXECUTABLE:= bin/place
SKEW_CHECK_FILES := $(SOURCES) $(OBJECTS) $(EXECUTABLE) Makefile check_density_target.pl

.PHONY: all clean fix-clock-skew

all: fix-clock-skew $(EXECUTABLE)

fix-clock-skew:
	@now=$$(date +%s); \
	for f in $(SKEW_CHECK_FILES); do \
		[ -e "$$f" ] || continue; \
		mtime=$$(stat -c %Y "$$f"); \
		if [ "$$mtime" -gt "$$now" ]; then \
			echo "Normalizing future timestamp: $$f"; \
			touch -d "@$$now" "$$f"; \
		fi; \
	done

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -f src/*.o $(EXECUTABLE)
