CXX ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude

BUILD := build
SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:src/%.cpp=$(BUILD)/obj/%.o)
LIB := $(BUILD)/libphys.a

.PHONY: all test demo recordings view clean

all: $(LIB) $(BUILD)/tests $(BUILD)/stack_demo

$(BUILD)/obj/%.o: src/%.cpp $(wildcard include/phys/*.h)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJS)
	ar rcs $@ $^

$(BUILD)/tests: tests/tests.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/stack_demo: demos/stack_demo.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test: $(BUILD)/tests
	./$(BUILD)/tests

demo: $(BUILD)/stack_demo
	./$(BUILD)/stack_demo

SCENES := tower dominoes pyramid
RECORDINGS := $(SCENES:%=$(BUILD)/%.json)
PORT ?= 8765

$(BUILD)/%.json: $(BUILD)/stack_demo
	./$(BUILD)/stack_demo $* --record $@

recordings: $(RECORDINGS)

# Records the demo scenes and serves the viewer at http://localhost:$(PORT)/viewer/
view: recordings
	@if python3 -c "import socket,sys; s=socket.socket(); sys.exit(s.connect_ex(('127.0.0.1',$(PORT)))!=0)"; then \
		echo "Port $(PORT) is already in use. If the viewer is already running, open"; \
		echo "http://localhost:$(PORT)/viewer/ or pick another port: make view PORT=8766"; \
		exit 1; \
	fi
	@echo "Viewer: http://localhost:$(PORT)/viewer/"
	python3 -m http.server $(PORT) --bind 127.0.0.1

clean:
	rm -rf $(BUILD)
