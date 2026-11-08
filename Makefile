CXX ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude

BUILD := build
SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:src/%.cpp=$(BUILD)/obj/%.o)
LIB := $(BUILD)/libphys.a

.PHONY: all test demo clean

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

clean:
	rm -rf $(BUILD)
