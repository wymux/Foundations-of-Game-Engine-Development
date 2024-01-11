all: src/main

src/main: src/main.o
	mkdir $(dir $@) -p
	g++ -o $@ $^
