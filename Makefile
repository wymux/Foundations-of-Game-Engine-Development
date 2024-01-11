all: run

src/main: src/main.o
	mkdir $(dir $@) -p
	g++ -o $@ $^

run: src/main
	./src/main
