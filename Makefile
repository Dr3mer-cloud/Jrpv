FLAGS = -std=c99 -Wall -Wshadow -Wextra -Werror -mfpmath=sse -fstack-protector-all -W -Wunused -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wold-style-definition -Wdeclaration-after-statement -Wbad-function-cast -Wnested-externs -Wuninitialized -Winit-self -O2

FLAGS_D = -std=c99 -Wall -Wshadow -Wextra -Werror -mfpmath=sse -fstack-protector-all -W -Wunused -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wold-style-definition -Wdeclaration-after-statement -Wbad-function-cast -Wnested-externs -Wuninitialized -Winit-self -g


all: pr pr_d

pr: main.o solve.o tests.o
	gcc main.o solve.o tests.o -lm -O2 -o pr

main.o: header.h tests.h solve.h main.c
	gcc $(FLAGS) -c main.c -o $@

solve.o: solve.h solve.c
	gcc $(FLAGS) -c solve.c -o $@

tests.o: solve.h tests.h tests.c
	gcc $(FLAGS) -c tests.c -o $@


pr_d: main_d.o solve_d.o tests_d.o
	gcc main_d.o solve_d.o tests_d.o -lm -g -o pr_d

main_d.o: header.h tests.h solve.h main.c
	gcc $(FLAGS_D) -c main.c -o $@

solve_d.o: solve.h solve.c
	gcc $(FLAGS_D) -c solve.c -o $@

tests_d.o: solve.h tests.h tests.c
	gcc $(FLAGS_D) -c tests.c -o $@

clean:
	rm *.o pr*
