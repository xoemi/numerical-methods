FLAGS = -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -std=gnu99 -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wold-style-definition -Wdeclaration-after-statement -Wbad-function-cast -Wnested-externs 

all: a1.out a2.out a3.out

f.o: f.c f.h
	gcc -c $(FLAGS) f.c 

task1.o: task1.c f.c task1.h f.h
	gcc -c $(FLAGS) task1.c f.c

task2.o: task2.c f.c task2.h f.h
	gcc -c $(FLAGS) task2.c

task3.o: task3.c f.c task3.h f.h
	gcc -c $(FLAGS) task3.c

main1.o: main1.c task1.h
	gcc -c $(FLAGS) main1.c

main2.o: main2.c task2.h
	gcc -c $(FLAGS) main2.c

main3.o: main3.c task3.h
	gcc -c $(FLAGS) main3.c

a1.out: main1.o task1.o f.o
	gcc main1.o task1.o f.o -o a1.out

a2.out: main2.o task2.o f.o
	gcc main2.o task2.o f.o -o a2.out

a3.out: main3.o task3.o f.o
	gcc main3.o task3.o f.o -o a3.out
