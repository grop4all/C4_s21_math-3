CFLAGS = -Wall -Wextra -Werror -g
LIBS = -lcheck
CLIBS = -lsubunit -lrt -lpthread -lm

all: clean test

s21_math.o: s21_math.c s21_math.h
	gcc -c $(CFLAGS) s21_math.c

s21_math.a: s21_math.o
	ar r s21_math.a s21*.o
	ranlib s21_math.a
	rm -rf *.o

test: s21_math.a
	gcc s21_tests.c -o check.out  ./s21_math.a $(LIBS) $(CLIBS)
	./check.out

clean:
	rm -rf *.o *.a *.h.gch check.out
	rm s21_tests.g*
	rm *.out
	rm *.info
	rm -rf coverage

gcov_report : s21_math.a
	gcc s21_tests.c -o gcov_rep_check.out ./s21_math.a $(LIBS) $(CLIBS) -fprofile-arcs -ftest-coverage
	./gcov_rep_check.out
	lcov -c -d . -o coverage.info
	genhtml coverage.info -o coverage
	open coverage/index.html


lint:
	cp ./../materials/linters/.clang-format .
	clang-format -n *.c *.h
	clang-format -i *.c *.h
	clang-format -i *.c *.h
	rm .clang-format

