CFLAGS = -Wall -Wextra -Werror -g -std=c11
LIBS = -lcheck
CLIBS = -lm
FLAGs = -Wall -Werror -Wextra -std=c11 -g
ifeq ($(shell uname), Linux)
	CLIBS=-lsubunit -lrt -lpthread -lm
endif

all: test


s21_math.o: s21_math.c s21_math.h
	@gcc $(FLAGs) -c s21_math.c 

s21_math.a: s21_math.o
	@gcc -c s21_math.c
	@ar rcs s21_math.a s21_math.o

test: s21_math.a
	gcc s21_tests.c -o check.out  ./s21_math.a $(LIBS) $(CLIBS)
	./check.out

clean:
	- rm -rf *.o *.a *.h.gch check.out s21_math.gcno gcov_test gcov_test.dSYM s21_math.gcda 
	- rm s21_tests.g*
	- rm *.info
	- rm -rf coverage

gcov_report : s21_math.a
	@gcc $(FLAGs) --coverage s21_tests.c s21_math.c  s21_math.a $(LIBS) $(CLIBS) -o gcov_test
	@chmod +x *
	@./gcov_test
	@lcov -t "gcov_test" -o gcov_test.info --no-external -c -d .
	@genhtml -o report/ gcov_test.info
	- @open ./report/index.html


lint:
	cp ./../materials/linters/.clang-format .
	clang-format -n *.c *.h
	clang-format -i *.c *.h
	clang-format -i *.c *.h
	rm .clang-format


