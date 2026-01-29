# Makefile for "simprot"

NAME = simprot
CC = g++
CFLAGS = -O3
LDLIBS = -lpopt -lm

all:	$(NAME).cpp eigen.h random.c random.h
	$(CC) -o $(NAME) $(NAME).cpp random.c -I/opt/local/include -L/usr/local/lib $(LDLIBS)

clean:
	rm $(NAME) *.o *~
.PHONY:	clean
