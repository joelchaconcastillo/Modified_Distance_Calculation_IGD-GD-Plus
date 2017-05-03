# Makefile for compiling a the idraw
CC=g++
LD=g++
RM=rm -f
#Production flags
CFLAGS=-Wno-unused-result -std=c++11
LDFLAGS=
UBUNTUFLAGS=-lm $(LDFLAGS) -lstdc++ -std=c++11
OPENSUSEFLAGS=-lm $(LDFLAGS)
DEBUG=-O2
OBJS:=main.o 
MAIN= MPlus
all:$(MAIN)
ubuntu:
	$(MAKE) "OPENSUSEFLAGS="
opensuse:
	$(MAKE) "UBUNTUFLAGS="
$(MAIN):$(OBJS)
	$(LD) $(DEBUG) $(OPENSUSEFLAGS) $(OBJS) -o $(MAIN) $(UBUNTUFLAGS)
%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(DEBUG) -c $<
%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEBUG) -c $<	
%.o: %.c
	$(CC) $(CFLAGS) $(DEBUG) -c $<	
%.o: %.cpp 
	$(CC) $(CFLAGS) $(DEBUG) -c $<	

#$(CC) wfg/WFG1/WFG1.h $(CFLAGS) $(DEBUG) -c -o wfg/WFG1/WFG1.o
debug:
	$(MAKE) "DEBUG=-g"
clean:
	$(RM) $(OBJS) *.*~
cleanall:
	$(RM) $(OBJS) $(MAIN) *.*~
	
