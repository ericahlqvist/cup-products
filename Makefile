# # Generic Makefile for PARI programs -- arm64 running darwin (aarch64 kernel) 64-bit version
# #
# #  This file was created by Configure. Any change made to it will be
# #  lost when Configure is run.
# #
# # make all will create
# #  extgcd-dyn (linked dynamically with libpari)
# #  extgcd-sta (linked statically)
# #  libextgcd.so (to be used by "install" under GP)
# #
# # Under GP: install("extgcd", "GG&&", "gcdex", "./libextgcd.so") enables
# # you to subsequently use gcdex to call extgcd (see the reference manual).
# #

# # change this TARGET to compile your own programs
# TARGET = main-pol
# SHELL  = /bin/sh

# DBGFLAGS   = -g -Wall
# CFLAGS     = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer
# EXTRACFLAGS=
# #CFLAGS    = $(DBGFLAGS)

# # Various linkers use different flags to force static compilation. Choose
# # the one which is relevant for your installation.
# #
# # Solaris ld (global)
# #STATIC    = -dn

# # Solaris ld (toggle: no shared object accepted until -B dynamic is seen
# #STATIC    = -B static

# # gcc
# #STATIC    = -static

# CC         = /usr/bin/gcc
# CPPFLAGS   = -I. -I/usr/local/include -I/opt/homebrew/Cellar/pari/2.15.5/include
# LD         = /usr/bin/gcc
# LDFLAGS    = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer    -Wl,-search_paths_first 
# MODLD      = /usr/bin/gcc
# MODLDFLAGS = -bundle -undefined dynamic_lookup $(CFLAGS) $(DLCFLAGS)
# EXTRAMODLDFLAGS = 
# EXTRALIBS  =

# RUNPTH     = 
# DLCFLAGS   = -fPIC
# LIBS       = -L/opt/homebrew/Cellar/pari/2.15.5/lib -lpari

# RM = rm -f


# OBJS = $(TARGET).o
# # DYN = lib$(TARGET).dylib
# ALL = $(TARGET)-sta # $(TARGET)-dyn $(DYN)

# # dft: $(TARGET)-dyn

# all: $(ALL)

# sta: $(TARGET)-sta

# # dyn: $(TARGET)-dyn

# # dynlib: $(DYN)

# # $(DYN): $(OBJS)
# # 	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)

# $(TARGET)-sta: $(OBJS)
# 	$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(STATIC) $(LIBS)

# # $(TARGET)-dyn: $(OBJS)
# # 	$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(RUNPTH) $(LIBS)

# %.o: %.c
# 	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
# clean:
# 	-$(RM) *.o $(ALL) main-pol-sta
#--------------------------------------------------------------
#--------------------------------------------------------------
#------ Non-parallelized----------------------------------
#--------------------------------------------------------------
# Compilation target.
TARGET = main-pol

# Compiler.
CC = clang

# Compiler flags. /Users/eric/Documents/Matematik/pari/GPDIR/include /Users/eric/Documents/Matematik/pari/GPDIR/lib
CPPFLAGS   = -I. -I/usr/local/include -I/Users/eric/Documents/Matematik/pari/GPDIR/include #-I/opt/homebrew/Cellar/pari/2.17.1/include
CFLAGS = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -pipe -flto=thin -march=native
LDFLAGS = -Wl,-O3 
LIBS = -lm -L/Users/eric/Documents/Matematik/pari/GPDIR/lib -lpari  #-L/opt/homebrew/Cellar/pari/2.17.1/lib -lpari 
#STATIC = -static

# Compilation.
$(TARGET): $(TARGET).c
	$(CC) -o $@ $< $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

clean:
	-$(RM) *.o $(ALL) main-pol

# #--------------------------------------------------------------
# #------ Parallelized----------------------------------
# #--------------------------------------------------------------
# # Compilation target.
# TARGET = main-pol

# # Compiler.
# CC = clang

# # Compiler flags.
# CPPFLAGS   = -I. -I/usr/local/include -I/opt/homebrew/Cellar/pari/2.15.5/include
# CFLAGS = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -pipe -flto=thin -march=native -pthread
# LDFLAGS = -Wl,-O3 -pthread
# LIBS = -lm -L/opt/homebrew/Cellar/pari/2.15.5/lib -lpari 
# #STATIC = -static

# # Compilation.
# $(TARGET): $(TARGET).c
# 	$(CC) -o $@ $< $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

# clean:
# 	-$(RM) *.o $(ALL) $(TARGET)