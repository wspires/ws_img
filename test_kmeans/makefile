##################### makefile ################################################
#
#    Author:     Wade Spires
#    Program:    test_kmeans
#    Date:       2005/08/04
#    Description:
# 	  Makes test_kmeans.
#    Cleans directory of object files and executables.
#
###############################################################################

###############################################################################
# Variable declaration
###############################################################################

# name of program
NAME = test_kmeans

# compilers/archivers to use
C  = gcc
CC = g++
AR = /usr/bin/ar

# compiler flags
FLAGS =
FLAGS += -g
FLAGS += -O2
FLAGS += -Wall

# my own libraries
TOOLS_DIR = /home/wade/cpp/ws_tools
IMAGE_DIR = ..

# directories to search for header files
INCLUDES = 
INCLUDES += -I$(TOOLS_DIR)
INCLUDES += -I$(IMAGE_DIR)
#INCLUDES += -I`pkg-config --cflags gsl`  # for GNU Scientific Laboratory (GSL)
#INCLUDES += -I<include_directory>

# directories to search for library files
LINK_DIRS =
LINK_DIRS += -L$(TOOLS_DIR)
LINK_DIRS += -L$(IMAGE_DIR)
#LINK_DIRS += -L<link_directory>

# libraries to link with
LIBS =
LIBS += -lm
LIBS += -lws_img  # libimage must precede libtools since image depends on tools
LIBS += -lws_tools
LIBS += -lgsl -lgslcblas  # for GNU Scientific Laboratory (GSL)
LIBS += -ljpeg
LIBS += -lavformat -lavcodec -lavutil -lz  # FFmpeg
LIBS += -lfftw3
#LIBS += -lqui -lqt-mt  # for QT
#LIBS += -lpthread
#LIBS += -l<library>

# loader flags
LD_FLAGS =
#LD_FLAGS +=

# complete set of options to pass to linker
LINK = $(LINK_DIRS) $(LIBS) $(LD_FLAGS)

# name of file containing main()
MAIN = main

# header files in program
HEADERS =

# source code in program
SOURCES = 
SOURCES += $(MAIN).cpp

# object code to generate
OBJECTS =

RM = /bin/rm -f

###############################################################################
# Rules for compiling
###############################################################################

# compile each source file into object code
.c.o:
		$(C)  -c $(FLAGS) $< $(INCLUDES)

.cc.o:
		$(CC) -c $(FLAGS) $< $(INCLUDES)

.SUFFIXES: .cpp .o
.cpp.o:
		$(CC) -c $(FLAGS) $< $(INCLUDES)

# link all object modules into executable
$(NAME): $(MAIN).o $(OBJECTS)
		$(CC) $(FLAGS) -o $(NAME) $(MAIN).o $(OBJECTS) $(LINK) $(INCLUDES)

###############################################################################
# Rules for other stuff
###############################################################################

# create static library (excludes $(MAIN).o from library)
lib: $(OBJECTS)
	$(AR) rcs lib$(NAME).a $(OBJECTS)
#	ranlib lib$(NAME).a  # `ar s` is same as `ranlib`

# create dependency list by examining header files that are included
depend:
	makedepend -- $(FLAGS) -- $(SOURCES) $(INCLUDES) -s'# DO NOT DELETE THIS LINE -- `makedepend` depends on it.'

# remove object files, executables, and libraries
clean:
	$(RM) ${OBJECTS}
	$(RM) ${MAIN}.o
	$(RM) ${NAME}
	$(RM) lib${NAME}.a

# DO NOT DELETE THIS LINE -- `makedepend` depends on it.
