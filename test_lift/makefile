##################### makefile ################################################
#
#    Author:     Wade Spires
#    Program:    test_lift
#    Date:       2005/04/06
#    Description:
# 	  Makes test_lift.
#    Cleans directory of object files and executables.
#
###############################################################################

###############################################################################
# Declare variables
###############################################################################

NAME = test_lift

C  = gcc
CC = g++
AR = /usr/bin/ar

FLAGS =
FLAGS += -g
FLAGS += -O2
FLAGS += -Wall

TOOLS_DIR = /home/wade/cpp/ws_tools
IMAGE_DIR = ..

INCLUDES = 
INCLUDES += -I$(TOOLS_DIR)
INCLUDES += -I$(IMAGE_DIR)
INCLUDES += `pkg-config --cflags gsl`
#INCLUDES += -I<include_directory>

LD_FLAGS =
#LD_FLAGS +=

LINK_DIRS =
LINK_DIRS += -L$(TOOLS_DIR)
LINK_DIRS += -L$(IMAGE_DIR)
#LINK_DIRS += -L<link_directory>

LIBS =
LIBS += -lm
LIBS += -lws_tools
LIBS += -lws_img
LIBS += -lgsl -lgslcblas
LIBS += -ljpeg
LIBS += -lfftw3
#LIBS += -l<library>

LINK = $(LINK_DIRS) $(LIBS) $(LD_FLAGS)

# name of file containing main()
MAIN = main

HEADERS =

SOURCES = 
#SOURCES += 

OBJECTS =
#OBJECTS += 

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

# create source files' dependency list of header files
depend:
	makedepend -- $(FLAGS) -- $(SOURCES) $(INCLUDES) -s'# DO NOT DELETE THIS LINE -- `makedepend` depends on it.'

# remove object files
clean:
	$(RM) ${OBJECTS}
	$(RM) ${NAME}

# DO NOT DELETE THIS LINE -- `makedepend` depends on it.

