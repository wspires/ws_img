##################### makefile ################################################
#
#    Author:     Wade Spires
#    Program:    sift
#    Date:       2005/08/04
#    Description:
# 	  Makes sift.
#    Cleans directory of object files and executables.
#
###############################################################################

###############################################################################
# Declare variables
###############################################################################

NAME = sift

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
#INCLUDES += -I`pkg-config --cflags gsl`  # for GNU Scientific Laboratory (GSL)
#INCLUDES += -I/home/wade/lib/boost-1.30.2
#INCLUDES += -I/home/wade/cpp/matrix
#INCLUDES += -I<include_directory>

LD_FLAGS =
#LD_FLAGS +=

LINK_DIRS =
LINK_DIRS += -L$(TOOLS_DIR)
LINK_DIRS += -L$(IMAGE_DIR)
#LINK_DIRS += -L<link_directory>

LIBS =
LIBS += -lm
LIBS += -lws_img  # libimage must precede libtools since image depends on tools
LIBS += -lws_tools
LIBS += -lgsl -lgslcblas  # for GNU Scientific Laboratory (GSL)
LIBS += -ljpeg
LIBS += -lavformat -lavcodec -lavutil -lz  # FFmpeg
LIBS += -lfftw3
#LIBS += -lqui -lqt-mt  # for QT
#LIBS += -lpthread -lsocket -lnsl
#LIBS += -l<library>

LINK = $(LINK_DIRS) $(LIBS) $(LD_FLAGS)

# name of file containing main()
MAIN = main

HEADERS =

SOURCES = 

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

# create static library
lib: $(OBJECTS)
	$(AR) rcs lib$(NAME).a $(OBJECTS)
#	ranlib lib$(NAME).a  # `ar s` is same as `ranlib`

# create source files' dependency list of header files
depend:
	makedepend -- $(FLAGS) -- $(SOURCES) $(INCLUDES) -s'# DO NOT DELETE THIS LINE -- `makedepend` depends on it.'

# remove object files
clean:
	$(RM) ${OBJECTS}
	$(RM) ${NAME}

# DO NOT DELETE THIS LINE -- `makedepend` depends on it.

