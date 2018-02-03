##################### makefile ################################################
#
#    Author:     Wade Spires
#    Program:    ws_img
#    Date:       2005/08/04
#    Description:
# 	  Makes example.
#    Cleans directory of object files and executables.
#
###############################################################################

###############################################################################
# Declare variables
###############################################################################

NAME = ws_img

C  = gcc
CC = g++
AR = ar

FLAGS =
#FLAGS += -g     # debugging info
FLAGS += -O2    # speed/size optimizations
#FLAGS += -O3    # speed optimizations
#FLAGS += -Os    # size optimizations--may help speed by reducing cache misses
FLAGS += -Wall  # enable all warnings

TOOLS_DIR   = /home/wade/cpp/ws_tools
CURRENT_DIR = $(shell /bin/pwd)
SIFT_DIR    = $(CURRENT_DIR)/sift

# whether GSL, jpeg, or FFmpeg are installed (1) or are not installed (0)
HAVE_GSL    = 1
HAVE_JPEG   = 1
HAVE_FFMPEG = 1
HAVE_FFTW   = 1

INCLUDES = 
INCLUDES += -I$(TOOLS_DIR)

# GNU Scientific Library (GSL): disable range-checking for GSL matrices
ifeq ($(HAVE_GSL),1)
	INCLUDES += -DHAVE_GSL=1
	#INCLUDES += `pkg-config --cflags gsl`
	INCLUDES += -DHAVE_INLINE=1 -DGSL_RANGE_CHECK=0
endif
ifeq ($(HAVE_FFTW),1)
	INCLUDES += -DHAVE_FFTW=1
endif

INCLUDES += -DHAVE_LIBJPEG=$(HAVE_JPEG)

LD_FLAGS =
#LD_FLAGS +=

LINK_DIRS =
LINK_DIRS += -L$(TOOLS_DIR)
#LINK_DIRS += -L<link_directory>

LIBS =
LIBS += -lm
LIBS += -lws_tools
ifeq ($(HAVE_GSL),1)
	LIBS += -lgsl -lgslcblas
endif
ifeq ($(HAVE_JPEG),1)
	LIBS += -ljpeg
endif
ifeq ($(HAVE_FFMPEG),1)
	LIBS += -lavformat -lavcodec -lavutil -lz
endif
ifeq ($(HAVE_FFTW),1)
	LIBS += -lfftw3
endif
#LIBS += -l<library>

LINK = $(LINK_DIRS) $(LIBS) $(LD_FLAGS)

HEADERS =
HEADERS += Image.hpp
HEADERS += Image_Coordinate.hpp
HEADERS += Image_Window.hpp
HEADERS += Matrix.hpp
HEADERS += Matrix_Region.hpp
HEADERS += Vector.hpp
HEADERS += Linear_Algebra.hpp
HEADERS += Lift_Scheme.hpp
HEADERS += Wavelet.hpp
HEADERS += Wavelet_Subband.hpp
HEADERS += Edge_Detector.hpp
HEADERS += Jpeg_IO.hpp
HEADERS += Harris.hpp
HEADERS += Gauss.hpp
HEADERS += Gabor.hpp
HEADERS += Wavelet_Coefficients.hpp
ifeq ($(HAVE_FFMPEG),1)
	HEADERS += Video_IO.hpp
endif
HEADERS += Integral_Image.hpp
HEADERS += K_Means.hpp
HEADERS += Pnm_IO.hpp

SOURCES = 
SOURCES += Image.cpp
SOURCES += Matrix.cpp
SOURCES += Matrix_Region.cpp
SOURCES += Vector.cpp
SOURCES += Linear_Algebra.cpp
SOURCES += Lift_Scheme.cpp
SOURCES += Wavelet.cpp
SOURCES += Wavelet_Subband.cpp
SOURCES += Edge_Detector.cpp
SOURCES += Jpeg_IO.cpp
SOURCES += Harris.cpp
SOURCES += Gauss.cpp
SOURCES += Gabor.cpp
SOURCES += Wavelet_Coefficients.cpp
ifeq ($(HAVE_FFMPEG),1)
	SOURCES += Video_IO.cpp
endif
SOURCES += Integral_Image.cpp
SOURCES += K_Means.cpp
SOURCES += Pnm_IO.cpp

OBJECTS =
OBJECTS += Image.o
OBJECTS += Matrix.o
OBJECTS += Matrix_Region.o
OBJECTS += Vector.o
OBJECTS += Linear_Algebra.o
OBJECTS += Lift_Scheme.o
OBJECTS += Wavelet.o
OBJECTS += Wavelet_Subband.o
OBJECTS += Edge_Detector.o
OBJECTS += Jpeg_IO.o
OBJECTS += Harris.o
OBJECTS += Gauss.o
OBJECTS += Gabor.o
OBJECTS += Wavelet_Coefficients.o
ifeq ($(HAVE_FFMPEG),1)
	OBJECTS += Video_IO.o
endif
OBJECTS += Integral_Image.o
OBJECTS += K_Means.o
OBJECTS += Pnm_IO.o

RM = /bin/rm -f

SUBDIR = 'test_matrix'
SUBDIR += 'test_kmeans'
SUBDIR += 'test_lift'
SUBDIR += 'test_wavelet'
SUBDIR += 'test_sift'
SUBDIR += 'example'

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

# create static library (make libsift, too)
$(NAME): $(OBJECTS)
	make -C $(SIFT_DIR)
	$(AR) rcs lib$(NAME).a $(OBJECTS) `ls $(SIFT_DIR)/*.o`
	@for test in $(SUBDIR); do \
		make -C $$test; \
	done;

# link all object modules into executable
all: $(OBJECTS)
	$(AR) rcs lib$(NAME).a $(OBJECTS)
	@for test in $(SUBDIR); do \
		make -C $$test; \
	done;

###############################################################################
# Rules for other stuff
###############################################################################

# create static library (make libsift, too)
lib: $(OBJECTS)
	make -C $(SIFT_DIR)
	$(AR) rcs lib$(NAME).a $(OBJECTS) `ls $(SIFT_DIR)/*.o`
#	ranlib lib$(NAME).a  # `ar s` is same as `ranlib`

# create source files' dependency list of header files
depend:
	makedepend -- $(FLAGS) -- $(SOURCES) $(INCLUDES) -s'# DO NOT DELETE THIS LINE -- `makedepend` depends on it.'

# remove object files
clean:
	$(RM) ${OBJECTS}
	$(RM) lib${NAME}.a
	@for test in $(SUBDIR); do \
		make -C $$test clean; \
	done;

# DO NOT DELETE THIS LINE -- `makedepend` depends on it.

