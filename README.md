# ws_img
Library for reading, writing, and processing pnm (PBM, PGM, and PPM) and JPEG images. Includes functions for linear algebra and computer vision.

Availability of a UNIX or GNU/Linux system is assumed.
Here is a brief description for each of these files:
	README -- This file
	go -- Compiles and runs program
	Image.cpp -- Source file for Image class; contains I/O a few functions
	Image.hpp -- Header file for Image class; contains processing functions 
	Vector.hpp -- Class representing a vector
	Matrix.{h,c}pp -- Class representing a matrix
	Linear_Algebra.{h,c}pp -- Functions for solving linear systems
	makefile -- Rules and dependencies for compiling image library
	example/example.cpp -- Example program using the Image class
	example/makefile -- Rules and dependencies for compiling example program
	Doxyfile -- Doxygen configuration file for generating documentation

Before compiling anything, type 'make depend'. This only has to be done once,
but must be done. Many supposed error messages will appear, but this is okay.

The class was developed on a GNU/Linux system, so it should compile on any
UNIX-like system with GCC installed. For compiling in UNIX/Linux, do the
following:
	- To compile the code, type 'make'.
	- To compile and run the code, type 'go'.
To change the images used, change the lines 'in_img=in.pgm' and
'out_img=out.pgm' in the files go to the input and output images one wishes to
use. To change the name of the source file used, such as to modify it for
another use, change the following lines in makefile to the new values:
	NAME = example
Also, in go, change the word 'com=example' to the new name being used.

To use the Image class, it is easiest to simply refer to the file example.cpp
for an example of how it is used in real code. Creating Image objects, reading
and writing PGM files, and accessing individual pixels are all shown in this
file. For compiling, refer to 'example/makefile', especially the following
lines:
	TOOLS_DIR = /home/wade/cpp/tools
	IMAGE_DIR = /home/wade/cpp/vision/Image

	INCLUDES += -I$(TOOLS_DIR)
	INCLUDES += -I$(IMAGE_DIR)

	LINK_DIRS += -L$(TOOLS_DIR)
	LINK_DIRS += -L$(IMAGE_DIR)

	LIBS += -ltools
	LIBS += -limage
These lines give the names and locations of the libraries and headers that must
be linked against during compilation.

Image.{cpp,hpp}, Matrix.{cpp,hpp}, etc. are the actual implementation files for
the Image class and are of limited use to a user.  

The file Doxyfile is a configuration file used by the program doxygen to
generate html documentation for the Image class. If doxygen is installed,
simply type 'doxygen' to create a directory html/ containing extensive
documentation and point one's browser to html/index.html to view it.
