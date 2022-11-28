Gradient Product Transform for Symmetry Detection and Medial Axis Extraction
============================================================================

The Gradient Product Transform (GPT) is an algorithm that assigns a symmetry
score to each pixel in a grayscale image. This score values can be used for
detecting symmetric objects or for extracting medial axes, e.g. from blood
vessels.

This code implements the algorithm described in (cited as "IPOL paper" below):

> C. Dalitz, J. Wilberg: "The Gradient Product Transform:
> An Image Filter for Symmetry Detection."  Image Processing On Line, 9
> (2019), pp. 413–431
> https://doi.org/10.5201/ipol.2019.270

Please cite this article when using the code. The article was published on
IPOL with version 1.1 of this code. For changes since then, see the file
*CHANGES*.

Compilation
-----------

Building the code requires cmake and a standard C++98 (or later) compiler.
We have tested the code with gcc 5.4.0, LLVM 9.0.0, and MSVC 14.1.

The number of kernels used for parallel computation can be set by the
variable NUM_THREADS in CMakeLists.txt. Parallelization requires OpenMP
to be installed. While OpenMP is supported by default on Linux, it must
must be installed manually on OSX with "brew install libomp". Moreover,
at least cmake version 3.12 is required for OpenMP support on OSX and for
MSVC, whilst on other platforms cmake 3.5 is sufficient.

Before running cmake, you should adjust the variable NUM_THREADS in
CMakeLists.txt to the number of kernels on your system.

Starting from the root directory (i.e., the directory, in which this
Readme file is located), the executable "gptsymmetry" is created with:

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make

To test the program, you can apply it to the test image included in this
source package as follows:

    gptsymmetry -r 80 -o result.png ../testimage.png

This will mark the detected symmetries in the result image file "result.png".


Usage
-----

The command line program *gptsymmetry* reads a PNG image, prints the detected
symmetries to stdout and optionally writes images showing a skeleton of
symmetry axes and of the detected symmetric objects.

Calling *gptsymmetry* without any or with an unknown option (e.g. "-?")
will print a usage message. The meaning of the parameters controlling
the algorithm is explained in the IPOL paper.

The input file must be in PNG format. With the option "-o <outfile>", the
detected symmetry regions are darwn into the image and written to <outfile>.
Unless the option "-notrace" is given, three additional images are written:

 - ``trace-symmetry.png``  
   the GPT symmetry score converted to grayscale
 - ``trace-ridges.png``  
    skeletons of symmetry axes
 - ``trace-localmaxpoints.png``
    the symmetry regions classified into axial (red) an rotational (cyan)

The detected symmetry regions are printed to stdout in the following format:

    rot;x;y;rx;ry;S;s_norm
    0;142;317;5;23;5044167.206566;0.903770
    0;141;325;5;14;4800286.097710;0.947354
    1;295;198;32;33;4569864.069128;0.784917
    0;106;298;38;42;3361987.012759;0.720138

where "rot" indicates rotational (1) or axial (0) symmetry, (x,y) is the
center and (rx,ry) the size of the symmetry region, S is the GPT score,
and s_norm is the normalized GPT score between zero and one. Thus, to report
only rotational symmetries, you can pipe the output to grep as follows:

    gptsymmetry -r 100 -notrace test.png | grep '^1'


Source Files
------------

 - ``png++/``
   this is a copy of libpng++ version 0.2.9 from https://www.nongnu.org/pngpp/
   with two modifications to error.hpp:
    - fix for bug #46312 applied
	- strerror_r redifened as strerror_s for MSVC compiler

 - ``main.cpp``
   Main program that implmenting Algorithms 1-3
   (sections 2 & 3 in the IPOL paper)

point.h
   Implements types and function related to 2D points.

image.[h|cpp]
   Implements three image types (Gray, RGB, Float)

image_processing.[h|cpp]

symmetry_transform.[h|cpp]
   Implementats functions for individual steps in Algorithms 1-3.



Authors & Copyright
-------------------

Christoph Dalitz, Jens Wilberg, 2019
Institute for Pattern Recognition
Niederrhein University of Applied Sciences
Krefeld, Germany

The software includes the library png++ by Alexander Shulgin, available
from https://www.nongnu.org/pngpp/. See the directory src/png++ for details.
