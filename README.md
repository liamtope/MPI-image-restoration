# MPI-image-restoration
*Restoration of an image which has been passed through an edge detection algorithm. Part of the EPCC Archer MPI course.*


This repository contains the source code files for the Archer MPI course, reconstructing an image from an edge-detected 
input image. The source files are the following:

  *  ```image.f90```: Main program for reconstructing the image
  *  ```mp_image.f90```: Module containing the MPI routines called in the main program
  *  ```decomp_1d.f90```: Module containing subroutines for 1D domain decomposition
  *  ```decomp_2d.f90```: Module containing subroutines for 2D domain decomposition
  *  ```pgmio.f90```: Module containing subroutines for handling ```.pgm``` file I/O (provided by EPCC)

When running the program, the file name of the input image should be used as a command line argument. No other 
command line arguments are required.
