# Futilities: unsorted Fortran utility modules

---

### Brief summary of their purpose

* bspline: Multidimensional spline interpolations (check source file
  for original author).

* cuba: Fortran interface to the [**cuba**](https://feynarts.de/cuba/)
  multidimensional integrator.

* eomrelax & sorelax: Successive over-relaxation integrator of partial
  differential equations (eomrelax encodes the equations of motion,
  sorelax is the relaxator).

* fifo: First In First Out memory buffer, you can empty of fill the
  queue from the top or from the bottom (used by the scheduler module).

* iofifo: Saving and reading data for the fifo module.

* scheduler: A basic Message Passing Interface (MPI) queue
  scheduler. It distributes an equal list of jobs to various MPI
  processes. When some of these processes have finished their load, they
  can start stealing work to the others.

* index: A little module to shuffle array indexing.

* iopara: A module to wrap-up some parallel writting function of the
  Message Passing Interface (MPI) library.

* flatten: Flattening a multidimensional array into 1D and reverting the operation.

* fwcs: A Fortran interface to the World Coordinate System library [**wcslib**](https://www.atnf.csiro.au/people/mcalabre/WCS/)

* wcswrap.c: A C-wrapper to some functions of the **wcslib**, for the fwcs module


* iofits: A Fortran interface to some functions of the [**fitsio**](https://heasarc.gsfc.nasa.gov/fitsio/) library.

* iotools: Some wrappers to do debugging and simple I/O in Fortran.

* linalg: A linear algebra module using the [**lapack**](https://netlib.org/lapack/) library

* linfft: One dimensional Fourier transform based on the [**fftw**](http://fftw.org/) library

* lm: A non-linear fitting module, using the Levenberg-Marquardt method (core from netlib, BSD-license)

* sepbvp: A module to solving boundary value problems, encapsulating the Cash & Mazzia code (see source file).

* sundials: A module to call the differential equation solvers of the [**sundials**](https://computing.llnl.gov/projects/sundials) library.

---

### Licensing

Some modules encapsulate some public code and they inherit their
licence from their ancestor. For the others, this is GPLv3, even if not
mentioned in the source file yet.