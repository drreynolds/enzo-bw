Dear User,

Here is a self-contained tarball with the code and makefiles for
Kraken, Jaguar and BlueDrop.  I'll add the others as necessary.  There
are many changes that I have yet to make to this code but I guess I'll
do them after you have made your changes and we have a stable code
that can be used for the INCITE runs. Unfortunately, this will result
in 2 disparate ENZO codes. But we are under the gun to get a big run
going ASAP and I am happy that you want to work with this code as I
think it has been well tested.

So...

The tar file contains a single Makefile and the machine-dependent files for
Kraken (PGI and Intel), Jaguar (PGI) and BlueDrop (IBM XL).  You shouldn't
need to change more than a couple of the Makefile parameters.

Some notes:

There are bugs in PGI and Intel which affect hybrid ENZO.
With PGI you need to use 10.3.0 or EARLIER.

On Kraken,

module swap pgi/10.6.0 pgi/10.3.0

On Jaguar,

pgi/10.3.0 should be default.

With Intel, either you need Intel V12.1 or you will need to
compile these two routines WITHOUT OpenMP:

OMP_NEW_multi_cool.src
OMP_solve_rate_cool.src

The Intel 11.1 compiler will suffer an internal error otherwise!

The directories CRAY_MPT* each contain two files:

ENZO_mpicxx.h
mpi.h


You will need to move one set into the source directory if you 
build the code on Kraken or Jaguar.  Remove these from the source
directory when using other platforms (or the compilation will die).
These files get around a bug in the Cray C++ wrappers.
These are versions for MPT 4.x and 5.x, but either will work.

I have maintained the source with and without OpenMP.
This isn't necessary any more, so please make your changes to
the routines with the OMP_ prefix.  I'll do the back-stitch to
the non-OMP if we ever need it.  Mainly, I preserved this option
because the various OpenMP compilers are flaky.  Beware the Intel
compiler in particular (and note the difference between Intel OpenMP
syntax and the standard supported by IBM, Cray and PGI!). 
Note that there are also C syntax differences for handling function
return codes in OpenMP.  I use a shared vector to gather the individual
return codes from the C routines.  These need checks on the return code
vectors and I I will fix those up later.  (BTW, I have a replacement for 
Matt's / David's macros for error handling in non-OpenMP situations.
I think mine is better because it allows for a choice in behaviour
rather than a specific action).

Possibly the biggest code difference lies in the initialization which
requires the 3D ICs to be decomposed prior to running ENZO.
This is necessary when using more than a few thousand tasks. 
There are also some details in the I/O for large core counts
on Lustre file systems (i.e. Kraken, but may be applicable to Jaguar).

(BTW, I dislike the "kitchen sink" philosophy in ENZO - my preference
would be for separating test and analysis codes from the main code and
instead linking these against a core "enzo.a" as necessary).

To avoid the dangerous macro definitions in the old typedefs.h,
I MANUALLY decide on the question of whether fields are density-like
and so on, and set this with the assignment of the field id number.

SetFieldDensityType.C
SetFieldColorType.C
PrintFieldDensityType.C
PrintFieldColorType.C

This isn't 100% complete: there is still some unsafe code in the
various Grid_SolveHydroEquations where assumptions are made about
the sequential assignment of the Baryon fields in the intitialization code.
As you will be adding more, perhaps we should fix this - its simple to do.

I want to preserve the three (four?) hydro versions - if we
have to choose one then I would choose the original Fortran version
since it is the fastest and uses the least memory, does the
fewest mallocs, etc. and parallelizes well. (And I have never had
any trouble with it at all on 64-bit or 32-bit computers).

I use RAD_HYDRO as the macro to isolate the RHD code, as you did
in your original code.  It is not clear to me why you abandoned
that in E2.0 and (apparently) mixed it up with JW's ray stuff.
If anything, I think ENZO needs _more_ separation of the various
physics packages, not less!

Since you will be modifying the ENZO chemistry, please take
note of the OpenMP critical regions (the routines with
OMP_CRITICAL_) where critical regions are necessary because
Pascal (I think) abused the (global) contents of .h files.
That was the hardest thing to detect and fix because it caused
race conditions.

Can we keep the mass units as they are in my version?
I would just prefer to use something we *know* works...

I am sure there are other things lurking.  Probably more bugs.
I have tests for 32-bit unigrid cosmology, 64-bit AMR, 64-bit RHD, etc.
I also have an 800^3 14 Mpc/h RHD test which I have run on Kraken
and I am setting up 1600^3 28 Mpc/h and 3200^3 56 Mpc/h (the actual
target mesh). 

No doubt testing will show up some interesting "issues"...

Regards,
Robert Harkness
