# The smithlab_cpp library

This library contains code that has been used in the Smith lab for too
many years, and that we depend on for several of our data analysis
tools. Many of those tools use older versions of this source code in
subdirectories of other repos.

## Requirements

- A C++ compiler that knows C++17. The GNU `g++` compiler works well
  for this after version GCC9, and the default since GCC11. You can
  get this from apt, conda and brew, but beware on macos where it's
  not really installed by default even though it pretends to be.
- The [zlib library](https://zlib.net), which we use for I/O of files
  in gzip format. You likely have this on your system. You can get
  this from apt, conda and brew, and it's likely installed already.
- Optional: The [HTSLib library](http://htslib.org), which we use for
  I/O of SAM and BAM format files. You can get this from apt, conda
  and brew.

## Building and installing the smithlab_cpp library

Assuming you downloaded the release X.X the tarball
`libsmithlab_cpp-X.X.tar.gz`, you would do the following:
```
$ tar -xvf libsmithlab_cpp-X.X.tar.gz
$ cd libsmithlab_cpp-X.X
$ ./configure
$ make
$ make install
```
If you do not want to contaminate your system's directories with
our code, you can modify the 3rd step above to:
```
$ ./configure --prefix=/some/unimportant/directory
```
If you want to build this code to use our htslib wrapper, you will need
to run like this:
```
$ ./configure --enable-hts
```
You must also have HTSlib installed in some standard place on your
system. If you have it installed in some other place, then you will
need to set variables (CPPFLAGS and LDFLAGS) when running the
configure script.

## Using the source directly from the repo

If you clone the repo and attempt to use the source directly, you are
likely to run into more problems than if you use a release. I will do
everything I can to provide support for the releases, but I may not
help if you have problems using the source repo directly.

This README.md file is written just as we are turning `smithlab_cpp`
into a library and not a collection of source files. If you want to
use it the way it has been used from 2010-2019, then you can use the
`Makefile` in this repo without running the `./configure` script:
```
$ make OptionParser.o
g++ -Wall -std=c++17 -c -o OptionParser.o OptionParser.cpp
```
Note: if you run the `./configure` script it will overwrite the
`Makefile` indicated above. If that happens, just get a new one.  The
`./configure` script must be obtained using autotools if you cloned
this repo; if you downloaded this code as a "release" then the
configure script should be present already.

## TODO

This code needs lots of changes. I'm listing them here for the present
and hope to take care of each with separate issues on GitHub. The
result should be less total code overall in smithlab_cpp.

- `QualityScore.*pp` likely should be removed, as we only use
  sequencing quality scores in specific places, and in those places
  have chosen to re-implement anything that would be here.
- `smithlab_os.*` Any use of character arrays should be replaced with
  strings for filenames. Implementation of many functions in the cpp
  file is sloppy.
- `smithlab_utils.*pp`: lots to replace here. Many functions seem
  redundant with functions elsewhere. Not sure of we need the smithlab
  namespace. Likely the `copy_if` function should be removed. We
  should test of the alphabet conversion functions are all needed.
- There might be redundancy between several functions that span
  `GenomicRegion.*pp`, `chromosome_utils.*pp` and `smithlab_os.*pp`,
  especially in relation to reading files that contain genomes.
