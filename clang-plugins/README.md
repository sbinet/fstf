clang-plugins
=============

A repository of interesting ``CLang`` based plugins.

The first thing is to get a complete ``LLVM+CLang`` installation:

```sh
#!/bin/sh

set -e

svn co http://llvm.org/svn/llvm-project/llvm/branches/release_32 llvm-src
cd llvm-src
cd tools
svn co http://llvm.org/svn/llvm-project/cfe/branches/release_32 clang
cd ../../llvm-src/projects
svn co http://llvm.org/svn/llvm-project/compiler-rt/branches/release_32 compiler-rt
cd ../../
mkdir llvm
cd llvm
../llvm-src/configure \
   --prefix=/some/path \
   --libdir=/some/path/lib \
   --enable-shared \
   --enable-libffi \
   --enable-targets=all \
   --enable-optimized \
   --disable-expensive-checks \
   --disable-debug-runtime \
   --disable-assertions
make REQUIRES_RTTI=1 -j8
make install
```
