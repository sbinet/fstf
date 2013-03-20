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
# drop "svn" suffix from version string
sed -i 's/3\.2svn/3.2/g' ../llvm-src/configure

# configure
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

## clang-fmt

``clang-fmt`` is a backport of ``clang-format`` from the (not yet
released) ``clang-3.3`` version.

It is available from [clang-fmt](http://github.com/sbinet/clang-fmt)

### installation

```sh
# hwaf usual dance
$ hwaf init work
$ hwaf setup work
$ cd work

# real install
$ hwaf pkg co git://github.com/sbinet/clang-fmt
$ hwaf configure
$ hwaf
```

### example

```sh
$ cat > foo.cxx
namespace Ns {class bar{
int m_int    ; public:
bar ( ) ; void some_meth(const int & an_int ) ;};} //> namespace Ns

$ clang-fmt --style=Google foo.cxx
namespace Ns {
class bar {
  int m_int;
 public:
  bar();
  void some_meth(const int& an_int);
};
}  //> namespace Ns
```
