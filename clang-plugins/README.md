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

## clang-cms

``clang-cms`` is a plugin from the CMS experiment to check for various
multi-thread codesmells.
It is available from
[Twiki](https://twiki.cern.ch/twiki/bin/view/Main/ClangCms) and has
been hwaf-ified over [there](http://github.com/sbinet/clang-cms)

### installation

```sh
# hwaf usual dance
$ hwaf init work
$ hwaf setup work
$ cd work

# real install
$ hwaf pkg co -b=wip/hwaf git://github.com/sbinet/clang-cms
$ hwaf configure
$ hwaf

# test everything
$ hwaf shell
[hwaf] $ cd src/clang-cms/test
[hwaf] $ scan-build -load-plugin libClangCms.so -enable-checker threadsafety make -B
scan-build: Using '/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/clang' for static analysis
/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/c++-analyzer mutable_member.cpp
/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/c++-analyzer const_cast.cpp
/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/c++-analyzer const_cast_away.cpp
/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/c++-analyzer global_static.cpp
/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/c++-analyzer static_local.cpp
/afs/cern.ch/work/b/binet/dev/clang-plugins/llvm/3.2/x86_64-slc6-gcc47-opt/bin/c++-analyzer global_static_edm.cpp
global_static_edm.cpp:20:1: warning: Non-const variable 'g_static' is static and might be thread-unsafe
static int g_static;
^
const_cast.cpp:12:23: warning: const_cast was used, this may result in thread-unsafe code
global_static_edm.cpp:23:1: warning: Non-const variable 'g_static_edm_namespace' is static and might be thread-unsafe
    std::string & r = const_cast< std::string & >( r_const );
                      ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static edm::InputSourcePluginFactory g_static_edm_namespace;
^
1 warning generated.
2 warnings generated.
```


## clang-refactor

``clang-refactor`` is a plugin to seamlessly refactor ``C++`` code at
the ``AST`` level

It is available from
[github](http://github.com/lukhnos/refactorial) and has
been hwaf-ified over [there](http://github.com/sbinet/clang-refactor)

### installation

```sh
# hwaf usual dance
$ hwaf init work
$ hwaf setup work
$ cd work

# real install
$ hwaf pkg co -b=wip/hwaf git://github.com/sbinet/clang-refactor clang-refactor
$ hwaf configure
$ hwaf

# test
$ hwaf shell
[hwaf] $ cd src/clang-refactor/tests/DeclTypeRename
[hwaf] $ ./test.sh

```
-- Configuring done
-- Generating done
-- Build files have been written to: /home/binet/tmp/clang/work/src/clang-refactor/tests/DeclTypeRename
Scanning dependencies of target foo
[100%] Building CXX object CMakeFiles/foo.dir/foo.cpp.o
Linking CXX executable foo
[100%] Built target foo
Running CMake to regenerate build system...
-- Configuring done
-- Generating done
-- Build files have been written to: /home/binet/tmp/clang/work/src/clang-refactor/tests/DeclTypeRename
Warning: No files selected. Operating on all files.
TypeRenameTransform
Processing: /home/binet/tmp/clang/work/src/clang-refactor/tests/DeclTypeRename/foo.cpp.
Ignoring: /usr/include/.*
renames: class A::Foo -> Foobar
renames: class A::(\w) -> \1\1
--- foo.orig.h	2013-03-20 13:46:03.000000000 +0100
+++ foo.h	2013-03-20 13:53:39.000000000 +0100
@@ -9,34 +9,34 @@
     int x;
   };
     
-  class P {};
-  class X : virtual P {};
-  class Y : virtual P {};
-  class Foo : public X, virtual Y {
+  class PP {};
+  class XX : virtual PP {};
+  class YY : virtual PP {};
+  class Foobar : public XX, virtual YY {
   protected:
     int x;
   public:
-    Foo() : x(0) {}
-    Foo(int px) : x(px) {}
+    Foobar() : x(0) {}
+    Foobar(int px) : x(px) {}
     
-     ~Foo ()  {  }
+     ~Foobar ()  {  }
   };  

[...]
-Foo foo3(Foo a, Foo* b, Foo& c, Foo***)
+Foobar foo3(Foobar a, Foobar* b, Foobar& c, Foobar***)
 {
-  return Foo(5);
+  return Foobar(5);
 }
Scanning dependencies of target foo
[100%] Building CXX object CMakeFiles/foo.dir/foo.cpp.o
Linking CXX executable foo
[100%] Built target foo
