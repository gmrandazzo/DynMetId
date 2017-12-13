DynMetId
=============

DynMetId an automatic annotation tool for LC-MS metabolomics

Info:
----------
Mail me at: gmrandazzo@gmail.com

License
============

DynMetId is distributed under GPLv3 license, this means that:

- you can use this library where you want doing what you want.
- you can modify this library and commit changes.
- you can not use this library inside a commercial software.

To know more in details how the licens work please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

DynMetId is currently property of Giuseppe Marco Randazzo which is also the
current package maintainer.


Dependencies
============

The required dependencies to use DynMetId are:
- mysql (mariadb, mysql)
- cpp compiler (gcc or clang for osx)
- cmake

Install
=======

Compile from source
-------------------

  cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
  make
  sudo make install

