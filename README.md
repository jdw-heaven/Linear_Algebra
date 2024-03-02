# These files is based on LZU lesson --- Linear_Algebra_2

## sometimes neeeds to use pytorch, libtorch and lapack.

## The steps to cmake:

```
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/absolute/path/to/libtorch ..
cmake --build . --config Release
```