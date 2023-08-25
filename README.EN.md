# Upper atmosphere GOST R 25645.166-2004

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[Документация на русском](README.md)

![under_construction](under_construction.png)

*(In progress, not tested yet!)* Function to calculate upper Earth atmosphere density by GOST R 25645.166-2004 model in C++, Python and Matlab. GOST R 25645.166-2004 model is russian equivalent to NRLMSISE-00, Jacchia-Bowman JB-2008 models and others.

## Use in C++ projects

1) It's a header-only library. Just place `atmosGOST_R_25645_166_2004.h` and `atmosGOST_R_25645_166_2004.cpp` somewhere into your project directory
2) `#include "atmosGOST_R_25645_166_2004.h"`
3) Calculate density by calling `atmosGOST_R_25645_166_2004()` function with appropriate parameters (look comments to function in atmosGOST_R_25645_166_2004.cpp or into GOST R 25645.166-2004 document itself).

## TODO list

- [ ] Python and Matlab functions.
- [ ] Performance optimization.
- [ ] Documentation.
- [ ] Internationalization.
- [ ] CMake infrastructure.
- [ ] Tests.
