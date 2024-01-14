# Polynomial Operations Library

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

A C++ library, that supports basic and advanced operations on polynomials. It also supports polynomials in the Z/pZ ring.

## Features

- **Type support**
  - You can choose the type of coefficients yourself (pass it as a template parameter)
  - The library also supports my ModInteger class (refer to https://github.com/vgzowski/Studying/tree/main/Modular), and allows you to do absolutely precise calculations with integers modulo some number

- **Basic Operations:**
  - Addition (`+`)
  - Subtraction (`-`)
  - Multiplication (`*`)
  - Division (`/`)
  - Input / output ('>>', '<<' operators)

- **Specialized Operations:**
  - Derivative of a polynomial
  - Integral of a polynomial
  - Multiplying by `x^k`
  - Logarithm of a polynomial
  - Exponential of a polynomial
  - Interpolation
  - Multi-point evaluation

- **Advanced Functions:**
  - Power function for polynomials
  - Square root function for polynomials
 
- **Requirements**
  - In order for this library to work, you must also include, and if needed, change the path in Polynomial.h file, the FFT library
  - Please refer to https://github.com/vgzowski/Studying/tree/main/FFT
 
- **Usage**
  - Usage can be found at /src/testing.cpp
