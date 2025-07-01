# Time-Dependent Schrödinger Equation Solver (TDSE)

This project numerically solves the Time-Dependent Schrödinger Equation (TDSE) using C++ and ArrayFire, with GPU acceleration.

## 🧠 Features

* Solves TDSE in 1D for a given potential
* GPU acceleration with ArrayFire
* Customizable parameters for wavefunctions and potentials

## 📸 Screenshots

(Add visual output if applicable, like wavefunction plots or gifs)

## ⚙️ Dependencies

- [ArrayFire](https://arrayfire.com/) (Install via package manager or from source)
- CMake ≥ 3.10
- A C++17 compiler
- CUDA (if using GPU acceleration)
- `json.hpp` (included in `/include/json/`)

## 🛠️ Build Instructions

```bash
git clone https://github.com/ahwinney/TDSE.git
cd TDSE
mkdir build && cd build
cmake ..
make
```


## 🔌 Dependencies

You’ll need:

- [ArrayFire](https://arrayfire.com/) (Install via package manager or from source)
- CMake ≥ 3.10
- A C++17 compiler
- CUDA (if using GPU acceleration)
- `json.hpp` (included in `/include/json/`)


