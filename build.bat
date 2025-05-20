@echo off

g++ -fopenmp -O3 -ffast-math -static-libstdc++ -funroll-loops -o main main.cpp


if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)