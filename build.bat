@echo off

g++ -fopenmp -O2 -ffast-math -static-libstdc++ -funroll-loops -g -o main main.cpp


if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)