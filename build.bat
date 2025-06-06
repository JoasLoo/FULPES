@echo off

g++ -O3 -ffast-math -funroll-loops -fprofile-generate -o main main.cpp
main
g++ -O3 -ffast-math -funroll-loops -fprofile-use -o main main.cpp
main

if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)