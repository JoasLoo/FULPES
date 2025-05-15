@echo off

g++ -O3 -pg -g -ffast-math -static-libstdc++ -funroll-loops -o main main.cpp
main
gprof main gmon.out > analysis.txt


if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)