@echo off

g++ -O3 -ffast-math -static-libstdc++ -gcflags "all=-m -m -d=ssa/check_bce/debug" -funroll-loops -o main main.cpp

if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)