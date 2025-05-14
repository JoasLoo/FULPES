@echo off

g++ -O3 -ffast-math -static-libstdc++ main.cpp -o main.exe

if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)