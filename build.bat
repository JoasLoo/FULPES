@echo off


g++ -O3 -ffast-math -static-libstdc++ -funroll-loops -fprofile-generate -o main main.cpp
main
g++ -O3 -ffast-math -static-libstdc++ -funroll-loops -fprofile-use -o main main.cpp


if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)