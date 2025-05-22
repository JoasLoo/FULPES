@echo off

<<<<<<< Updated upstream
g++ -O3 -ffast-math -static-libstdc++ -funroll-loops -g -o main main.cpp
=======
g++ -O3 -ffast-math -static-libstdc++ -funroll-loops -o main main.cpp
>>>>>>> Stashed changes


if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)