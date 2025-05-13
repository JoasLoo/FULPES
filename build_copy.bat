@echo off
set PYTHON_DIR=C:\Users\Joas\AppData\Local\Programs\Python\Python313

g++ main_copy.cpp -o main.exe -I%PYTHON_DIR%\include -L%PYTHON_DIR%\libs -lpython313

if %ERRORLEVEL% EQU 0 (
    echo Build successful!
) else (
    echo Build failed!
)