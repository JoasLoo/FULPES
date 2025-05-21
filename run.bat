@echo off
.\main.exe

if %ERRORLEVEL% EQU 0 (
    echo Run successful!
) else (
    echo Restarting...
    run.bat
)
