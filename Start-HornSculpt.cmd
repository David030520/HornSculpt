@echo off
setlocal
set "hornsculpt_program=%~dp0bin\HornSculpt.exe"
if not exist "%hornsculpt_program%" if exist "%~dp0out\launch-path.txt" set /p "hornsculpt_program="<"%~dp0out\launch-path.txt"
if not exist "%hornsculpt_program%" (
    echo HornSculpt.exe wurde nicht gefunden.
    echo Bitte zuerst das Anwendungspaket installieren.
    pause
    exit /b 1
)
start "" "%hornsculpt_program%" %*
