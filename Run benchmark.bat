@echo off
cd /d "%~dp0"
where docker >nul 2>&1
if errorlevel 1 (
    echo Install and start Docker Desktop, then run this file again.
    pause
    exit /b 1
)
docker info >nul 2>&1
if errorlevel 1 (
    echo Start Docker Desktop and wait for it to be ready.
    pause
    exit /b 1
)
docker compose run --rm --build benchmark
if errorlevel 1 (
    echo Reproduction failed. Review the message above.
    pause
    exit /b 1
)
echo Results are in the reproduced-results folder.
pause
