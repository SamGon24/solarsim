@echo off
rem ============================================================================
rem  SolarSim - one-shot build script for Windows
rem
rem    build.bat            Release build
rem    build.bat Debug      Debug build
rem
rem  Tries, in order: CMake (Visual Studio or MinGW generator), the MSVC
rem  compiler cl.exe (run from a "Developer Command Prompt for VS"), and
rem  MinGW-w64 g++ (e.g. from MSYS2 or WinLibs) on the PATH.
rem ============================================================================
setlocal EnableDelayedExpansion
cd /d "%~dp0"

set MODE=%1
if "%MODE%"=="" set MODE=Release

set SRCS=
for %%d in (core astro physics analysis render io platform app) do (
  for %%f in (src\%%d\*.cpp) do set SRCS=!SRCS! %%f
)
set SRCS=!SRCS! src\main.cpp

where cmake >nul 2>nul
if %ERRORLEVEL%==0 (
  echo ==^> CMake build ^(%MODE%^)
  cmake -S . -B build -DCMAKE_BUILD_TYPE=%MODE%
  if errorlevel 1 goto try_cl
  cmake --build build --config %MODE% --parallel
  if errorlevel 1 goto try_cl
  if exist build\%MODE%\solarsim.exe (set BIN=build\%MODE%\solarsim.exe) else (set BIN=build\solarsim.exe)
  goto done
)

:try_cl
where cl >nul 2>nul
if %ERRORLEVEL%==0 (
  echo ==^> MSVC cl.exe build
  if not exist build-cl mkdir build-cl
  set CLOPT=/O2
  if /I "%MODE%"=="Debug" set CLOPT=/Od /Zi
  cl /nologo /std:c++17 !CLOPT! /EHsc /utf-8 /W3 /DNOMINMAX /D_CRT_SECURE_NO_WARNINGS /Isrc /Fobuild-cl\ !SRCS! /Fe:solarsim.exe
  if errorlevel 1 goto fail
  set BIN=solarsim.exe
  goto done
)

where g++ >nul 2>nul
if %ERRORLEVEL%==0 (
  echo ==^> MinGW g++ build
  set GOPT=-O2
  if /I "%MODE%"=="Debug" set GOPT=-O0 -g
  g++ -std=c++17 !GOPT! -D__USE_MINGW_ANSI_STDIO=1 -Isrc !SRCS! -o solarsim.exe -static -static-libgcc -static-libstdc++ -pthread
  if errorlevel 1 goto fail
  set BIN=solarsim.exe
  goto done
)

echo error: no C++ compiler found.
echo Install Visual Studio 2019+ (Desktop development with C++) and run this from a
echo "Developer Command Prompt", or install MSYS2/MinGW-w64 and put g++ on the PATH.
exit /b 1

:fail
echo Build failed.
exit /b 1

:done
echo.
echo Build complete: %BIN%
echo Try:  %BIN% selftest     %BIN% validate     %BIN% view     %BIN% help
echo Tip: use Windows Terminal for 24-bit colour; in old consoles add --color ascii --ascii
endlocal
