@echo off
REM Windows batch file to run gRNA finder
"%~dp0crispr_env\Scripts\python.exe" "%~dp0run.py" grna %*
