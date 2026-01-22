@echo off
REM Windows batch file to run indel analysis
"%~dp0crispr_env\Scripts\python.exe" "%~dp0run.py" analysis %*
