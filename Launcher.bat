@echo off
chdir Projects
call conda activate base
python master.py
pause