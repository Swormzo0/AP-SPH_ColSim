@echo off
set OPTIONS=2 6 0 12 0 0 0
rem set OPTIONS=1 2 180 -.2 0 1 0
rem set OPTIONS=2 1 0 1 0 -.2 0
rem set OPTIONS=8 5 307.75 4.5 0 0 0
echo %OPTIONS% | AP_SPH_CollisionSim 
pause
rem > Linear.xls

