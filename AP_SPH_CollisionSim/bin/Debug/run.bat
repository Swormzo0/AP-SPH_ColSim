@echo off
rem set OPTIONS=1 1 300 5 0 0 0
rem set OPTIONS=1 2 180 -.2 0 1 0
set OPTIONS=2 1 0 1 0 -.2 0
rem set OPTIONS=8 5 307.75 45 0 0 0
echo %OPTIONS% | AP_SPH_CollisionSim > log.xls
rem echo %OPTIONS% | AP_SPH_CollisionSimQuad > logQ.xls
