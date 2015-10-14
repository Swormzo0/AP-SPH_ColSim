@echo off
rem set OPTIONS=1 1 300 5 0 0 0
rem set OPTIONS=1 2 180 -.2 0 1 0
rem OPTIONS=2 1 0 1 0 -.2 0
set OPTIONS=5 2 330 3 0 -7 2
echo %OPTIONS% | AP_SPH_CollisionSimRoughVDebug > 2DMoving.txt
