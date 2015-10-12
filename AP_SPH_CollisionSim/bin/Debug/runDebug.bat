@echo off
rem set OPTIONS=1 1 300 5 0 0 0
rem set OPTIONS=1 2 180 -.2 0 1 0
rem set OPTIONS=2 1 0 1 0 -.2 0
set OPTIONS=8 5 308 4.5 0 0 0
echo 0 %OPTIONS% | AP_SPH_CollisionSim > logA.txt
echo 1 %OPTIONS% | AP_SPH_CollisionSim > logB.txt

