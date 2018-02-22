# Makefile for particle filling and sedimention
# 2017/07/12
# Mu LiAn LanZhou University



.PHONY: clean

main.exe: Accumulation.o Main.o GlobalVariable.o MeshGrid.o FileOperate.o Sedimentation.o ScheduleBar.o Distribution.o
	ifort GlobalVariable.o MeshGrid.o FileOperate.o Accumulation.o  Sedimentation.o  ScheduleBar.o Distribution.o Main.o -o main.exe
Main.o: Accumulation.o Sedimentation.o FileOperate.o
	ifort -c Main.f90


Accumulation.o: GlobalVariable.o MeshGrid.o FileOperate.o ScheduleBar.o Accumulation.f90 Distribution.o
	ifort -c Accumulation.f90


Sedimentation.o: FileOperate.o Sedimentation.f90
	ifort -c Sedimentation.f90


FileOperate.o: GlobalVariable.o FileOperate.f90
	ifort -c FileOperate.f90


GlobalVariable.o: GlobalVariable.f90
	ifort -c GlobalVariable.f90
	

MeshGrid.o: MeshGrid.f90
	ifort -c MeshGrid.f90


ScheduleBar.o: ScheduleBar.f90
	ifort -c ScheduleBar.f90

Distribution.o: Distribution.f90 GlobalVariable.o
	ifort -c Distribution.f90



clean:
	-rm *.mod *.o
