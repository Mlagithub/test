
.PHONE: clean


bj.exe: GlobalVariable.o FileOperate.o MeshGrid.o Collision.o Main.o Electricity.o Evolution.o ScheduleBar.o
	ifort  GlobalVariable.o FileOperate.o MeshGrid.o Collision.o Electricity.o Evolution.o ScheduleBar.o Main.o -o bj.exe 

Main.o: Evolution.o 
	ifort  -c Main.f90 

Evolution.o: GlobalVariable.o FileOperate.o MeshGrid.o Collision.o ScheduleBar.o
	ifort  -c Evolution.f90

GlobalVariable.o: 
	ifort  -c GlobalVariable.f90 

FileOperate.o: GlobalVariable.o Electricity.o
	ifort  -c FileOperate.f90

Electricity.o: GlobalVariable.o
	ifort  -c Electricity.f90

MeshGrid.o: GlobalVariable.o
	ifort  -c MeshGrid.f90 

Collision.o: GlobalVariable.o
	ifort  -c Collision.f90

ScheduleBar.o:
	ifort  -c ScheduleBar.f90




clean:
	rm  *.mod *.o
