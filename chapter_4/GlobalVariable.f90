Module GlobalVariable
	Implicit None

    Integer,parameter ::					MyRk=Selected_Real_Kind(6,36), MyIK=Selected_Int_Kind(8)
    Real(kind = MyRk), Parameter :: 		Pi = 3.14159 
    Real(kind = MyRk), Parameter :: 		Miu = 0.5 !// Friction coefficient
    Real(kind = MyRk), Parameter :: 		G = -9.8 !// Acceleration 
    Real(kind = MyRk), Parameter :: 		Rho_Sand = 2650 !// Density of sand
    Real(kind = MyRk), Parameter :: 		Epsl = 0.35
    Real(kind = MyRk), Parameter :: 		Miv = 0.3 !// Poisson ratio
    Real(kind = MyRk), Parameter :: 		Ym = 1.0E7 !// Young modulus 

	Integer( kind=MyIk ) ::					N_Total !// Total number of sand particles
	Integer( kind=MyIk ) ::					N_Sorts !// Nmber of radius sorts
	Integer( kind=MyIk ) ::					N_Incident !// Number of sand incident
	Real( Kind=MyRk ) ::					Height !// Sand bed height 
	Real( Kind=MyRk ) ::					Length !// Sand bed length
	Real( Kind=MyRk ) ::					N_Second !// Total time of simulation
	Real( Kind=MyRk ) ::					DeltaT !// Time step
	Real( Kind=MyRk ) :: 					Edgeright !// Right boder of sand bed 
	Real( Kind=MyRk ) :: 					Edgetop !// Upper boder of sand bed 
	Real( Kind=MyRk ) :: 					MaxRadius !// Maximum of radius
	Real( Kind=MyRk ) :: 					MinRadius !// Minimum of radius
	Real( Kind=MyRk ) :: 					H_SandBed !// Height or sand bed 
	Real(kind = MyRk), Allocatable :: 		Radius(:)



	Type :: Particle
	    Integer( kind=MyIk ) :: 			i !// Particle number
	    Integer( kind=MyIk ) :: 			M !// Row number of particle in mesh 
	    Integer( kind=MyIk ) :: 			N !// Column number of particle in mesh 
	    Real( Kind=MyRk ) :: 				X !// X coordinate
	    Real( kind=Myrk ) :: 				Y !// Y coordinate
	    Real( kind=Myrk ) :: 				Radius !// Radius of particle
	    Real( kind=Myrk ) :: 				Vx !// X direction velocity
	    Real( kind=Myrk ) :: 				Vy !// Y direction velocity
	    Real( kind=Myrk ) :: 				V !// Velocity of particle 
	    Real( kind=Myrk ) :: 				Mass !// Particle Mass
	    Real( kind=Myrk ) :: 				W !// Velocity of particle 
	    Real( kind=Myrk ) :: 				Inertia !//Transshipment inertia

	    Integer( kind=MyIk ) :: 			Co !// Number of collisions
	    Real( kind=Myrk ) :: 				Alpha
	    Real( kind=Myrk ) :: 				Q !// Particle electricity 
	    Real( kind=Myrk ) :: 				Pd !// Electron transfer probability
	End Type Particle
	Type( Particle ), Allocatable ::		Pt(:), Pt_old(:)


End Module GlobalVariable
