!// Global variable 
Module GlobalVariable
	Implicit None


	! Physical constants
	Integer,Parameter :: 					Myrealk=4, Myintk=4
	Real( kind=Myrealk ),Parameter :: 		Pi=3.14159, G=9.8


	! Variables
	Integer :: 								N_Total ! Total number of particles
	Integer :: 								N_Sorts ! Number of particle sorts
	Real( kind=Myrealk ) :: 				Accuracy ! Interate accuracy
	Real( kind=Myrealk ) :: 				Bl ! The left border of the area
	Real( kind=Myrealk ) :: 				Br ! The right border of the area
	Real( kind=Myrealk ) :: 				Bb ! The lower border of the area
	Real( kind=Myrealk ) :: 				Bt ! The upper border of the area
	Real( kind=Myrealk ) :: 				MaxRadius, MinRadius ! The maximum and minimum of the particle radius
	Real( kind=Myrealk ) :: 				Rho ! Particle density
	Real( kind=Myrealk ) :: 				Ym ! Young's Modulus
	Real( kind=Myrealk ) :: 				Miv ! Poisson ratio
	Real( kind=Myrealk ), Allocatable :: 	Radius(:), Ratio(:) ! Particle radius and ratio of radius 


	! Define a type which describes the physical properties of the particle
	Type :: Particle
		Integer ::							i ! Particle ID
		Integer ::							M, N ! Particle mesh location
		Real( kind=Myrealk ) ::				X ! Particle X corrdinate
		Real( kind=Myrealk ) ::				Y ! Particle Y corrdinate
		Real( kind=Myrealk ) ::				Radius ! Particle radius
		Real( kind=Myrealk ) ::				Mass ! Particle mass
	End Type
	Type( Particle ), Allocatable :: 		Pt(:) ! Array of type Particle


End Module GlobalVariable
