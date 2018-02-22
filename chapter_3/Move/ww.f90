!//--------------------------------------------------------------------//
!// This module is the main algorithm of particle accumulation.
!//		To generate a particle accumulation with specific parameter
!//		you should call ParticleFilling().
!//--------------------------------------------------------------------//

Module Accumulation
	Use GlobalVariable
	Use MeshGrid
	Use FileOperate
	Implicit None

Contains 

!// Call this subroutine to generate a particle accumulation
Subroutine ParticleFilling( )
	Use ScheduleBar
	Implicit None
	Integer ::					i, Njc, M, N, Mi, Err
	Real( kind=Myrealk ) :: 	Top, Rate
	Character( len=512 ) ::		FileName
	Type( CLS_CMD_Progress ) :: Progress
	Character( len=100 ) ::		A(100)


	!// Initialize the variable
	Call SetVariable()


	!// Initialize schedual bar
	Call Progress%Set( N=N_Total,L=30 )
	Progress%Prefix = "Accumulate Particle: "
	Progress%m = "#"
	Progress%o = "."

	!// Save program parameters to configuration file
	Call SaveConfigure()

	
	!// Divide the grid
	MeshSize = 2.0*MaxRadius
	Call Get_Mesh_2D( N_Total,Bl,Br,Bb,Bt,MinRadius,MeshSize )


	!// Two methods of accumulate
	A(1) = '|  Please select the label of filling method'
	A(2) = '|      1>    Falling from a point'
	A(3) = '|      2>    Slide from the edge'
	Call ShowInformation( A(1:3),3,3 )
	Do 
		Read( *,*,Iostat=Err ) Mi
		If( Err/=0 .or. Mi<0 .or. Mi>2 )then
			Print*,'Try again'
		Else
			Exit
		End If
	End Do 


	Select Case( Mi )
	!// First method
	Case( 1 )
		Top = 2.0*MaxRadius
		Do i = 1, N_Total
			Call GetPosition( Top,i,Mi )
			call Progress%Put( i,.True. )

			! Move particle left and right
			Call Down( i,Njc ) !First

			! Rotate particle along particle Njc
			If( Njc>N_Total)then
				If( Top<Pt(i)%Y+Pt(i)%Radius ) Top = Pt(i)%Y+Pt(i)%Radius
				Call MeshUpdate( i,Pt(i)%X,Pt(i)%Y,M,N,Bl,Bb )
				Pt(i)%M = M
				Pt(i)%N = N
			Elseif( Njc>0 )then
				Call Rotate( i,Njc )
				If( Top<Pt(i)%Y+Pt(i)%Radius ) Top = Pt(i)%Y+Pt(i)%Radius
				Call MeshUpdate( i,Pt(i)%X,Pt(i)%Y,M,N,Bl,Bb )
				Pt(i)%M = M
				Pt(i)%N = N
			Else
				Print*,'Something maybe wrong'
			End If
		End Do 
		

	!// Second method
	Case( 2 )
		Top = 2.0*MaxRadius
		Pt(1)%X = Pt(1)%Radius + Bl
		Pt(1)%Y = Pt(1)%Radius + Bb
		Call MeshUpdate( 1,Pt(1)%X,Pt(1)%Y,Pt(1)%M,Pt(1)%N,Bl,Bb )
		Do i = 2, N_Total
			Call GetPosition( Top,i,Mi )
			call Progress%Put( i,.True. )

			! Move particle left and right
			Call Down( i,Njc ) !First
			Call Left( i,Njc ) !Second

			! Rotate particle along particle Njc
			If( Njc>N_Total)then
				If( Top<Pt(i)%Y+Pt(i)%Radius ) Top = Pt(i)%Y+Pt(i)%Radius
				Call MeshUpdate( i,Pt(i)%X,Pt(i)%Y,M,N,Bl,Bb )
				Pt(i)%M = M
				Pt(i)%N = N
			Elseif( Njc>0 )then
				Call Rotate1( i,Njc )
				If( Top<Pt(i)%Y+Pt(i)%Radius ) Top = Pt(i)%Y+Pt(i)%Radius
				Call MeshUpdate( i,Pt(i)%X,Pt(i)%Y,M,N,Bl,Bb )
				Pt(i)%M = M
				Pt(i)%N = N
			Else
				Print*,'Something maybe wrong'
			End If
		End Do 

	End Select


	!// Calculate filling rate
	Call FillRate( Top,Rate )
	Write( *,* )
	Print*,'Fill rate is :', Rate


	!// Save result 
	FileName = 'Particle'
	Call SaveResult( FileName,N_Total )
End Subroutine ParticleFilling




!// Generate a particle
Subroutine GetPosition( Top,i,PosLabel )
	Implicit None
	Integer ::						i, j, PosLabel
	Real( kind=Myrealk ) ::			P, Lij, Top
	Character( len=512 ) ::			FileName


	!// Two methods of accumulate 
	Pt(i)%i = i
	Select Case( PosLabel )
	Case( 1 )
		Pt(i)%x = 0.5*(Br-Bl) + Bl
	Case( 2 )
		Call Random_Number(P)
		Pt(i)%x =  Bl + P*(Br-Bl)
	End Select
	Pt(i)%y = Top+Pt(i)%Radius


	!// Check if particle over the upper boder
	If( Pt(i)%Y-Pt(i)%Radius>Bt )then
		Print*
		Print*, 'Fill number is ',i-1
		FileName = 'Particle'
		Call SaveResult( FileName,i-1 )
		Stop
	End If

End Subroutine GetPosition




!// Move particle to the left until it contacts with boder or other particles
Subroutine Left( i,Njc )
	Implicit None
	Integer ::							i, j, Jc, Njc
	Real( kind=Myrealk ) ::				Lij, Xstep
	Integer ::							i1, i2, i3, l, M, N

	Xstep = 1.0E-3*Pt(i)%Radius
	Jc = 0
	Njc = 0
	Do while ( Jc<1 )
		Pt(i)%X = Pt(i)%X - Xstep
		If( Pt(i)%X-Pt(i)%Radius<Bl )then
			Pt(i)%X = Pt(i)%Radius + Bl
			Njc = N_Total + 1
			Exit
		End If
		M = Floor( (Pt(i)%Y-Bb)/MeshSize ) + 1
		N = Floor( (Pt(i)%X-Bl)/MeshSize ) + 1
		Do i1 = M-1,M
			Do i2 = N-1,N+1
				l = Mesh(0,i2,i1)
				If( l<1 ) Cycle
				Do i3 = 1, l
					j = Mesh(i3,i2,i1)
					If(j==i .or. j>N_Total ) cycle
					Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
					If( Pt(i)%Radius+Pt(j)%Radius>Lij )then
						Jc = Jc + 1
						Njc = j
						Exit 
					End If
				End Do 
			End Do 
		End DO 
	End Do 
	
End Subroutine Left




!// Move particle to the right until it contacts with boder or other particles
Subroutine Right( i,Njc )
	Implicit None
	Integer ::							i, j, Jc, Njc
	Real( kind=Myrealk ) ::				Lij, Xstep
	Integer ::							i1, i2, i3, l, M, N

	Xstep = 1.0E-3*Pt(i)%Radius
	Jc = 0
	Njc = 0
	Do while ( Jc<1 )
		Pt(i)%X = Pt(i)%X + Xstep
		If( Pt(i)%X+Pt(i)%Radius>Br )then
			Pt(i)%X = Br - Pt(i)%Radius
			Njc = N_Total + 1
			Exit
		End If
		M = Floor( (Pt(i)%Y-Bb)/MeshSize ) + 1
		N = Floor( (Pt(i)%X-Bl)/MeshSize ) + 1
		Do i1 = M-1,M
			Do i2 = N-1,N+1
				l = Mesh(0,i2,i1)
				If( l<1 ) Cycle
				Do i3 = 1, l
					j = Mesh(i3,i2,i1)
					If(j==i .or. j>N_Total ) cycle
					Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
					If( Pt(i)%Radius+Pt(j)%Radius>Lij )then
						Jc = Jc + 1
						Njc = j
						Exit 
					End If
				End Do 
			End Do 
		End DO 
	End Do 
	
End Subroutine Right




!// Move particle down until it contacts with boder or other particles
Subroutine Down( i,Njc )
	Implicit None
	Integer ::							i, j, Jc, Njc
	Integer ::							i1, i2, i3, l, M, N
	Real( kind=Myrealk ) ::				Lij, Ystep

	Ystep = 1.0E-3*Pt(i)%Radius
	Jc = 0
	Njc = 0
	Do while ( Jc<1 )
		Pt(i)%Y = Pt(i)%Y - Ystep
		If( Pt(i)%Y-Pt(i)%Radius<Bb )then
			Pt(i)%Y = Pt(i)%Radius + Bb
			Njc = N_Total + 2
			Exit
		End If
		M = Floor( (Pt(i)%Y-Bb)/MeshSize ) + 1
		N = Floor( (Pt(i)%X-Bl)/MeshSize ) + 1
		Do i1 = M-1,M
			Do i2 = N-1,N+1
				l = Mesh(0,i2,i1)
				If( l<1 ) Cycle
				Do i3 = 1, l
					j = Mesh(i3,i2,i1)
					If(j==i .or. j>N_Total) cycle
					Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
					If( Pt(i)%Radius+Pt(j)%Radius>Lij )then
						Jc = Jc + 1
						Njc = j
						Exit 
					End If
				End Do 
			End Do 
		End DO 
	End Do 
	
End Subroutine Down




!// Rotate particle around a contact one until contact with other particles
Subroutine Rotate( i,Njc )
	Implicit None
	Integer ::							i, j, k, Njc, Njcold
	Integer ::							i1, i2, i3, l, M, N
	Real( kind=Myrealk ) ::				Lij, Tstep, Thetaij, Lik, Xtemp

	l = 0
	Njcold = Njc


	If( Pt(i)%X>=Pt(Njc)%X )then
		Tstep = -0.001
		Xtemp = Bl
	Else
		Tstep = 0.001
		Xtemp = Br
	End If


	Lij = sqrt( (Pt(i)%X-Pt(Njc)%X)**2 + (Pt(i)%Y-Pt(Njc)%Y)**2 )
	Thetaij = Acos( (Pt(i)%X-Pt(Njc)%X)/Lij )
	If( Pt(i)%Y>Pt(Njc)%Y )then
		Thetaij = Thetaij
	Else
		Thetaij = 2.0*Pi - Thetaij 
	End If


	Do while( .True. )
		If( Tstep<0.0 )then
			If( Pt(i)%X<Xtemp .and. Abs(Pt(i)%Y-Pt(Njc)%Y)/Abs(Pt(i)%X-Pt(Njc)%X)<1.0 ) Exit
		Else
			If( Pt(i)%X>Xtemp .and. Abs(Pt(i)%Y-Pt(Njc)%Y)/Abs(Pt(i)%X-Pt(Njc)%X)<1.0 ) Exit
		End If
		If( Pt(Njcold)%Y+Pt(Njcold)%Radius<Pt(Njc)%Y+Pt(Njc)%Radius ) Exit
		Thetaij = Thetaij + Tstep
		Pt(i)%x = Pt(Njc)%X + Lij*Cos(Thetaij)
		Pt(i)%y = Pt(Njc)%Y + Lij*sin(Thetaij)
		If( Pt(i)%x-Pt(i)%Radius<Bl )then
			Pt(i)%x = Bl + Pt(i)%Radius
			Njc = N_Total + 1
			Exit
		End If
		If( Pt(i)%x+Pt(i)%Radius>Br )then
			Pt(i)%x = Br - Pt(i)%Radius
			Njc = N_Total + 3
			Exit
		End If
		If( Pt(i)%Y-Pt(i)%Radius<Bb )then
			Pt(i)%Y = Bb + Pt(i)%Radius
			Njc = N_Total + 2
			Exit
		End If
		M = Floor( (Pt(i)%Y-Bb)/MeshSize ) + 1
		N = Floor( (Pt(i)%X-Bl)/MeshSize ) + 1
		Loop1:Do i1 = M-1,M+1
			Do i2 = N-1,N+1
				l = Mesh(0,i2,i1)
				If( l<1 ) Cycle
				Do i3 = 1, l
					k = Mesh(i3,i2,i1)
					If(k==i .or. k==Njc .or. k>N_Total) cycle
					Lik = sqrt( (Pt(i)%X-Pt(k)%X)**2 + (Pt(i)%Y-Pt(k)%Y)**2 )
					If( Pt(i)%Radius+Pt(k)%Radius>Lik )then
						Njc = k
						Njcold = Njc
						Lij = sqrt( (Pt(i)%X-Pt(Njc)%X)**2 + (Pt(i)%Y-Pt(Njc)%Y)**2 )
						Thetaij = Acos( (Pt(i)%X-Pt(Njc)%X)/Lij )
						Xtemp = Pt(Njc)%X
						Exit Loop1
					End If
				End Do 
			End Do 
		End DO Loop1
	End Do 

End Subroutine Rotate




Subroutine Rotate1( i,Njc )
	Implicit None
	Integer ::							i, j, k, Njc
	Integer ::							i1, i2, i3, l, M, N
	Real( kind=Myrealk ) ::				Lij, Tstep, Thetaij, Lik, Xtemp

	l = 0
	Lij = sqrt( (Pt(i)%X-Pt(Njc)%X)**2 + (Pt(i)%Y-Pt(Njc)%Y)**2 )
	Thetaij = Acos( (Pt(i)%X-Pt(Njc)%X)/Lij )
	If( Pt(i)%Y>Pt(Njc)%Y )then
		Thetaij = Thetaij
	Else
		Thetaij = 2.0*Pi - Thetaij 
	End If

	If( Pt(i)%X>=Pt(Njc)%X )then
		Tstep = -0.001
		Xtemp = Bl
	Else
		Tstep = 0.001
		Xtemp = Br
	End If


	Do while( .True. )
		If( Tstep<0.0 )then
			If( Pt(i)%X<Xtemp ) Exit
		Else
			If( Pt(i)%X>Xtemp ) Exit
		End If
		Thetaij = Thetaij + Tstep
		Pt(i)%x = Pt(Njc)%X + Lij*Cos(Thetaij)
		Pt(i)%y = Pt(Njc)%Y + Lij*sin(Thetaij)
		If( Pt(i)%x-Pt(i)%Radius<Bl )then
			Pt(i)%x = Pt(i)%Radius + Bl
			Njc = N_Total + 1
			Exit
		End If
		If( Pt(i)%x+Pt(i)%Radius>Br )then
			Pt(i)%x = Br - Pt(i)%Radius
			Njc = N_Total + 3
			Exit
		End If
		If( Pt(i)%Y-Pt(i)%Radius<Bb )then
			Pt(i)%Y = Pt(i)%Radius + Bb
			Njc = N_Total + 2
			Exit
		End If
		M = Floor( (Pt(i)%Y-Bb)/MeshSize ) + 1
		N = Floor( (Pt(i)%X-Bl)/MeshSize ) + 1
		Loop1:Do i1 = M-1,M
			Do i2 = N-1,N+1
				l = Mesh(0,i2,i1)
				If( l<1 ) Cycle
				Do i3 = 1, l
					k = Mesh(i3,i2,i1)
					If(k==i .or. k==Njc .or. k>N_Total) cycle
					Lik = sqrt( (Pt(i)%X-Pt(k)%X)**2 + (Pt(i)%Y-Pt(k)%Y)**2 )
					If( Pt(i)%Radius+Pt(k)%Radius>Lik )then
						Njc = k
						j = k
						Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
						Thetaij = Acos( (Pt(i)%X-Pt(j)%X)/Lij )
						Xtemp = Pt(Njc)%X
						Exit Loop1
					End If
				End Do 
			End Do 
		End DO Loop1
	End Do 

End Subroutine Rotate1




!// Calculation of filling rate
Subroutine FillRate( Top,Rate )
	Implicit None
	Integer ::					j
	Real( kind=Myrealk ) :: 	V1, V2, Xmax=0.0, Ymax=0.0, Rate, Top

    Do j = 1, N_Total
		If( Pt(j)%Y<0.2*Top )then
			V1 = V1 + 3.14*Pt(j)%Radius * Pt(J)%Radius
		End If
    End Do
    V2 = 0.2*Top*(Br-Bl)
	Rate = V1/V2

End Subroutine FillRate




!// Configure file should like this:
!//--------------Not contains this line----------------
!// 1.0E-6				Accuary
!// 100					N_Total
!// 4					N_Sorts
!// 2.65E3				Rho
!// 0.0					Boder_lower
!// 10.0				Boder_upper
!// 0.0					Boder_left
!// 10.0				Boder_right
!// 1.0	2.0	3.0	4.0		Radius
!//	0.1	0.3	0.5	0.1		Ratio
!// 1					1.0
!// 2					3.0
!// ...					...
!// ...					...
!// ...					...
!// 100					4.0
!//--------------Not contains this line----------------

Subroutine SetVariable()
	Use Distribution
	Implicit None
	Integer ::								i, j, Err, Fi, Ri
	Real( kind=Myrealk ),Allocatable ::		Rad(:)
	Real( Kind=Myrealk ) ::					R, Rmin, Rmax, Mean, Sigma
	Character ::							C, Conf
	Character( len=100 ) ::					A(20)

	A(1) = '|  You can interactively enter the contents of the variable from'
	A(2) = '|  terminal or use a configuration file which is named as Configure'
	A(3) = '|  in the current directory. Please input Y/y to use a configure '
	A(4) = '|  file, or N/n to input from terminal'
	Call ShowInformation( A(1:4),4,3 )
	Do 
		Read( *,*,Iostat=Err ) Conf
		If( Err==0 ) Exit
	End Do 

	! Set global variables
	If( Conf=='Y' .or. Conf=='y' )then
		! Use a configuration file, which must be prepared in advance in a
		! certain format( Please refer to the template file
		! Open Configure file
		Open( 11,File='Configure',Status='Old', Iostat=Err )
		If( Err/=0 )then
			Print*,'Can not find Configure file'
			Read(*,*)
		End If


		! Read Accuracy
		Read( 11,*, Iostat=Err ) Accuracy
		If( Err /= 0 .or. Accuracy<=0.0 )then
			Print*,'Can not find Configure file'
			Read(*,*)
			Stop
		End If


		! Read N_Total
		Read( 11,* , Iostat=Err ) N_Total
		If( Err /= 0 .or.N_Total<=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Total is:', N_Total
			Read(*,*)
			Stop
		End IF
		Allocate( Pt(N_Total) )


		! Read N_Total
		Read( 11,*, Iostat=Err ) N_Sorts
		If( Err /= 0 .or. N_Sorts<=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Sorts is:', N_Sorts
			Read(*,*)
			Stop
		End IF


		! Read RHO
		Read( 11,*, Iostat=Err ) Rho
		If( Err /= 0 .or. Rho<=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'RHO is:', Rho
			Read(*,*)
			Stop
		End IF


		! Read Ym
		Read( 11,*, Iostat=Err ) Ym
		If( Err /= 0 .or. Ym<=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Ym is:', Ym
			Read(*,*)
			Stop
		End IF


		! Read Miv
		Read( 11,*, Iostat=Err ) Miv
		If( Err /= 0 .or. Miv<=0 .or. Miv>1.0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Miv is:', Miv
			Read(*,*)
			Stop
		End IF


		! Read lower border of the area
		Read( 11,*, Iostat=Err ) Bb
		If( Err /= 0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bb is:', Bb
			Read(*,*)
			Stop
		End IF


		! Read upper border of the area
		Read( 11,*, Iostat=Err ) Bt
		If( Err /= 0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bt is:', Bt
			Read(*,*)
			Stop
		Elseif( Bt<Bb )then
			Print*,'Bt should be larger than Bb'
			Print*,'Bb Bt:',Bb,Bt
			Read(*,*)
			Stop
		End IF



		! Read left border of the area
		Read( 11,*, Iostat=Err ) Bl
		If( Err /= 0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bb is:', Bl
			Read(*,*)
			Stop
		End IF


		! Read right border of the area
		Read( 11,*, Iostat=Err ) Br
		If( Err /= 0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bt is:', Br
			Read(*,*)
			Stop
		Elseif( Br<Bl )then
			Print*,'Bt should be larger than Bb'
			Print*,'Bb Bt:',Bl,Br
			Read(*,*)
			Stop
		End IF


		Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
		! Read radius of particle and ratio of erery radius
		If( N_Sorts/=N_total )then
			Ratio(:) = 0.0
			Radius(:) = 0.0
			Read( 11,* ) Radius(:)
			Read( 11,* ) Ratio(:)
		Else
			Read( 11,* ) c
			Read( 11,* ) C
		End If


		! Read particle radius
		Do i = 1, N_Total
			Read( 11,*,Iostat=Err ) Pt(i)%i,Pt(i)%Radius
			If( Err/=0 )then
				Print*,'Configure file maybe wrong'
				Print*,'N_Total of particle in the file is not enough'
				Read(*,*)
				Stop
			End If
			If( N_Sorts==N_Total ) Radius(i)=Pt(i)%Radius
		End Do 
		! Close configuration file
		Close( 11 )


	! Get variables from terminal 
	Else
		! Read Accuracy
		Print*
		Print*,'Input the Accuracy(<1.0E-3)'
		Do 
			Read( *,*,Iostat=Err ) Accuracy
			If( Err /= 0 .or. Accuracy>1.0E-3)then
				Print*,'Try again'
			Else
				Exit
			End If
		End Do 


		! Read N_Total
		Print*,'Input the total number of particles(Int)'
		Do 
			Read( *,*,Iostat=Err ) N_Total
			If( Err /= 0 )then
				Print*,'Try again'
			Else
				Exit
			End If
		End Do 
		! Allocate Pt
		Allocate( Pt(N_Total) )



		! Read left border of aera
		Print*,'Input the left border of area(Unit/M)'
		Do 
			Read( *,*,Iostat=Err ) Bl
			If( Err/=0 )then
				Print*,'Try again'
			Else
				Exit
			End If
		End Do 


		! Read right border of aera
		Print*,'Input the right border of area(Unit/M)'
		Do 
			Read( *,*,Iostat=Err ) Br
			If( Err/=0 )then
				Print*,'Try again'
			Elseif( Br<Bl )then
				Print*,'Br should be larger than Bl, Try again'
			Else
				Exit
			End If
		End Do 


		! Read lower border of aera
		Print*,'Input the lower border of area(Unit/M)'
		Do 
			Read( *,*,Iostat=Err ) Bb
			If( Err/=0 )then
				Print*,'Try again'
			Else
				Exit
			End If
		End Do 


		! Read upper border of aera
		Print*,'Input the upper border of area(Unit/M)'
		Do 
			Read( *,*,Iostat=Err ) Bt
			If( Err/=0 )then
				Print*,'Try again'
			Elseif( Bt<Bb )then
				Print*,'Bt should be larger than Bb, Try again'
			Else
				Exit
			End If
		End Do 


		! Get particle material properties
		! Read density of particle
		Print*,'Input the density of particle(Real)'
		Do 
			Read( *,*,Iostat=Err ) Rho
			If( Err/=0 .or. Rho<0.0 )then
				Print*,'Try again'
			Else
				Exit
			End If
		End Do 


		! Read Young's modulus
		Print*,'Please input Young modulus of particle'
		Do 
			Read( *,*,Iostat=Err ) Ym
			If( Err/=0 .or. Ym<0.0 )then
				Print*,'Ym: ',Ym
				Print*,'Try again'
			Else
				Exit
			End If
		End Do


		! Read Possion's Ratio 
		Print*,'Please input Possion ratio of particle'
		Do 
			Read( *,*,Iostat=Err ) Miv
			If( Err/=0 .or. Miv<0.0 .or. Miv>1.0 )then
				Print*,'Miv: ',Miv
				Print*,'Try again'
			Else
				Exit
			End If
		End Do


		! Set radius and ratio of particle
		Allocate ( Rad(N_Total) )
		Rad(:) = 0.0
		A(1)='|    Please choice the label of distribution of particle radius'
		A(2)='|  '
		A(3)='|        1>    Uniform Distribution radius'
		A(4)='|        2>    Random Distribution radius'
		A(5)='|        3>    Normal Distribution radius'
		!A(5)='// 4>    Log-Normal Distribution radius'
		Call ShowInformation( A(1:5),5,3 )
		Do
			Read( *,*,Iostat=Err ) Fi
			If( Err/=0 .or. Fi<1 .or. Fi>3 )then
				Print*,'Try again'
			Else
				Exit
			End IF
		End Do 
		If( Fi>1 )then
			A(1) = '|    Please choice the label of sorts number of particle radius'
			A(2) = '| '
			A(3) = '|        1>    Several specific radii'
			A(4) = '|        2>    Random size from given maximum to minimum radii or Mean and Sigma'
			Call ShowInformation( A(1:4),4,3 )
			Do
				Read( *,*,Iostat=Err ) Ri
				If( Err/=0 .or. Ri<1 .or. Ri>2 )then
					Print*,'Try again'
				Else
					Exit
				End IF
			End Do 
			Select Case( Ri )
			Case( 1 )
				! Read N_sorts
				print*,'Input the sorts of particles(Int) '
				Do 
					Read( *,*,Iostat=Err ) N_Sorts
					If( Err/=0 .or. N_Sorts<0 )then
						Print*,'Try again'
					Elseif( N_Sorts>N_Total )then
						Print*,'N_sorts > N_total Try again'
					Else
						Exit
					End IF
				End Do 
				Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
				Radius(:) = 0.0
				Ratio(:) = 0.0
				Print*,'Input the radius of particle(Unit/M)'
				Do i = 1, N_Sorts
					Print*,'Radius',i
					Do 
						Read( *,*,Iostat=Err ) Radius(i)
						If( Err/=0 .or. Radius(i)<0.0 )then
							Print*,'Try again'
						Else
							Exit
						End If
					End Do 
				End Do
				Ratio(:) = Real(N_Sorts)/1.0
				Select Case ( Fi )
				Case( 2 )
					Call RandomDistribution( Nt=N_Total,Ns=N_Sorts,Rs=Radius,R=Rad,Ratio=Ratio )
				Case( 3 )
					print*,'Input the Mean value of Normal distribution(Real)'
					Do 
						Read( *,*,Iostat=Err ) Mean
						If( Err/=0 )then
							Print*,'Try again'
						Else
							Exit
						End IF
					End Do 
					print*,'Input the variance(sigma) value of Normal distribution(Real)'
					Do 
						Read( *,*,Iostat=Err ) Sigma
						If( Err/=0 .or. Sigma<0.0 )then
							Print*,'Try again'
						Else
							Exit
						End IF
					End Do 
					Call NormalDistribution( Mean=Mean,Sigma=Sigma,Nt=N_Total,Ns=N_Sorts,Rs=Radius,R=Rad,Ratio=Ratio )
				End Select
				Do i = 1, N_Total
					Pt(i)%i = i
					Pt(i)%Radius = Rad(i)
				End Do 
			Case( 2 )
				N_sorts = N_total
				Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
				Select Case ( Fi )
				Case( 2 )
					print*,'Input the Minimum of particle radius(Raal)'
					Do 
						Read( *,*,Iostat=Err ) Rmin
						If( Err/=0 .or. Rmin<0.0 )then
							Print*,'Try again'
						Else
							Exit
						End IF
					End Do 
					print*,'Input the Maximum of particle radius(Raal)'
					Do 
						Read( *,*,Iostat=Err ) Rmax
						If( Err/=0 .or. Rmax<Rmin )then
							Print*,'Try again'
						Else
							Exit
						End IF
					End Do 
					Call RandomDistribution( Nt=N_Total,Ns=N_Sorts,Rmax=Rmax,Rmin=Rmin,R=Rad,Ratio=Ratio )
				Case( 3 )
					print*,'Input the Mean value of Normal distribution(Real)'
					Do 
						Read( *,*,Iostat=Err ) Mean
						If( Err/=0 )then
							Print*,'Try again'
						Else
							Exit
						End IF
					End Do 
					print*,'Input the variance(sigma) value of Normal distribution(Real)'
					Do 
						Read( *,*,Iostat=Err ) Sigma
						If( Err/=0 .or. Sigma<0.0 )then
							Print*,'Try again'
						Else
							Exit
						End IF
					End Do 
					Call NormalDistribution( Mean=Mean,Sigma=Sigma,Nt=N_Total,Ns=N_Sorts,R=Rad,Ratio=Ratio )
				End Select
				Do i = 1, N_Total
					Radius(i) = Rad(i)
					Pt(i)%i = i
					Pt(i)%Radius = Rad(i)
				End Do 
			End Select
		Else
			N_sorts = 1
			Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
			print*,'Input the radius of uniform radius(real) '
			Do 
				Read( *,*,Iostat=Err ) Radius(1)
				If( Err/=0 .or. Radius(1)<0.0 )then
					Print*,'Try again'
				Else
					Exit
				End IF
			End Do 
			Do i = 1, N_Total
				Pt(i)%i = i
				Pt(i)%Radius = Radius(1)
			End Do 
			Ratio(1) = 1.0
		End IF


	End If 

	! Maximum and minimum value of radius
	MaxRadius = Maxval( Radius )
	MinRadius = Minval( Radius )

End Subroutine SetVariable

End Module Accumulation




!// This module contains subroutines calculate particle distribution
Module Distribution
	Implicit None
	Integer,Parameter ::								Myrealk=4, Myintk=4

Contains

!// Uniform distribution
Subroutine UniformDistribution( Radii,Nt,R )
	Implicit None
	Integer ::												i
	Real( kind=Myrealk ), Intent( In ) ::					Radii
	Integer ::												Nt !// Total number of particle
	Integer ::												Ns !// Total number of sorts
	Real( kind=MyrealK ), Allocatable ::					R(:) !// Results

	Allocate( R(Nt) )
	R(:) = Radii

End Subroutine UnIformDistribution




!// Random distribution
Subroutine RandomDistribution( Nt,Ns,Rs,Rmax,Rmin,R,Ratio )
	Implicit None
	Integer ::														i, j
	Integer, Intent( In ) ::										Nt !// Total number of particle
	Integer, Intent( In ) ::										Ns !// Total number of sorts
	Real( kind=MyrealK ), Intent( Inout ), Allocatable ::			R(:) !// Results
	Real( kind=Myrealk ), Intent( In ), Optional ::					Rmax, Rmin !// Maximum and minimum of radius
	Real( kind=Myrealk ), Intent( In ), Allocatable, Optional ::	Rs(:) !// SpecIfic radius
    Real( kind=Myrealk ), Intent( Out), Optional, Allocatable ::	Ratio(:) !// Particle Ratio
	Real( kind=MyrealK ), Allocatable ::							P(:) !// Random number

	
	Allocate( P(Nt) )
	If( .Not. Allocated(R) )then
		Allocate( R(Nt) )
	End If


	Call RANDOM_SEED()
	Call Random_Number( P )


	If( Ns==Nt )then
		!// All particles have dIfferent radius, random size from Rmin to Rmax
		If( .Not. Present(Rmax) .or. .Not. Present(Rmin) )then
			Print*,'Please input the maximum and minimum of particle radius'
			Stop
		End If
		Do i = 1, Nt
			R( i ) = Rmin + (Rmax-Rmin)*P(i)
		End Do 
		If( Present(Ratio) )then
			If( .Not. Allocated(Ratio) )then
				Allocate( Ratio(Ns) )
			End If
			Ratio(:) = 1.0/Real(Nt)
		End If
	Else
		!// All particle have specIfic radius set in the Rs array
		If( .Not. Present(Rs) )then
			Print*,'Argument is not enough--RandomDistribution'
			stop
		End If
		Do i = 1, Nt
			j = Ceiling( P(i)/(1.01/Ns) )
			R(i) = Rs(j)
		End Do
		If( Present(Ratio) )then
			If( .Not. Allocated(Ratio) )then
				Allocate( Ratio(Ns) )
			End If
			Ratio(:) = 1.0/Real(Ns)
		End If
	End If

End Subroutine RandomDistribution




!// Normal distribution
Subroutine NormalDistribution( Mean,Sigma,Nt,Ns,Rs,R,Ratio )
	Implicit None
	Integer ::														i, j, Flag=0
	Integer ::														Nt !// Total number of particle
	Integer ::														Ns !// Total number of sorts
	Real( kind=MyrealK ), Allocatable ::							R(:) !// Results
    Real( kind=Myrealk ), Parameter ::								Pi = 3.14159
    Real( kind=Myrealk ) ::											U(2), y1, y2 
	!Real( kind=Myrealk ), Intent( In ), Optional ::					Rmax, Rmin !// Maximum and minimum of radius
	Real( kind=Myrealk ), Intent( In ), Allocatable, Optional ::	Rs(:) !// SpecIfic radius
	Real( Kind=Myrealk ), Intent( In ) ::							Mean !// Mean value of Normal-distribution
	Real( Kind=Myrealk ), Intent( In ) ::							Sigma !// Standard Sigma of Normal-distribution
    Real( kind=Myrealk ), Intent( Out), Optional, Allocatable ::	Ratio(:) !// Particle Ratio

	If( Present(Ratio) )then
		If( .Not. Allocated(Ratio) ) Allocate( Ratio(Ns) )
	End If
	If( .Not. Allocated(R) ) Allocate( R(Nt) )
	Call RANDOM_SEED()
	If( Nt==Ns )then
		!// All particles have dIfferent radius, random size from Rmin to Rmax
		i = 0
		Do 
			i = i + 1
			If( I>Nt ) Exit
			Call Random_Number( U )
			If( Flag==0 ) then 
				y1 = sqrt( -2.0*log(U(1)) )*cos( 2.0*pi*U(2))  
				R(i) = Mean + Sigma*y1 
				If( R(i)<0.0 ) i = i - 1
				Flag = 1 
			else 
				y2 = sqrt( -2.0*log(U(1)) )*sin( 2.0*pi*U(2) ) 
				R(i) = Mean + Sigma*y2 
				If( R(i)<0.0 ) i = i - 1
				Flag = 0 
			EndIf  
		End Do 
	Else
		!// All particle have specIfic radius set in the Rs array
		Call Normal_Ratio( Mean,Sigma,Ns,Ratio )
		Call Normal_Radius( Nt,Ns,Rs,Ratio,R )
	End If

End Subroutine NormalDistribution




Subroutine Normal_Ratio( Mean,Sigma,N,Ratio )
    Implicit None
    Integer ::												i
    Integer, Intent( In ) ::								N
    Real( kind=Myrealk ) ::									Start, Finish, JianGe, S, Step
    Real( kind=Myrealk ), Intent( In ) ::					Sigma, Mean
    Real( kind=Myrealk ), Allocatable ::					Ratio(:)
    
    Start = -4.0*Sigma + Mean
    Finish = 4.0*Sigma + Mean
    JianGe = (Finish-Start)/Float(N)
    S = Start
    i = 1
    Step = 1.0E-5
    Ratio(:) = 0.0
    Do 
        Ratio(i) = Ratio(i) + Step*1.0/Sqrt(2.0*3.1415*Sigma**2)*Exp(-(S-Mean)**2/2.0/Sigma/Sigma)
        S = S + Step
        If(S>JianGe*Float(i)+Start)then
            i = i + 1
            If( i>N )then
                Exit
            End If
        End If
    End Do 
    Do i = 1, N
        Ratio(i) = Ratio(i)/Sum(Ratio)
    End Do 

End Subroutine Normal_Ratio




Subroutine Normal_Radius( Nt,Ns,Rs,Ratio,RadiusList )
    Implicit None
    Integer ::										j, i, Nt, Ns
	Integer, Allocatable ::							N(:)
    Real( kind=MyrealK ), Allocatable ::			RadiusList(:), Sj
    Real( kind=Myrealk ), Allocatable ::			Ratio(:), Rs(:)
	Character ::									a*3, aa*6

	Allocate( N(Ns) )
	N(:) = 0
    RadiusList(:) = 0.0
    Call RANDOM_SEED
    Do i = 1, Ns    
        N(i) = Int(Real(Nt)*Ratio(i))    
        If( N(i)==0 )then    
            N(i) = 1    
        End IF    
    End Do     
    N(Int(float(Ns)*0.5)) = N(Int(float(Ns)*0.5)) + Nt - Sum(N)    
	Do i = 1, Ns	
		Write(a,'(g0)') i	
		Write(aa,*) "R0",Trim(a),"="	
		Write(*,"(a,i8)") aa, N(i)	
	End Do	
    If ( N(Int(float(Ns)*0.5)) < 0 )then    
        Print*,'Particle number wrong'
        Read(*,*)    
    End If    
    j = 1    
    Do i = 1, Ns    
        RadiusList(j:N(i)+j-1) = Rs(i)    
        j = j + N(i)    
    End Do    
    Call UpsetList(RadiusList)          


ENd Subroutine Normal_Radius




Subroutine UpsetList( A )
    Implicit None
    Real( kind=Myrealk ) ::		Temp, R
    Real( kind=Myrealk ) ::		A(:)
    Integer ::					i, j, Narray
    
    Narray = Size(A)
    Call RANDOM_SEED()
    Do i = Narray, 2, -1
        Call Random_Number(r)   
        j = Int(r*Narray)+1
        Temp = A(j)
        A(j) = A(i)
        A(i) = Temp
    End Do 

End Subroutine UpsetList 


End Module Distribution




!// This module contains all functions of file operate and screen output 
Module FileOperate
	Use GlobalVariable
	Implicit None

Contains 

!// Description of the project  
Subroutine Description()
	Implicit None
	Integer ::								i
	Character(len=100) ::					A, B, C(20)

	A = '|-------------------------------------------------------------------------------------|'
	C(1) = '| '
	C(2) = '|  Particle parking simulation. '
	C(3) = '| '
	C(4) = '|  Tightly stacked balls in a 2D rectangular. Similar to Y.T.Feng(2003)  '
	C(5) = '|  but they do not consider structural balance and compression of the '
	C(6) = '|  particles. This algorithm use a new method which simulated gravity '
	C(7) = '|  of the structure '
	C(8) = '| '
	C(9) = '|  Author:	MU LiAn LanZhou University '
	C(10) = '|  Date:	2017/07/06 '
	C(11) = '|  Mail:	1791430496@qq.com '
	C(12) = '|  Edition:	1.0 '
	Write(*,*)
	Write(*,*)
	Write(*,"(A100)") Adjustl(A)
	Do i = 2, 12
		B = C(i)
		Write(*,"(A100)") Adjustl(B)
	End Do 
	Write(*,"(A100)") Adjustl(A)
	Write(*,*)
	Write(*,*)

End Subroutine Description




!// Terminal output
Subroutine ShowInformation( C,N,FormatLabel )
	Implicit None
	Integer ::								i, N, FormatLabel
	Character(len=100) ::					A, B
	Character(len=100) ::					C(:)

	A = '|-------------------------------------------------------------------------------------|'
	Select Case( FormatLabel )
	Case( 0 )
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
	Case( 1 )
		Write(*,"(A100)") Adjustl(A)
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
	Case( 2 )
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
		Write(*,"(A100)") Adjustl(A)
	Case( 3 )
		Write(*,"(A100)") Adjustl(A)
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
		Write(*,"(A100)") Adjustl(A)
	End Select 

End Subroutine ShowInformation



	
!// Write result to specific file named FileName
Subroutine SaveResult( FileName,N )
	Implicit None
	Integer ::						i, N
	Character ::					FileName*512

	FileName = Trim(Adjustl(FileName))
	Open( 99,File=FileName )
	Do i = 1, N
		Write( 99,'(i13,2x,3(E13.6,2x))' )  i,Pt(i)%X,Pt(i)%Y,Pt(i)%Radius
	End Do 
	Close( 99 )

End Subroutine SaveResult




!// Close specific file
Subroutine CloseFile( Filepointer )
	Integer ::			Filepointer

	Close( Filepointer )

End Subroutine




!// Delete specific file
Subroutine DeleteFile( Filepointer )
	Integer ::			Filepointer

	Close( Filepointer,Status='Delete' )

End Subroutine DeleteFile




!// Write program parameter into configuration file
Subroutine SaveConfigure()
	Implicit None
	Integer ::				j
	Character ::						WFormat*30

	Open( 11, File='Configure' )
	Write( 11,"(E10.3,2x,A30)" ) Accuracy, 'Newton ineration Accuracy'
	Write( 11,"(i10,2x,A30)" ) N_Total, 'Particle N_Total'
	Write( 11,"(i10,2x,A30)" ) N_Sorts, 'Particle sorts'
	Write( 11,'(E10.3,2x,A30)' ) RHO,'Density of particle'
	Write( 11,'(E10.3,2x,A30)' ) Ym,'Young Modulus'
	Write( 11,'(E10.3,2x,A30)' ) Miv,'Poisson ratio'
	Write( 11,"(E10.3,2x,A30)" ) Bb,'Bottom edge'
	Write( 11,'(E10.3,2x,A30)' ) Bt,'Top edge'
	Write( 11,'(E10.3,2x,A30)' ) Bl,'Left edge'
	Write( 11,'(E10.3,2x,A30)' ) Br,'Right edge'
	If( N_Sorts==N_Total )then
		Write( 11,* ) 'Too many Radius'
		Write( 11,* ) 'To many Ratio'
	Else
		Write( WFormat,* ) N_Sorts
		WFormat = '('//Trim(Adjustl(WFormat))//'(E10.3),2x,A30)'
		Write( 11,WFormat ) Radius,'Radius'
		Write( 11,WFormat ) Ratio, 'Ratio'
	End If
    Do j = 1, N_Total
		Write(11, "(I6, 2(f13.6))") Pt(j)%i, Pt(j)%Radius
    End Do
	Close(11)

End Subroutine SaveConfigure


End Module FileOperate




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




Program Main
	Use Accumulation
	Use Sedimentation
	Use FileOperate
	Implicit None
	Character( len=100 ) ::				A(100)
	Integer ::							i, Err

	Call Description()

	A(1) = '|    Please select the label of the following module'
	A(2) = '|'
	A(3) = '|        1>    Particle accumulation'
	A(4) = '|        2>    Particle sedimentation'
	A(5) = '|        3>    Particle accumulation and sedimentation'
	Call ShowInformation( A(1:5),5,3 )
	Do
		Read( *,*,Iostat=Err ) i
		If( Err/=0 .or. i<0 .or. i>3 )then
			Print*,'Try again'
		Else
			Exit
		End If
	End Do 


	Select Case( i )
	Case( 1 )
		Call ParticleFilling()
	Case( 2 )
		Call ParticleSettlement()
	Case( 3 )
		Call ParticleFilling()
		Call ParticleSettlement()
	End Select 



End Program Main

Module MeshGrid
	Implicit None
	Integer :: 							N_ColumnOfMesh ! Column of mesh 
	Integer :: 							N_RowOfMesh ! Row of mesh
	Integer,Allocatable ::				Mesh(:,:,:) ! Mesh information 
	Real( kind=4 ) ::					MeshSize ! Mesh size

Contains


!// Mesh the two-dimensional area by the given size
Subroutine Get_Mesh_2D( N,Boder_Left,Boder_Right,Boder_lower,Boder_upper,Rmin,MeshSize )
	Implicit None
	Integer ::						i, j, N
	Real( kind=4 ),Parameter ::		PI=3.14159
	Real( kind=4 ) ::				Rmin
	Real( kind=4 ) ::				MeshSize ! Mesh size
	Real( kind=4 ) :: 				Boder_Left ! The left border of the area
	Real( kind=4 ) :: 				Boder_Right ! The right border of the area
	Real( kind=4 ) :: 				Boder_lower ! The lower border of the area
	Real( kind=4 ) :: 				Boder_upper ! The upper border of the area

	!// Row and column number of mesh grid
	N_ColumnOfMesh = Int( (Boder_Right-Boder_Left)/MeshSize )
	MeshSize = (Boder_Right-Boder_Left)/Real(N_ColumnOfMesh)
	N_RowOfMesh = Int( (Boder_Upper-Boder_lower)/MeshSize )
	i = Int( MeshSize**2/Pi/Rmin**2 ) + 1
	Allocate( Mesh(0:i,0:N_ColumnOfMesh+1,0:N_RowOfMesh+1) )

	!// Initialize the mesh array 
	!// Set the out grid of area to a specific particle
	Mesh(:,:,:) = 0
	!// Left grid 
	Mesh(1,0,:) = N+1
	Mesh(0,0,:) = 1
	!// Right grid 
	Mesh(1,N_ColumnOfMesh+1,:) = N+3
	Mesh(0,N_ColumnOfMesh+1,:) = 1
	!// Lower grid 
	Mesh(1,:,0) = N+2
	Mesh(0,:,0) = 1
	!// Upper grid 
	Mesh(1,:,N_RowOfMesh+1) = N+4
	Mesh(0,:,N_RowOfMesh+1) = 1

End Subroutine Get_Mesh_2D




!// Update particle mesh information
Subroutine MeshUpdate( i,X,Y,M,N,Boder_Left,Boder_lower )
	Implicit None
	Integer ::						i, M, N, j
	Real( kind=4 ) ::				X, Y
	Real( kind=4 ) :: 				Boder_Left ! The left border of the area
	Real( kind=4 ) :: 				Boder_lower ! The lower border of the area
	
	M = Floor( (Y-Boder_lower)/MeshSize ) + 1
	N = Floor( (X-Boder_Left)/MeshSize ) + 1
	If( M>N_RowOfMesh )then
		Print*,'MeshGrid:',i,M
		stop
	End If
	If( N>N_ColumnOfMesh )then
		Print*,'MeshGrid:',i,N
		stop
	End If
	j = Mesh(0,n,m)
	Mesh(j+1,n,m) = i
	Mesh(0,n,m) = Mesh(0,n,m) + 1

End Subroutine MeshUpdate


End Module MeshGrid




Module ScheduleBar
	Implicit None
	private
	Logical, parameter, public ::	CMD_PROGRESS_ABSOLUTE = .True.
	Type, public ::					CLS_CMD_Progress
		Integer , private ::		N, lens, i
		Character ::				M = "*", O = "."
		Character(len=64) ::		Prefix
	Contains
		Procedure ::				Set
		Procedure ::				Put
	End Type CLS_CMD_Progress

contains

Subroutine Set( this, N, L )
	Class( CLS_CMD_Progress ) ::	this
	Integer , Intent( IN )	::		N , L
	this % N    = N
	this % lens = L
	this % i = 0
	this % Prefix = " Progress: " !//
End Subroutine Set




Subroutine Put( this , k, bAbsol )
	Class( CLS_CMD_Progress ) ::		this
	Integer , Intent( IN ) ::			k
	Logical , optional ::				bAbsol
	Character( len=1 ) ::				br
	Integer ::							jm

	this%i = this%i + k 
	If( present(bAbsol) )then
		if( bAbsol ) this%i = k
	End If
	If( this%i>this%n ) this%i = this%n
	jm = Nint( Real(this%i*this%lens)/Real(this%n) )
	If( jm<1 ) jm=1
	If( this%i<this%n )then
		br = char(13)
	Else
		br = char(10)
	End If
	Write( *,"(5A,F6.2,2A\)" ) Trim(this%Prefix),'[', &
		Repeat( this%m,jm ), Repeat( this%o,this%lens-jm ), ']', &
		this%i*100.0/this%n, "%", br
End Subroutine Put


End Module ScheduleBar

!//--------------------------------------------------------------------//
!// This module is the main algorithm of particle sedimentation
!//		To generate a particle sedimentation you should call 
!//		ParticleSettlement()
!//--------------------------------------------------------------------//

Module Sedimentation
	Use GlobalVariable
	Use FileOperate
	Implicit None

	Real( kind=Myrealk ) :: 				SmallValue ! A small value for control iteration
	Integer(kind=Myintk),Allocatable ::		Neighbor(:,:) ! Neighbor of particles 
	Real( kind=Myrealk ), Allocatable :: 	X(:), Y(:) ! The value of Newton iteration of particle location
	Real( kind=Myrealk ), Allocatable :: 	Xold(:), Yold(:) ! The value of Newton iteration of particle location
	Real( kind=Myrealk ), Allocatable :: 	Ux(:), Uy(:), U2x(:), U2y(:)

Contains


Subroutine ParticleSettlement()
	Implicit None
	Integer ::								i, j, l, N
	Logical ::								Convergence
	Character ::							FileName*512
	Real( kind=Myrealk ) ::					Xtemp, Ytemp, E1, E2, E
	Real( kind=Myrealk ) ::					Rij, Kij, Lij, S1, S2, DeltaLij
	Real( kind=Myrealk ) :: 				Dldx2, Dldy2, Dldx, Dldy
	
	
	! Set global variable for sedimentation
	Call SetVariable()


	! Get initial particle neighbors and write to file
	FileName = 'NeighborOri'
	Call GetNeighbor( FileName )


	! Set initial value of iteration
	Call SetInitialValueOfIteration()
	

	! Loop for iteration until particle coordinates almost no longer change
	Do while( .True. ) 
		Xold = X
		Yold = Y
		Convergence = .True.
		Do i = 1, N_Total !Loop i
			E1 = Ym*Pt(i)%Radius
			Do l = 1, Neighbor(i,1)
				j = Neighbor(i,l+1)
				If( J==N_Total+1 .Or. J==N_Total+3 )then
					E2 = 1.0E12
					E = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
					Rij = Pt(i)%Radius
					Kij = 1.34*E*Rij**0.5 
					Ux(i) = Ux(i) - kij*Xold(i)
					U2x(i) = U2x(i) - Kij
					Cycle
				End IF	
				If( J==N_Total+2 )then
					E2 = 1.0E12
					E = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
					Rij = Pt(i)%Radius
					Kij = 1.34*E*Rij**0.5
					Uy(i) = Uy(i) + kij*Yold(i)
					U2y(i) = U2y(i) + Kij
					Cycle
				End IF	
				Lij = ( (Pt(i)%X-Pt(j)%x)**2+(Pt(i)%Y-Pt(j)%Y)**2 )**0.5
				E2 = Ym*Pt(j)%Radius
				E = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
				Rij = Pt(i)%Radius*Pt(j)%Radius /(Pt(i)%Radius+Pt(j)%Radius)
				Kij = 1.34*E*Rij**0.5
				S1 = (Xold(i)-Xold(j))**2 + 2.0*(Pt(i)%X-Pt(j)%X)*(Xold(i)-Xold(j))
				S2 = (Yold(i)-Yold(j))**2 + 2.0*(Pt(i)%Y-Pt(j)%Y)*(Yold(i)-Yold(j))
				DeltaLij = (Lij**2+S1+S2)**0.5 - Lij 
				Dldx = ( (Xold(i)-Xold(j)) + (Pt(i)%X-Pt(j)%X) )/(DeltaLij+Lij)
				Dldy = ( (Yold(i)-Yold(j)) + (Pt(i)%Y-Pt(j)%Y) )/(DeltaLij+Lij)
				Dldx2 = (DeltaLij+Lij-(Pt(i)%X-Pt(j)%X)*(Xold(i)-Xold(j))*Dldx)/(DeltaLij+Lij)**2
				Dldy2 = (DeltaLij+Lij-(Pt(i)%Y-Pt(j)%Y)*(Yold(i)-Yold(j))*Dldy)/(DeltaLij+Lij)**2
				Ux(i) = Ux(i) + kij*DeltaLij*Dldx
				Uy(i) = Uy(i) + Kij*DeltaLij*Dldy
				U2x(i) = U2x(i) + Kij*( DeltaLij*Dldx2+Dldx**2 )
				U2y(i) = U2y(i) + Kij*( DeltaLij*Dldy2+Dldy**2 )
			End Do 
			Uy(i) = Uy(i) - Pt(i)%Mass*G
			Xtemp = X(i) 
			Ytemp = Y(i)
			If( U2x(i)/=0.0 )then
				X(i) = X(i) + Ux(i)/U2x(i)
			Else
				X(i) = X(i)
			End If
			If( U2y(i)/=0.0 )then
				Y(i) = Y(i) + Uy(i)/U2y(i)
			Else
				Y(i) = Y(i)
			End If
			Write( *,'(A40,2(F13.6,2x),2x,i9,1a)',Advance='No' ) 'Iteration of X and Y:', X(i),Y(i), i, char(13)
			If( Abs(X(i)-Xtemp)>SmallValue .Or. Abs(Y(i)-Ytemp)>SmallValue )then
				Convergence = .False.
			End If
		End Do !Loop i
		If( Convergence )then
			Exit
		End If
	End Do 
	Write( *,'(A40)' ) 'Sedimentation has done'


	! Save iteration value and get new neighbor of particle 
	Do i = 1, N_Total
		Pt(i)%X = Pt(i)%X + X(i)
		Pt(i)%Y = Pt(i)%Y + Y(i)
	End Do 
	FileName = 'NeighborNew'
	Call GetNeighbor(FileName)


	! Write particle information into file
	FileName = 'Sedimentation'
	Call SaveResult( FileName,N_Total )


End Subroutine ParticleSettlement




Subroutine SetVariable()
	Implicit None
	Integer ::			i, Err
	Character ::		C

	
	If( .Not. Allocated( Pt ) )then
		! Open Configure file
		Open( 11,File='Configure',Status='Old', Iostat=Err )
		If( Err/=0 )then
			Print*,'Can not find Configure file'
			Read(*,*)
		End If


		! Read Accuracy
		Read( 11,*, Iostat = Err ) Accuracy
		If( Err /= 0 )then
			Print*,'Can not find Configure file'
			Read(*,*)
			Stop
		End If


		! Read N_Total
		Read( 11,* , Iostat = Err ) N_Total
		If( Err/=0 .or. N_Total<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Total is:', N_Total
			Read(*,*)
			Stop
		End IF
		Allocate( Pt(N_Total) )


		! Read N_Sorts
		Read( 11,*, Iostat = Err ) N_Sorts
		If( Err/=0 .or. N_Sorts<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Sorts is:', N_Sorts
			Read(*,*)
			Stop
		End IF


		! Read RHO
		Read( 11,*, Iostat = Err ) Rho
		If( Err/=0 .or. RHO<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'RHO is:', Rho
			Read(*,*)
			Stop
		End IF


		! Get particle material properties
		! Read Young's modulus
		Read( 11,*, Iostat = Err ) Ym
		If( Err/=0 .or. Ym<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Yong modulus is:', Ym
			Read(*,*)
			Stop
		End IF


		! Read Possion's Ratio 
		Read( 11,*, Iostat = Err ) Miv
		If( Err/=0 .or. Miv<0.0 .or. Miv>1.0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Poisson ratio is:', Miv
			Read(*,*)
			Stop
		End IF


		! Read lower border of the area
		Read( 11,*, Iostat = Err ) Bb
		If( Err/=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bb is:', Bb
			Read(*,*)
			Stop
		End IF


		! Read upper border of the area
		Read( 11,*, Iostat = Err ) Bt
		If( Err/=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bt is:', Bt
			Read(*,*)
			Stop
		Elseif( Bt<Bb )then
			Print*,'Bt should be larger than Bb'
			Print*,'Bb Bt:',Bb,Bt
			Read(*,*)
			Stop
		End IF



		! Read left border of the area
		Read( 11,*, Iostat = Err ) Bl
		If( Err/=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bb is:', Bl
			Read(*,*)
			Stop
		End IF


		! Read right border of the area
		Read( 11,*, Iostat = Err ) Br
		If( Err/=0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Bt is:', Br
			Read(*,*)
			Stop
		Elseif( Br<Bl )then
			Print*,'Bt should be larger than Bb'
			Print*,'Bb Bt:',Bl,Br
			Read(*,*)
			Stop
		End IF


		Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
		! Read radius of particle and ratio of erery radius
		If( N_Sorts/=N_total )then
			Ratio(:) = 0.0
			Radius(:) = 0.0
			Read( 11,* ) Radius(:)
			Read( 11,* ) Ratio(:)
		Else
			Read( 11,* ) C
			Read( 11,* ) C
		End If


		! Close configuration file
		Close( 11 )


		! Initialize array
		Allocate( X(N_Total), Y(N_Total) )
		Allocate( Xold(N_Total), Yold(N_Total) )
		Allocate( Ux(N_Total), U2x(N_Total), Uy(N_Total), U2y(N_Total) )
		X = 0.0
		Y = 0.0
		Xold = 0.0
		Yold = 0.0
		Ux = 0.0
		U2x = 0.0
		Uy = 0.0
		U2y = 0.0


		! Read particle information from file
		Open( 11,File='Particle',Status='Old', Iostat=Err )
		If( Err/=0 )then
			Print*,'Open file Particle wrong'
			stop
		End If
		Rewind( 11 )
		Do i = 1, N_Total
			Read( 11,* ) Pt(i)%i,Pt(i)%x,Pt(i)%Y,Pt(i)%Radius
			Pt(i)%Mass = Rho*4.0*Pi*Pt(i)%Radius**3/3.0
		End Do 
		Close( 11 )

	Else

		! Initialize array
		Allocate( X(N_Total), Y(N_Total) )
		Allocate( Xold(N_Total), Yold(N_Total) )
		Allocate( Ux(N_Total), U2x(N_Total), Uy(N_Total), U2y(N_Total) )
		X = 0.0
		Y = 0.0
		Xold = 0.0
		Yold = 0.0
		Ux = 0.0
		U2x = 0.0
		Uy = 0.0
		U2y = 0.0

	End If


	MinRadius = Minval( Radius )
	MaxRadius = Maxval( Radius )
	SmallValue = MinRadius*1.0E-1
	i = Int( Pi*(MaxRadius+MinRadius)/MinRadius ) + 4
	Allocate( Neighbor(N_Total,i)  )

End Subroutine SetVariable




Subroutine GetNeighbor( FileName )
	Implicit None
	Integer :: 					i, j, k, N, Col, Row, l, Njc
	Real( kind=Myrealk ) ::		Lij, Cd
	Character( len=512 ) :: 	Ch, FileName


	FileName = Trim(Adjustl(FileName))
	Open( 111,File=FileName )
	Neighbor = 0
	Njc = 0
	Do i = 1, N_Total
		Do j = 1, N_Total
			If( i==j )then
				Cycle
			End If
			Lij = Sqrt( (Pt(i)%x-Pt(j)%x)**2 + (Pt(i)%y-Pt(j)%y)**2 )
			If( Lij < Pt(i)%Radius+Pt(j)%Radius )then
				n = Neighbor(i,1)
				Neighbor(i,n+2) = j
				Neighbor(i,1) = Neighbor(i,1)+1
			ElseIf( Abs(Lij-Pt(i)%Radius-Pt(j)%Radius)<SmallValue )then
				n = Neighbor(i,1)
				Neighbor(i,n+2) = j
				Neighbor(i,1) = Neighbor(i,1)+1
			End If
		End Do
		If( Abs(Pt(i)%x-Pt(i)%Radius-Bl)<SmallValue )then
			N = Neighbor(i,1)
			Neighbor(i,N+2) = N_Total+1
			Neighbor(i,1) = Neighbor(i,1) + 1
		End If
		If( Abs(Pt(i)%Y-Pt(i)%Radius-Bl)<SmallValue )then
			N = Neighbor(i,1)
			Neighbor(i,N+2) = N_Total+2
			Neighbor(i,1) = Neighbor(i,1) + 1
		End If
		If( Abs(Pt(i)%x+Pt(i)%Radius-Br)<SmallValue )then
			N = Neighbor(i,1)
			Neighbor(i,N+2) = N_Total+3
			Neighbor(i,1) = Neighbor(i,1) + 1
		End If
		N = Neighbor(i,1)
		if( N==0 )then
			Print*, i,'neighbor is 0'
		End If
		Write( Ch,* ) N+2
		Ch = '('//Trim(Adjustl(Ch))//'(i9))'
		Write( 111,ch ) i,N,Neighbor(i,2:N+1)
		Ch = ''
	End Do 
	Close( 111 )

End Subroutine GetNeighbor




Subroutine SetInitialValueOfIteration()
	Implicit None
	Integer ::									i
	Real( kind=MyrealK ), Allocatable :: 		Sj(:)
	Real( kind=MyrealK ) ::						Initial
	
	Allocate( Sj(N_Total) )
	Call Random_seed()
	Call Random_Number( Sj )
	Initial = MinRadius * 1.0E-3
	Do i = 1, N_Total
		X(i) = (Sj(i)-0.5)/Abs(Sj(i)-0.5) * Initial*Sj(i)
		!X(i) = -1.0 * Initial*Sj(i)
		Y(i) = -1.0 * Initial*Sj(i)
	End Do 
	Deallocate( Sj )

End Subroutine SetInitialValueOfIteration



End Module Sedimentation
