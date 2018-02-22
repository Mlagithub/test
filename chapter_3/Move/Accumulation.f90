
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


	Call SetVariable()

	Call Progress%Set( N=N_Total,L=30 )
	Progress%Prefix = "Accumulate Particle: "
	Progress%m = "#"
	Progress%o = "."

	Call SaveConfigure()
	MeshSize = 2.0*MaxRadius
	Call Get_Mesh_2D( N_Total,Bl,Br,Bb,Bt,MinRadius,MeshSize )

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

	Call FillRate( Top,Rate )
	FileName = 'Particle'
	Call SaveResult( FileName,N_Total )

	Write( *,* )
	Print*,'Fill rate is :', Rate

End Subroutine ParticleFilling


!// Generate a particle
Subroutine GetPosition( Top,i,PosLabel )
	Implicit None
	Integer ::						i, j, PosLabel
	Real( kind=Myrealk ) ::			P, Lij, Top
	Character( len=512 ) ::			FileName

	Pt(i)%i = i
	Select Case( PosLabel )
	Case( 1 )
		Pt(i)%x = 0.5*(Br-Bl) + Bl
	Case( 2 )
		Call Random_Number(P)
		Pt(i)%x =  Bl + P*(Br-Bl)
	End Select
	Pt(i)%y = Top+Pt(i)%Radius
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
