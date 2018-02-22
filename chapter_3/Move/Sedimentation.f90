!//-----------------------------------------------------------------------//
!// 
!
!
!
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
