Module Global_variable
	Implicit None
	Integer, Parameter ::					Myrealk = 4
	Integer ::								N_Total, N_Distribution, N_Radius, N_Second
	Integer, Allocatable ::					N_Sand(:), Mesh(:,:,:), SandInMesh(:,:,:)
	Real( kind=Myrealk ), Parameter :: 		Pi=3.14159, Kn=1.0E2, Kt=1.0E2, Rho=2.65E3, Epsl=0.35, G=-9.8, KnEdge=1.0E2  !Question
	Real( kind=Myrealk ), Parameter :: 		Miu=0.5, SmallValue=1.0E-8, Rinfinity=1.0E0
	Real( kind=Myrealk ), Allocatable ::	Radius(:), Grid(:,:), Percentage(:)
	Real( kind=Myrealk ) ::					L_Bed, H_Bed, Height, Length, MeshSize
	Type :: Particle
		Integer :: 						i
		Integer :: 						M
		Integer ::						N
		Real( kind=Myrealk ) :: 		x
		Real( kind=Myrealk ) :: 		y
		Real( kind=Myrealk ) :: 		vx
		Real( kind=Myrealk ) :: 		vy
		Real( kind=Myrealk ) :: 		radius
		Real( kind=Myrealk ) :: 		w
		Real( kind=Myrealk ) :: 		Mass
		Real( kind=Myrealk ) :: 		Inertia
	End Type
	Type( Particle ), Allocatable ::	Sand(:)

End Module Global_variable

Module MAlgorithm
	Use Global_variable
	Implicit None
	Integer ::              N_RowOfMesh, N_ColumnOfMesh

Contains
Subroutine Info_SandBed()
	Implicit None

	Print*,'------------------Infomation--------------------'
	Print*
	Print*,'This algorithm build a sand bed contains different particles'
	Print*,'Length Unit is M'
	Print*,'Array Radius must save value in order'
	Print*,'------------------Infomation--------------------'
	Print*
	!Print*,'Input the total number of sand bed '
	!Read(*,*),N_Total
	!Allocate( Sand(N_Total) )
	!Print*,'Input the number of particle sorts'
 	!Read(*,*),N_Radius
	!Allocate( Radius(N_Radius), Percentage(N_Radius), N_Sand(N_Radius) )
	!Radius = 0
	!Print*,'Input radius of particle'
	!Read(*,*), ( Radius(i),i=1,N_Radius )
	!Print*,'Input length and height of sand bed in order'
	!Read(*,*), Length, Height
	!Print*,'Which distribution of sand bed do you need ?'
	!Print*,'Normal LogNormal Uniform ? 1/2/3'
	!Read(*,*),N_Distribution

	N_Second = 20
	N_Total = 1000
	Allocate( Sand(N_Total+1) )
	N_Radius = 9
	Allocate( Radius(N_Radius), Percentage(N_Radius), N_Sand(N_Radius) )
	Radius = (/0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5/)
	Radius(:) = Radius(:)*1.0E-3
	Length = 0.045
	N_Distribution = 1
	Call Get_Particle()

End Subroutine Info_SandBed

Subroutine Get_Particle()
	Implicit None
	Integer :: 			i, j, k

	Call Distribution()
	k = 0
	Do i = 1, N_Radius
		k = k + N_Sand(i)
		Do j = k-N_Sand(i)+1, k
			Sand(j) % i = j
			Sand(j) % radius = Radius(i)
			Sand(j) % Mass = Rho*Pi*Radius(i)**3*4.0/3.0
			Sand(j) % Inertia = 0.4 * Sand(j)%Mass * Sand(j)%Radius * Sand(j)%Radius
			!Write( 11,* ), j, Sand(j)%Radius, Sand(j)%Mass, Sand(j)%Inertia
		End  Do
	End do


End Subroutine Get_Particle

Subroutine Distribution()
	Implicit None
	Integer :: 					i, Nstep, j
	Real( kind=Myrealk ) :: 	Start, Finish, Jiange, Step, X, S

	Percentage(:) = 0.0
	Select case( N_Distribution )
	Case( 1 )
		Start = -3.0
		Finish = 3.0
		Jiange = ( -Start+Finish )/Real( N_Radius )
		Step = 1.0E-5
		Nstep = Int( Jiange/Step )
	Case( 2 ) ! Question Start and finish
		Start = -3.0
		Finish = 3.0
		Jiange = ( Start-Finish )/Real( N_Radius )
		Step = 1.0E-5
		Nstep = Int( Jiange/Step )
	Case( 3 )
		Percentage(:) = 1.0/Real( N_Radius )
		Return
	End Select
	Do i = 1, N_Radius
		Do j = 1, Nstep
		    x = Start+Real(i-1)*Jiange+Real(J-1)*Step
			Percentage(i) = Percentage(i)+Step*F( x )
		End Do
	End Do
	S = Sum(Percentage)
	Do i = 1, N_Radius
		Percentage(i) = Percentage(i)/S
	 	N_Sand(i) = Int( Real(N_Total)*Percentage(i) )
	End Do
	If ( Sum(N_Sand)/=N_Total )then
		N_Sand( Int(0.5*N_Radius)+1 ) = N_Sand( Int(0.5*N_Radius)+1 ) + N_Total - Sum(N_Sand)
	End If
    Print*,'Total number is', Sum(N_Sand)

End Subroutine Distribution

Real Function F( X )
	Implicit None
	Real( kind=Myrealk ) :: 	X

	Select Case( N_Distribution )
	Case( 1 )
		F = 1.0/(2.0/Pi)**0.5*Exp( -x*x*0.5 )
	Case( 2 )
		F = 1!Question
	End Select

End Function F

Subroutine Get_Grid( R_Grid )
	Implicit None
	Integer :: 					N_Length, N_Height, i, j, k
	Real( kind=Myrealk ) ::		R_Grid, Start

	Start = 3.0*Maxval(Radius)   !Question
	R_Grid =  Maxval(Radius)
	N_Length = Int( Length*0.8/R_Grid/2.0 )
	N_Height = Int( N_Total/N_Length ) + 1
	Allocate( Grid(N_Total,2) )
	Grid = 0
	k = 0
	Do i = 1, N_Height
		Do j = 1, N_Length
			k = k + 1
			If( K>N_Total )then
				Exit
			End if
			Grid(k,1) = Real(j-1)*R_Grid*2.0 + R_Grid
			Grid(k,2) = Start + Real(i-1)*2.0*R_Grid + R_Grid
		End Do
	End Do
	Height = Grid(N_Total,2) + R_Grid

End Subroutine Get_Grid

Subroutine Set_SandBed()
	Implicit None
	Integer :: 					i, A, B
	Real( kind=Myrealk ) ::		Temp, n1, n2, V

	Call Random_Seed()
	Do i = 1, 2*N_Total
		Call Random_Number( n1 )
		Call Random_Number( n2 )
		A = Int( n1*Real(N_Total) ) + 1
		B = Int( n2*Real(N_Total) ) + 1
		Temp = Sand(A) % Radius
		Sand(A) % Radius = Sand(B) % Radius
		Sand(B) % Radius = Temp
	End Do
	Call Flush(111)
	Call Random_Seed()
	Do i = 1, N_Total
		Call Random_Number(v)
		Sand(i) % x = Grid(i,1)
		Sand(i) % y = Grid(i,2)
		Sand(i) % vx = 0.1*(v-0.5)
		Sand(i) % vy = 0.0
		Sand(i) % w = 0.0
		Sand(i) % M = 0
		Sand(i) % N = 0
	End Do
	Sand( N_Total+1 )%i = N_Total + 1
	Sand( N_Total+1 )%X = -RInfinity
	Sand( N_Total+1 )%Y = -RInfinity
	Sand( N_Total+1 )%W = 0.0
	Sand( N_Total+1 )%M = 0
	Sand( N_Total+1 )%N = 0
	Sand( N_Total+1 )%Vx = 0.0
	Sand( N_Total+1 )%Vy = 0.0
	Sand( N_Total+1 )%Mass = Rho*Pi*Radius(i)**3*4.0/3.0
	Sand( N_Total+1 )%Radius  = RInfinity
	Sand( N_Total+1 )%Inertia = 0.0

End Subroutine Set_SandBed

Subroutine Get_Mesh( R_Grid )
	Implicit None
	Integer :: 				i, j, ii, jj
	Real( kind=MyrealK ) ::	R_Grid, Rmin, Rmax

	Rmin = Minval( Radius )
	Rmax = Maxval( Radius )
	MeshSize = Rmax*2.0
	N_RowOfMesh = Int( (Grid(N_Total,2)+R_Grid)/MeshSize ) + 2
	N_ColumnOfMesh = Int( Length/MeshSize ) + 2
	Length = Real(N_ColumnOfMesh-2)*MeshSize
	Height = Real(N_RowOfMesh-2)*MeshSize
	i = Int(  MeshSize*MeshSize/Pi/Rmin/Rmin ) + 1
	Allocate( Mesh(0:N_RowOfMesh-1,0:N_ColumnOfMesh-1,i), SandInMesh(N_Total,5,2) )
	SandInMesh(:,:,:) = 0
	!Initialize Mesh, Set (N_Total+1)nth particle in Mesh out-ring.
	Mesh(:,:,:) = 0
	Mesh(0,:,1) = 1
	Mesh(0,:,2) = N_Total + 1
	Mesh(N_RowOfMesh-1,:,1) = 1
	Mesh(N_RowOfMesh-1,:,2) = N_Total + 1
	Mesh(:,0,1) = 1
	Mesh(:,0,2) = N_Total + 1
	Mesh(:,N_ColumnOfMesh-1,1) = 1
	Mesh(:,N_ColumnOfMesh-1,2) = N_Total + 1
	Print*,'Mesh size is ', N_RowOfMesh, N_ColumnOfMesh, i
	Print*,'Calculation length and height are ', Length, Height
!    Do ii = 0, N_RowOfMesh-1
!        Do jj = 0, N_ColumnOfMesh-1
!            Write( 11,* ) ii,jj,Mesh(ii,jj,1), Mesh(ii,jj,2), Mesh(ii,jj,3), Mesh(ii,jj,4), Mesh(ii,jj,5), Mesh(ii,jj,6)
!        End Do
!    End Do

End Subroutine Get_Mesh

!Init sandbed. Get SandMesh and Mesh of the start state
Subroutine InitSand()
	Implicit None
	Integer :: 					i, j, M, N, X
	Real( kind=Myrealk ) ::		Xmin, Xmax, Ymin, Ymax, Lxmin, Lxmax, Lymin, Lymax
	Real( kind=Myrealk ) ::		X1, X2, X3, X4

	Do i = 1, N_Total
        M = Int( Sand(i)%Y / MeshSize ) + 1
        N = Int( Sand(i)%X / MeshSize ) + 1
		Xmin = Sand(i) % x - Sand(i) % Radius
		Xmax = Sand(i) % x + Sand(i) % Radius
		Ymin = Sand(i) % y - Sand(i) % Radius
		Ymax = Sand(i) % y + Sand(i) % Radius
        Lxmin = Real( N-1 ) * MeshSize
        Lxmax = Real( N ) * MeshSize
        Lymin = Real( M-1 ) * MeshSize
        Lymax = Real( M ) * MeshSize
        If(Abs(Xmin-Lxmin)>SmallValue)then
        	x1 = (Xmin-Lxmin)/Abs((Xmin-Lxmin)) * 9.0
        Else
            X1 = 9.0
        End If
        If(Abs(Xmax-Lxmax)>SmallValue)then
            x2 = (Xmax-Lxmax)/Abs((Xmax-Lxmax)) * 2.0
        Else
            X2 = 2.0
        End If
        If(Abs(Ymin-Lymin)>SmallValue)then
            x3 = (Ymin-Lymin)/Abs((Ymin-Lymin)) * 1.0
        Else
            X3 = 1.0
        End If
        If(Abs(Ymax-Lymax)>SmallValue)then
            x4 = (Ymax-Lymax)/Abs((Ymax-Lymax)) * 4.0
        Else
            X4 = 4.0
        End If
		x = Int( x1 + x2 + x3 + x4 )
		Sand(i) % M = M
		Sand(i) % N = N
		SandInMesh( i,1,1 ) = 1
		SandInMesh( i,2,1 ) = M
		SandInMesh( i,2,2 ) = N
		j = Mesh(m,n,1)
		Mesh(m,n,j+2) = Sand(i)%i
		Mesh(m,n,1) = Mesh(m,n,1) + 1

        Select Case( x )
        Case( -16 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m,n-1,1)
            Mesh(m,n-1,j+2) = Sand(i)%i
            Mesh(m,n-1,1) = Mesh(m,n-1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m-1,n-1,1)
            Mesh(m-1,n-1,j+2) = Sand(i)%i
            Mesh(m-1,n-1,1) = Mesh(m-1,n-1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m-1,n,1)
            Mesh(m-1,n,j+2) = Sand(i)%i
            Mesh(m-1,n,1) = Mesh(m-1,n,1) + 1
        Case( 2 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m-1,n,1)
            Mesh(m-1,n,j+2) = Sand(i)%i
            Mesh(m-1,n,1) = Mesh(m-1,n,1) + 1
        Case( 6 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m,n+1,1)
            Mesh(m,n+1,j+2) = Sand(i)%i
            Mesh(m,n+1,1) = Mesh(m,n+1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m-1,n+1,1)
            Mesh(m-1,n+1,j+2) = Sand(i)%i
            Mesh(m-1,n+1,1) = Mesh(m-1,n+1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m-1,n,1)
            Mesh(m-1,n,j+2) = Sand(i)%i
            Mesh(m-1,n,1) = Mesh(m-1,n,1) + 1
        Case( 8 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m,n+1,1)
            Mesh(m,n+1,j+2) = Sand(i)%i
            Mesh(m,n+1,1) = Mesh(m,n+1,1) + 1
        Case( 16 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m,n+1,1)
            Mesh(m,n+1,j+2) = Sand(i)%i
            Mesh(m,n+1,1) = Mesh(m,n+1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m+1,n+1,1)
            Mesh(m+1,n+1,j+2) = Sand(i)%i
            Mesh(m+1,n+1,1) = Mesh(m+1,n+1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m+1,n,1)
            Mesh(m+1,n,j+2) = Sand(i)%i
            Mesh(m+1,n,1) = Mesh(m+1,n,1) + 1
        Case( 12 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m+1,n,1)
            Mesh(m+1,n,j+2) = Sand(i)%i
            Mesh(m+1,n,1) = Mesh(m+1,n,1) + 1
        Case( -6 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m,n-1,1)
            Mesh(m,n-1,j+2) = Sand(i)%i
            Mesh(m,n-1,1) = Mesh(m,n-1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m+1,n-1,1)
            Mesh(m+1,n-1,j+2) = Sand(i)%i
            Mesh(m+1,n-1,1) = Mesh(m+1,n-1,1) + 1
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m+1,n,1)
            Mesh(m+1,n,j+2) = Sand(i)%i
            Mesh(m+1,n,1) = Mesh(m+1,n,1) + 1
        Case( -14 )
            SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
            SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
            SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
            j = Mesh(m,n-1,1)
            Mesh(m,n-1,j+2) = Sand(i)%i
            Mesh(m,n-1,1) = Mesh(m,n-1,1) + 1
        End Select
	End Do
!	Do i = 1, N_Total
!        Write( 11,* ), SandInMesh(i,1,1), SandInMesh(i,2,:), SandInMesh(i,3,:), SandInMesh(i,4,:), SandInMesh(i,5,:)
!	End Do
!	Read( *,* )

End Subroutine InitSand

Subroutine SaveResult( Flag )
    Implicit None
    Integer ::          i, Flag
    Character ::        FileName*20

    Write( FileName,"(g0)" ) Flag
    FileName = Trim(Adjustl(FileName))//'.txt'
    Open( Flag,File=FileName )

    Do i = 1, N_Total
        Write( Flag,"(i4,2x,E13.6,2x,E13.6,2x,E13.6,2x,E13.6,2x,E13.6,2x)" )  &
            i, Sand(i)%Vx, Sand(i)%Vy, Sand(i)%X, Sand(i)%Y, Sand(i)%Radius
    End Do
	Call Flush( Flag )
	Close( Flag )

End Subroutine SaveResult

End Module MAlgorithm


Program Main
	Use Global_variable
    Use MAlgorithm
	Implicit None
	Integer :: 					i, j, k, P, M, N, i1, i2, N_Time, Time, X, Flag, Step, Nsave
	Real( kind=Myrealk ) ::  	Fx, Fy, Fw, Fnx, Fny, Ftx, Fty, Dn, Dt, Fn, Ft, Dnx, Dny, Dtx, Dty, Bn, Bt, R_Grid
	Real( kind=Myrealk ) ::  	Vi, Vj, CosAngle, SinAngle, CosAngleT, SinAngleT, VjNt, ViNt, Compression
	Real( kind=Myrealk ) ::		Xmin, Xmax, Ymin, Ymax, Lxmin, Lxmax, Lymin, Lymax, KineticEnergy
	Real( kind=Myrealk ) ::		X1, X2, X3, X4, DeltaT, VectorFlag, Vitangent, Vjtangent, Vt,mm,nn

	Call Info_SandBed()
	Call Get_Grid( R_Grid )
	Call Set_SandBed()
	Call Get_Mesh( R_Grid )
	Call InitSand()
	Call SaveResult( 0 )
	!start settle calculation
	DeltaT = 1.0E-6
	N_Time = Int(DeltaT**(-1)*N_Second)
	Step = DeltaT**(-1)*1E-2
	Nsave = 0
    Do Flag = 1, N_Time, Step!Loop Time
        Do Time = Flag, Flag+Step-1
            !Calculation of snd movement
            Do i = 1, N_Total!Loop calculation of movement
                k = SandInMesh( i,1,1 )
                If( k<1 )then
                    Print*,'SandInMesh',i,'maybe wrong'
					Read( *,* )
                    Exit
                End If
                Fx = 0.0
                Fy = 0.0
                Fw = 0.0
                Fnx = 0.0
                Fny = 0.0
                Ftx = 0.0
                Fty = 0.0
                Dnx = 0.0
                Dny = 0.0
                Dtx = 0.0
                Dty = 0.0
                Do i1 = 2,k+1!Loop SandInMesh
                    !Calculation of compression of sand i and sand j
                    M = SandInMesh( i,i1,1 )
                    N = SandInMesh( i,i1,2 )
                    p = Mesh( M,N,1 )
					If( P<1 )then
						Print*,'Mesh(M,N) has no sand, something is wrong'
						Read( *,* )
                    	Exit
					End If
                    Do i2 = 2, p+1
                        j = Mesh( M,N,i2 )
                        If( j==i )then
                            Cycle
                        End if
                        If( j==N_Total+1 )then
                            If( M==0 )then
                            	If( N==0 )then
									Fnx =  KnEdge * Abs(Sand(i)%x-Sand(i)%Radius) + Fnx
									Fny =  KnEdge * Abs(Sand(i)%y-Sand(i)%Radius) + Fny
								ElseIf( N==N_ColumnOfMesh-1 )then
									Fnx = - KnEdge * Abs(Sand(i)%x+Sand(i)%Radius-Length) + Fnx
									Fny =  KnEdge * Abs(Sand(i)%y-Sand(i)%Radius) + Fny
								Else
									Fnx = 0.0 + Fnx
									Fny =  KnEdge * Abs(Sand(i)%y-Sand(i)%Radius) + Fny
								End If
							ElseIf( M==N_RowOfMesh-1 )then
                            	If( N==0 )then
									Fnx =  KnEdge * Abs(Sand(i)%x-Sand(i)%Radius) + Fnx
									Fny = - KnEdge * Abs(Sand(i)%y+Sand(i)%Radius-Height) + Fny
								ElseIf( N==N_ColumnOfMesh-1 )then
									Fnx = - KnEdge * Abs(Sand(i)%x+Sand(i)%Radius-Length) + Fnx
									Fny = - KnEdge * Abs(Sand(i)%y+Sand(i)%Radius-Height) + Fny
								Else
									Fnx = 0.0 + Fnx
									Fny = - KnEdge * Abs(Sand(i)%y+Sand(i)%Radius-Height) + Fny
								End If
                            Else
                            	If( N==0 )then
									Fny = 0.0 + Fny
									Fnx =  KnEdge * Abs(Sand(i)%x-Sand(i)%Radius) + Fnx
								Elseif( N==N_ColumnOfMesh-1 )then
									Fny = 0.0 + Fny
									Fnx = - KnEdge * Abs(Sand(i)%x+Sand(i)%Radius-Length) + Fnx
								End If
                            End If
							Ftx = Ftx
							Fty = Fty
							Dtx = Dtx
							Dty = Dty
							Dnx = Dnx
							Dny = Dny
							Fw = Fw
						Else
                        Compression = Sand(i)%Radius+Sand(j)%Radius - ( (Sand(i)%X-Sand(j)%X)**2 + (Sand(i)%Y-Sand(j)%Y)**2 )**0.5
                        If( Compression>0.0 )then
                            !Calculation of Normal Force
                            Fn = Kn * Compression
                            CosAngle = ( Sand(i)%X-Sand(j)%X ) / ( (Sand(i)%X-Sand(j)%X)**2 + (Sand(i)%Y-Sand(j)%Y)**2 )**0.5
                            SinAngle = ( Sand(i)%Y-Sand(j)%Y ) / ( (Sand(i)%X-Sand(j)%X)**2 + (Sand(i)%Y-Sand(j)%Y)**2 )**0.5
                            Fnx = Fn * CosAngle + Fnx
                            Fny = Fn * SinAngle + Fny
                            !Calculation of Tangential Force
                            !Tangential vector  Nt = ( Sand(i)%Y-Sand(j)%Y, Sand(j)%x-Sand(i)%x )
                            !Velocity vector Vi = ( Sand(i)%Vx, Sand(i)%Vy ) Vj = ( Sand(j)%Vx, Sand(j)%Vy )
                            !Decompostion vector(Vi Vj) to tangent dircetion Nt
                            Vi = ( Sand(i)%Vx**2 + Sand(i)%Vy**2 )**0.5
                            Vj = ( Sand(j)%Vx**2 + Sand(j)%Vy**2 )**0.5
                            ViNt = ( Sand(i)%Vx*(Sand(i)%Y-Sand(j)%Y) + Sand(i)%Vy*(Sand(j)%X-Sand(i)%X) )
                            VjNt = ( Sand(j)%Vx*(Sand(i)%Y-Sand(j)%Y) + Sand(j)%Vy*(Sand(j)%X-Sand(i)%X) )
                            Vitangent = Vi * ViNt / ( Vi  * ((Sand(i)%X-Sand(j)%X)**2+(Sand(j)%Y-Sand(i)%Y)**2)**0.5 )
                            Vjtangent = Vj * VjNt / ( Vj  * ((Sand(i)%X-Sand(j)%X)**2+(Sand(j)%Y-Sand(i)%Y)**2)**0.5 )
                            Vt = Vjtangent - Vitangent
							If( Vt/=0.0 )then
                            	Ft = Miu * Fn * Vt/Abs(Vt)  !Question
							Else
								Ft = 0.0
							End If
                            CosAngleT = ( Sand(i)%Y-Sand(j)%Y ) / ( (Sand(i)%X-Sand(j)%X)**2+(Sand(j)%Y-Sand(i)%Y)**2 )**0.5
                            SinAngleT = ( Sand(j)%X-Sand(i)%X ) / ( (Sand(i)%X-Sand(j)%X)**2+(Sand(j)%Y-Sand(i)%Y)**2 )**0.5
                            Ftx = Ft * CosAngleT + Ftx
                            Fty = Ft * SinAngleT + Fty
                            !Calculatin of damping force( Tangential and normal direction )
                            Bn = ( 4.0*Kn*Sand(i)%Mass*log(Epsl)*log(Epsl)/(Pi*Pi-log(Epsl)*log(Epsl)) )**0.5
                            Bt = ( 4.0*Kt*Sand(i)%Mass )**0.5
                            Dn = -Bn * Compression / DeltaT
                            Dt = -Bt * Vt
                            Dnx = Dn * CosAngle + Dnx
                            Dny = Dn * SinAngle + Dny
                            Dtx = Dt * CosAngleT + Dtx
                            Dty = Dt * SinAngleT + Dty
							If( Vt/=0.0 )then
                            	Fw = Fw + Ft*Vt/Abs(Vt) + Dt
							Else
								Fw = Fw 
							End If
                        End If
                        End IF
                    End Do
				End Do !Loop SandInMesh
				If( Sand(i)%Radius>Sand(i)%Y )then
					Sand(i)%Vx = 0.9*Sand(i)%Vx
					Sand(i)%Vy = 0.9*Sand(i)%Vy
				End If
				Fx = Ftx + Fnx! + Dnx + Dtx
				Fy = Fty + Fny! + Dny + Dty
                Sand(i)%Vx = Sand(i)%Vx + Fx/Sand(i)%Mass * Deltat
                Sand(i)%Vy = Sand(i)%Vy + Fy/Sand(i)%Mass * Deltat + G*Deltat
                Sand(i)%X = Sand(i)%X + Sand(i)%Vx * DeltaT + 0.5*Fx/Sand(i)%Mass*DeltaT**2
                Sand(i)%Y = Sand(i)%Y + Sand(i)%Vy * DeltaT + 0.5*( Fy/Sand(i)%Mass+G )*DeltaT**2
                !Sand(i)%X = Sand(i)%X + Fx/Sand(i)%Mass * DeltaT**2
                !Sand(i)%Y = Sand(i)%Y + ( Fy/Sand(i)%Mass+G ) * DeltaT**2
                Sand(i)%W = Sand(i)%W + Fw*Sand(i)%Radius*DeltaT/Sand(i)%Inertia
            End Do !Loop calculation of movement
            !Update Mesh and SandInMesh
            Mesh( :,:,: ) = 0
            Mesh(0,:,1) = 1
            Mesh(0,:,2) = N_Total + 1
            Mesh(N_RowOfMesh-1,:,1) = 1
            Mesh(N_RowOfMesh-1,:,2) = N_Total + 1
            Mesh(:,0,1) = 1
            Mesh(:,0,2) = N_Total + 1
            Mesh(:,N_ColumnOfMesh-1,1) = 1
            Mesh(:,N_ColumnOfMesh-1,2) = N_Total + 1
            SandInMesh( :,:,: ) = 0
            Do i = 1, N_Total!Loop Update Mesh and SandInMesh
                M = Int( Sand(i)%Y / MeshSize ) + 1
                N = Int( Sand(i)%X / MeshSize ) + 1
                Lxmin = Real( N-1 ) * MeshSize
                Lxmax = Real( N ) * MeshSize
                Lymin = Real( M-1 ) * MeshSize
                Lymax = Real( M ) * MeshSize
                Xmin = Sand(i)%x - Sand(i)%Radius
                Xmax = Sand(i)%x + Sand(i)%Radius
                Ymin = Sand(i)%y - Sand(i)%Radius
                Ymax = Sand(i)%y + Sand(i)%Radius
                If(Abs(Xmin-Lxmin)>SmallValue)then
                    x1 = (Xmin-Lxmin)/Abs((Xmin-Lxmin)) * 9.0
                Else
                    X1 = 9.0
                End If
                If(Abs(Xmax-Lxmax)>SmallValue)then
                    x2 = (Xmax-Lxmax)/Abs((Xmax-Lxmax)) * 2.0
                Else
                    X2 = 2.0
                End If
                If(Abs(Ymin-Lymin)>SmallValue)then
                    x3 = (Ymin-Lymin)/Abs((Ymin-Lymin)) * 1.0
                Else
                    X3 = 1.0
                End If
                If(Abs(Ymax-Lymax)>SmallValue)then
                    x4 = (Ymax-Lymax)/Abs((Ymax-Lymax)) * 4.0
                Else
                    X4 = 4.0
                End If
                x = Int( x1 + x2 + x3 + x4 )
                Sand(i) % M = M
                Sand(i) % N = N
                SandInMesh( i,1,1 ) = 1
                SandInMesh( i,2,1 ) = M
                SandInMesh( i,2,2 ) = N
                j = Mesh(m,n,1)
                Mesh(m,n,j+2) = Sand(i)%i
                Mesh(m,n,1) = Mesh(m,n,1) + 1
                Select Case( x )
                Case( -16 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m,n-1,1)
                    Mesh(m,n-1,j+2) = Sand(i)%i
                    Mesh(m,n-1,1) = Mesh(m,n-1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m-1,n-1,1)
                    Mesh(m-1,n-1,j+2) = Sand(i)%i
                    Mesh(m-1,n-1,1) = Mesh(m-1,n-1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m-1,n,1)
                    Mesh(m-1,n,j+2) = Sand(i)%i
                    Mesh(m-1,n,1) = Mesh(m-1,n,1) + 1
                Case( 2 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m-1,n,1)
                    Mesh(m-1,n,j+2) = Sand(i)%i
                    Mesh(m-1,n,1) = Mesh(m-1,n,1) + 1
                Case( 6 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m,n+1,1)
                    Mesh(m,n+1,j+2) = Sand(i)%i
                    Mesh(m,n+1,1) = Mesh(m,n+1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m-1,n+1,1)
                    Mesh(m-1,n+1,j+2) = Sand(i)%i
                    Mesh(m-1,n+1,1) = Mesh(m-1,n+1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M-1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m-1,n,1)
                    Mesh(m-1,n,j+2) = Sand(i)%i
                    Mesh(m-1,n,1) = Mesh(m-1,n,1) + 1
                Case( 8 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m,n+1,1)
                    Mesh(m,n+1,j+2) = Sand(i)%i
                    Mesh(m,n+1,1) = Mesh(m,n+1,1) + 1
                Case( 16 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m,n+1,1)
                    Mesh(m,n+1,j+2) = Sand(i)%i
                    Mesh(m,n+1,1) = Mesh(m,n+1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N+1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m+1,n+1,1)
                    Mesh(m+1,n+1,j+2) = Sand(i)%i
                    Mesh(m+1,n+1,1) = Mesh(m+1,n+1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m+1,n,1)
                    Mesh(m+1,n,j+2) = Sand(i)%i
                    Mesh(m+1,n,1) = Mesh(m+1,n,1) + 1
                Case( 12 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m+1,n,1)
                    Mesh(m+1,n,j+2) = Sand(i)%i
                    Mesh(m+1,n,1) = Mesh(m+1,n,1) + 1
                Case( -6 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m,n-1,1)
                    Mesh(m,n-1,j+2) = Sand(i)%i
                    Mesh(m,n-1,1) = Mesh(m,n-1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m+1,n-1,1)
                    Mesh(m+1,n-1,j+2) = Sand(i)%i
                    Mesh(m+1,n-1,1) = Mesh(m+1,n-1,1) + 1
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M+1
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m+1,n,1)
                    Mesh(m+1,n,j+2) = Sand(i)%i
                    Mesh(m+1,n,1) = Mesh(m+1,n,1) + 1
                Case( -14 )
                    SandInMesh( i,SandInMesh(i,1,1)+2,1 ) = M
                    SandInMesh( i,SandInMesh(i,1,1)+2,2 ) = N-1
                    SandInMesh( i,1,1 ) = SandInMesh( i,1,1 ) + 1
                    j = Mesh(m,n-1,1)
                    Mesh(m,n-1,j+2) = Sand(i)%i
                    Mesh(m,n-1,1) = Mesh(m,n-1,1) + 1
                End Select
            End Do!Loop Update Mesh and SandInMesh
        End Do
        Print*,'Save information, Time is', Time-1
		Nsave = Nsave + 1
        Call SaveResult( Nsave )
    End Do!Loop Time
    Print*,'Settle has down. Enter any key to exit'
    Read( *,* )

End Program Main
