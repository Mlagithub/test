Module Evolution
	Use GlobalVariable
	Use FileOperate
	Use MeshGrid
	Use Collision
	Use ScheduleBar
	Use Wind
	Implicit None

Contains


!// Calculation sand bed splashing
Subroutine Salation( Sw_W,Sw_E )
	Implicit None
	Logical, Intent(In) ::			Sw_W !// Switch of weather couping wind and sand
	Logical, Intent(In) ::			Sw_E !// Switch of weather consider electric field force of sand movement
	Logical ::						Boder(2)
	Integer ::						Flag, Time, Bh=0, i
	Type( CLS_CMD_Progress ) ::		Progress



	!// Initialize 	
	Call Initialization( Sw_W,Sw_E )
	Call SaveResult( FileFlag=0,N=N_Total )



	!// Mesh grid
	MeshSize = 2.0*MaxRadius
	Call Get_Mesh_2D( N=N_Total,Boder_Left=0.0,Boder_Right=Edgeright,Boder_lower=0.0,&
					 & Boder_upper=Height,Rmin=MinRadius,MeshSize=MeshSize )
	Call MeshUpdate2D( Boder_Left=0.0,Boder_lower=0.0 )
	Boder = (/.True.,.False./)
	Call PeriodicBoundary2D( Boder )


	!// Schedule bar
	Call Progress%Set( N=Int(N_Second/DeltaT),L=30 )
	Progress%Prefix = "Step : "
	Progress%m = "#"
	Progress%o = "."


	Open( 111,File='Ustar.txt',Access='Append' )
	Do Flag = 1, Int(N_Second/DeltaT), Int(1.0E-2/DeltaT)
		Do Time = Flag, Flag+Int(1.0E-2/DeltaT)-1
			Pt_Old = Pt
			Call Orbit( Sw_W,Sw_E )
			Call MeshUpdate2D( Boder_Left=0.0,Boder_lower=0.0 )
			Call PeriodicBoundary2D( Boder )
		End Do 

		Call Progress%Put( k=Time-1, bAbsol=.True. )
		!// Couping wind and sand 
		If( Sw_W )then
			Call ModifyWind()
			Write(111,"(F13.6)") Ustar
		End If
	

		!// Save sand information 
		Bh = Bh + 1
		Call SaveResult( FileFlag=Bh,N=N_Total )
	End Do 




End Subroutine Salation



!// Initializate variable 
Subroutine Initialization( Sw_W,Sw_E )
	Implicit None
	Logical, Intent(In) ::			Sw_W !// Switch of weather couping wind and sand
	Logical, Intent(In) ::			Sw_E !// Switch of weather consider electric field force of sand movement
	Character ::					FileName*512, A
	Integer( kind=MyIk ) ::			i, Err, j
	Real( kind=MyRk ) ::			P(10)


	!// Read variable form Configuration file
	Call ReadConfigure()
	Call GetFileN( 'SandBed',i )
	If( i/=N_Total-N_Incident )then
		Print*,'SandBed and Configure dose not match'
		Stop
	End If



	!// Initializate sand information 
	Allocate( Pt(N_Total),Pt_old(N_Total) )
    Pt(:)%i = 0
    Pt(:)%M = 0
    Pt(:)%N = 0



	!// Initialize sand bed and incident sands
	!// Read in sand bed from local file SandBed
    Open( 11,file='SandBed',Status='old',Iostat=Err )!SandBed save particle location and radius only
	H_SandBed = 0.0
	Do i = N_Incident+1, N_Total
		Read( 11,*,Iostat=Err )  A, Pt(i)%X,  Pt(i)%Y, Pt(i)%Radius
		If( Pt(i)%Radius<1.0E-19 )then
			Print*,Pt(i)%Radius,i
			stop
		End If
		If( Err /=  0 )then
			Print*, 'Read file wrong, maybe file format is not not confirm'
			Stop
		End If
		Pt(i)%i = i
		Pt(i)%Vx = 1.0E-9
		Pt(i)%Vy = 1.0E-9
		Pt(i)%V = ( Pt(i)%Vx**2 + Pt(i)%Vy**2 )**0.5
		Pt(i)%Mass = (4.0*Rho_Sand*Pi*(Pt(i)%Radius**3))/3.0
		Pt(i)%Inertia = 8.0*Pt(i)%Radius**5*Pi*Rho_Sand/15.0
		Pt(i)%W = 0.0
		If( Pt(i)%Y + Pt(i)%Radius  > H_SandBed  ) H_SandBed = Pt(i)%Y + Pt(i)%Radius
	End Do


	!// Max and min radius
	MinRadius = Minval( Pt(N_Incident+1:N_Total)%Radius )
	MaxRadius = Maxval( Pt(N_Incident+1:N_Total)%Radius )
	Print*,H_SandBed
	!H_SandBed = Min( H_SandBed,0.006 )



	!// Randomly generated incident sands
	Call Random_Seed()
	i = 0
	Do while( .True. )
		i = i + 1
		If( i>N_Incident ) Exit
		Call Random_Number( P )
		If( N_Sorts==N_Total )then
			Pt(i)%Radius = P(1)*(MaxRadius-MinRadius) + MinRadius
		Else
			Pt(i)%Radius = Radius( Int(P(1)*N_Sorts)+1 )
		End If
		If( Pt(i)%Radius<1.0E-19 )then
			Print*,Pt(i)%Radius,i
			stop
		End If
		Pt(i)%i = i
		Pt(i)%X = P(2)*Edgeright
		Pt(i)%Y = EdgeTop + ( 5.0*P(5)+1 ) * Pt(i)%Radius
		Pt(i)%Vx = P(3)
		Pt(i)%Vy = -1.0 - P(4)
		Pt(i)%V = ( Pt(i)%Vx**2 + Pt(i)%Vy**2 )**0.5
		Pt(i)%Mass = (4.0*Rho_Sand*Pi*(Pt(i)%Radius**3))/3.0
		Pt(i)%Inertia = 8.0*Pt(i)%Radius**5*Pi*Rho_Sand/15.0
		Pt(i)%W = 0.0
		Do j = 1, i-1
			If( Pt(i)%Radius+Pt(j)%Radius > &
				Distance( x1=Pt(i)%X,y1=Pt(i)%Y,x2=Pt(j)%x,y2=Pt(j)%Y) )then
				i = i -1
				Exit
			End If
		End Do 
	End Do 




	!// Step length of simulation
	N_Second = 20
	DeltaT = 1.0E-6



	!// Area of calculation 
	Height = max( 0.5,H_SandBed*4.0 )
	Length = Edgeright



	!// Wind variable
	If( Sw_W )then
		Ustaro = 0.5
		Ustar = Ustaro
		Tau = Rho_Air*Ustar**2
		UstarT = 0.08*sqrt((-G*0.01E-3*(Rho_Sand-Rho_Air))/Rho_Air)
		TauT = UstarT**2 * Rho_Air
		Tau_Z0 = Tau - TauT
		Delta_Z = 0.01
		Delta_Z_Acc = 0.001
		Height_Z_Acc = 0.02
		N_Z_Acc = 20
		Call Initialize_Wind()
	End If
		

	!// Electric field
	Allocate( ComMax(N_Total,N_Total) )
	ComMax = 0.0
        Pt(:)%Co = 0
	Pt(:)%Q = 0.0
	Pt(:)%Alpha = 0.132
	Pdd = Pd( Ttemp=Temperature,Rh=Rh,Xfn=Xfn )
	Pt(:)%Pd = Pdd

	Call SaveConfigure()

End Subroutine Initialization




End Module Evolution
