 
!****************************************************************************************************************** 
!请注意  此版本完成的模块为： 周期边界条件  弹塑性受力模型  碰撞电只粗略地考虑 Alpha 为 0.5
!   
!   bt = 碰撞角度(位置）   tr = 小球半径   st = 撞击角度（速度）  ve1 = 入射速度  Gama(n,n) = 两小球的瞬时压缩量  r(n,n) = 两球心瞬时间距
!   Gama1(n,n) = 两小球上一时刻的压缩量    deltg(n,n) = 压缩量变化量
!
!****************************************************************************************************************** 
    
Module Qiu_li
    Implicit none
    Integer,public,parameter :: 				MyRealk=Selected_Real_Kind(6,36), MyIntegerK=Selected_Int_Kind(8)
    Integer( kind=MyIntegerK ) :: 				N_Time, N_Second, N_Total, N_RowOfMesh, N_ColumnOfMesh, N_Radius, N_Incident, Mx, Step, Time
    Integer( kind=MyIntegerK ) :: 				N_CellofQHeight, Bh=10, N_Z, N_Z_Acc
	Real(Kind = Myrealk) :: 					DeltaT, Length, Height, MeshSize, Compression, Pd1, MinRadius, MaxRadius, F_DragX, F_DragZ, E
	Real(Kind = Myrealk) :: 					Fx, Fy, Fn, Ft, Bn, Bt, Fnx, Fny, Ftx, Fty, Bnx, Bny, Dn, Dt, Dnx, Dny, Dtx, Dty, Fw
	Real(Kind = Myrealk) :: 					Gao,GaoE,GaoW, Kn, Kt, KnEdge, Height_ElectricField, CellOfElectricField, Temperature
	Real(Kind = Myrealk) :: 					Delta_Z, Delta_Z_Acc, Tau, Tau_Z0, Height_Z_Acc , Ustar
	Real(Kind = Myrealk) :: 					Edgeright, Edgetop, Ee
	Real(Kind = 8 ) :: 							Pai
    Real(kind = MyRealk), Parameter :: 			Pi=3.14159, Miu=0.5, G=-9.8, Den=2.65E3, Epsl=0.35, Miv=0.3
    Real(kind = MyRealk), Parameter :: 			CellOfQHeight = 0.01, Rho_Air=1.29, Ga=1.5E-5, Kaman=0.4, Z0=9.3E-4
	Real(kind = MyRealk), Allocatable :: 		Incident_R(:), Incident_V(:), Incident_A(:), ComMax(:,:), Q_Height(:), ElectricField(:)
	Real(kind = MyRealk), Allocatable :: 		TauA(:), TauP(:), Zwind(:), U_New(:), R(:), Probability(:), Qm(:), Eee(:)
	Integer( kind=MyIntegerK ), Allocatable ::	Mesh(:,:,:), SandInMesh(:,:)
	Type :: Particles
	    Integer( kind=MyIntegerK ) :: 			i, M, N, Co, OrderofRadius
	    Real(Kind = MyRealk) :: 				X, Y, Radius, Vx, Vy, V, Alpha, Q, Mass, Pd!, W, Inertia
	End Type Particles
	Type(Particles), Allocatable :: 			Pt(:), Pt_Old(:)
	
Contains

!//Solve function integral
Subroutine Integral(Inte,T)
	Implicit None
	Real( Kind=8 ) ::  				X0, X1, X, Bc, Inte, Pai1, T
	Real( kind=8 ), Parameter :: 	Zeta0=3.2E-8, Eps0=8.85E-12 ,Kb=1.38E-23, Q=1.6E-19
	Inte = 0.0
	X0 = 1.0
	X1 = 1.0E5
	Bc = 1.0E-2
	X = X0
	Do 
		X = X + Bc
		Pai1 = Exp((-Q*Q*X)/(4*Pi*Zeta0*Eps0*Kb*T))/X**2
		If(X >= X1) Exit
		Inte = Inte + Pai1*Bc
	End Do
End Subroutine Integral


!//Get rows of file
Subroutine Getfilen(Flag,Hs)
	Implicit None
	Integer( kind=MyIntegerK ) :: 	Hs, Flag, Err
	Character :: 					A

    Open(11,file ='SandBed',Status='old',Iostat=Err )!SandBed save particle location and radius only
	If( Err/=0 )then
		Print*,'Can not find file SandBed'
		Read(*,*)
		Stop
	End If
	Hs = 0
	Rewind(Flag)
	Do 
		Read(Flag, *, Iostat = Err) A
		If(Err /=  0) Exit
		Hs = Hs + 1
	End Do 
	Rewind(Flag)
	Close( 11 )

End Subroutine Getfilen


!//Read original module information(particle position and randius)
Subroutine Ke_Li_Bai_Fang()
    Implicit none
	Integer( kind=MyIntegerK ) ::				Err, i, a, NGrid, j
	Real( Kind=MyRealk ) :: 					Jd, Kuan, SmallValue, ri, Seed(3)
	Character ::								Ch*40
	
	!//---------------Read radius sorts number from configure--------------//
	Open( 999,File='Configure',Status='old',Iostat=Err )
	If( Err/= 0 )then
		Write(*,"(A40)") 'Cant not find file Configure '
	End If
	Rewind( 999 )
	Read( 999,* ) Ch
	Read( 999,* ) Ch
	Read( 999,* ) N_Radius
	Read( 999,* ) ri
	Read( 999,* ) Edgetop
	Read( 999,* ) ri
	Read( 999,* ) Edgeright
	Allocate( Eee(N_Radius),R(N_Radius),Probability(N_Radius) )
	Read( 999,* ) R(:)
	Read( 999,* ) Probability(:)
	close( 999 )
	!Eee = (/20,25,40,45/)
	!Eee = 30
	Eee = (/20,70,90,110,130,150,200/)
	MaxRadius = Maxval( R )
	MinRadius = Minval( R )
	!//---------------Set the sandbed particle(location and radius)---------//
	Gao = 0.0
	Kuan = 0.0
	SmallValue = MinRadius*0.01
    Open(11,file ='SandBed',Status='old',Iostat=Err )!SandBed save particle location and radius only
	Do i = N_Incident+1, N_Total
		Read(11,*,Iostat = Err)  a, Pt(i)%X,  Pt(i)%Y, Pt(i)%Radius
		If( Err /=  0 )then
			Print*, 'Read file wrong, maybe file format is not not confirm'
			Read(*,*)
			Stop
		End If
		Pt(i)%i = i
		Do j = 1,N_Radius
			If( Abs(Pt(i)%Radius-R(j))<SmallValue )then
				Pt(i)%OrderofRadius = j
				Exit
			End If
		End Do
		Pt(i)%Vx = 1.0E-9
		Pt(i)%Vy = 1.0E-9
		Pt(i)%V = ( Pt(i)%Vx**2 + Pt(i)%Vy**2 )**0.5
		Pt(i)%Mass = (4.0*den*Pi*(Pt(i)%Radius**3))/3.0
		!Pt(i)%Inertia = 0.4*Pt(i)%Mass*Pt(i)%Radius**2
		If( Pt(i)%X + Pt(i)%Radius  >  Kuan ) Kuan = Pt(i)%X + Pt(i)%Radius
		If( Pt(i)%Y + Pt(i)%Radius  >  Gao  ) Gao = Pt(i)%Y + Pt(i)%Radius
	End Do
	Close( 11 )
	Length = Kuan
	If( Gao>Height )then
		Print*,'Original sandbed maybe wrong. SandBed Height is too large'
	End If
	!Gao = 0.9*Gao
	GaoE = 0.8*Gao
	GaoW = 0.4*Gao
	!//---------------Set the incident particle-----------------------------//
	NGrid = Int( 0.8*Length/MaxRadius/10.0 )
	Call Random_Seed()
	Do i = 1, N_Incident
		Call Random_Number( Seed )
		Incident_R(i) = R(Int(Seed(1)/Real(N_Radius))+1)
		!Incident_R(i) = 11.0E-5
		Incident_V(i) = Seed(2)*3.0
		Incident_A(i) = Seed(3)*Pi*0.25
	End do 
	Do i = 1, N_Incident
		!Jd = Incident_A(i)
		Jd = 0.17*Pi
		Pt(i)%i = i
		Pt(i)%X = 0.1*Length + Real(i-Int(i/(NGrid+0.1))*NGrid-0.5,MyRealk)*MaxRadius*10.0
		Pt(i)%Y = Gao + Pt(i)%Radius + Real(Int(i/(NGrid+0.1)),MyRealk)*MaxRadius*10.0 + 2.0*MaxRadius
		Pt(i)%Vx = 2.0*Abs(Cos(Jd))
		Pt(i)%Vy = -3.0*Abs(Sin(Jd))
		Pt(i)%V = ( Pt(i)%Vx**2 + Pt(i)%Vy**2 )**0.5
		Pt(i)%Radius = Incident_R(i)
		Pt(i)%Mass = ( 4.0*den*Pi*(Pt(i)%Radius**3) )/3.0
		!Pt(i)%Inertia = 0.4*Pt(i)%Mass*Pt(i)%Radius**2
		Do j = 1,N_Radius
			If( Abs(Pt(i)%Radius-R(j))<SmallValue )then
				Pt(i)%OrderofRadius = j
				Exit
			End If
		End Do
	End Do
	!//---------------Set the incident particle-----------------------------//
	Open( 413,File='Original' )
	Do i = 1, N_Total
		write(413,"(i5,2x,5(E13.6,3x))" ) Pt(i)%i, Pt(i)%X, Pt(i)%Y, Pt(i)%Radius, Pt(i)%vx, Pt(i)%Vy
	End Do 
	close( 413 )
	Write(*,"(A40,i13,i13)") 'N_Radius and N_Total:', N_Radius, N_Total
	Write(Ch,*) N_Radius
	Ch = '(A40,'//Trim(Adjustl(Ch))//'(F13.6,2x))'
	Write(*,ch) 'Radius:',R
	Write(*,ch) 'Probability distribution:',Probability
	Write(*,"(A40,2(F13.6,2x))") 'Edgeright and Edgetop:',Edgeright,Edgetop
	Write(*,"(A40,2(F13.6,2x))") 'Height and length of sandbed:', Gao, Kuan

End subroutine Ke_Li_Bai_Fang


!//Variable Process
Subroutine SetVariable()
	Implicit None
	Integer( kind=MyIntegerK ) :: 				Hs, i
	Real( kind=MyRealk ) :: 					Jd, UstarT, TauT
	
	!//---------------Confirm kind value of the program-----------------------//
	Write(*,"(A40,2(i13))") 'Myrealk and MyIntegerK:', MyRealk, MyIntegerK
	If( MyRealk<0 .or. MyIntegerK<0 )then
		Write(*,"(A40,2(i13))") 'Myrealk and MyIntegerK: ', MyRealk, MyIntegerK
		Read(*,*)
	End IF
	!//---------------Sandbed height and length-------------------------------//
	Length = 0.0
	Height = 0.5
	!//---------------Get total number of sandbed-----------------------------//
    Call Getfilen(11, Hs)
	N_Incident = 100
	N_Total = Hs + N_Incident
	Allocate( Pt(N_Total), Pt_Old(N_Total) )
	Allocate( ComMax( N_Total,N_Total ) )
	Allocate( Incident_R(N_Incident), Incident_V(N_Incident), Incident_A(N_Incident)) 
	!//---------------Initialition of particle information--------------------//
	Call ArgumentofElectric()
	ComMax( :,: ) = 0.0
    Pt(:)%X = 0.0
    Pt(:)%Y = 0.0
    Pt(:)%Vx = 0.0
    Pt(:)%Vy = 0.0
    Pt(:)%V = 0.0
    Pt(:)%Mass = 0.0
    !Pt(:)%Inertia = 0.0
    !Pt(:)%W = 0.0
	Pt(:)%Alpha = 0.132
	Pt(:)%Q = 0.0
    Pt(:)%i = 0
    Pt(:)%M = 0
    Pt(:)%N = 0
    Pt(:)%Co = 0
	Pt(:)%Pd = Pd1
	Call Ke_Li_Bai_Fang()
	!//---------------Set step length of simulation---------------------------//
	N_Second = 20
	DeltaT = 1.0E-6
	N_Time = DeltaT**(-1)*N_Second
	Step = Int(1.0E-2/DeltaT) !100 steps every second 
	!//---------------Set wind variables-------------------------------------//
	Ustar = 0.5
    Tau = Rho_Air*Ustar**2
	UstarT = 0.08*sqrt((-G*0.01E-3*(Den-Rho_Air))/Rho_Air)
	TauT = UstarT**2 * Rho_Air
	Tau_Z0 = Tau - TauT
	Delta_Z = 0.01
	Delta_Z_Acc = 0.001
	Height_Z_Acc = 0.02
	N_Z_Acc = 20
	Call Set_Tau()
	!//---------------Calculation of electric argument------------------------//
	N_CellofQHeight = Int(Height/CellOfQHeight)
	Height_ElectricField = Gao
	CellOfElectricField = 0.01
	Allocate( ElectricField( Int((Height-GaoE)/CellOfElectricField) ) )
	ElectricField(:) = 0.0
	!//---------------Calculation of Q/M ------------------------------------//
	Allocate( Qm(Int((Height-Gao)/Delta_Z_Acc)) )
	Qm = 0.0

End Subroutine SetVariable

Subroutine Set_Tau()
	Implicit None
	Integer( kind=MyIntegerK ) :: i

	N_Z = Int( (Height-N_Z_Acc*Delta_Z_Acc-Gao)/Delta_Z )
	Allocate( Zwind(N_Z+N_Z_Acc),TauA(N_Z_Acc+N_Z), TauP(N_Z+N_Z_Acc),U_New(N_Z+N_Z_Acc) )
	TauA(:) = 0.0
	TauP(:) = 0.0
	U_New(:) = 0.0
	Do i = 1, N_Z_Acc
		Zwind(i) = Real( Delta_Z_Acc*i )
	End Do 
	Do i = 1+N_Z_Acc, N_Z+N_Z_Acc
		Zwind(i) = Real( Delta_Z*Float(i-20) ) + Zwind(N_Z_Acc)
	End Do 

End Subroutine Set_Tau
	
Subroutine Get_Mesh()
	Implicit None
	Integer( kind=MyIntegerK ) :: 					i,j

	MeshSize = MaxRadius*2.0
	N_ColumnOfMesh = Int( Length/MeshSize ) + 2
	MeshSize = Length/Real(N_ColumnOfMesh-2)
	N_RowOfMesh = Int( Height/MeshSize ) + 2
	i = Int( MeshSize**2/Pi/MinRadius**2 ) + 1
	Allocate( Mesh(0:N_RowOfMesh-1,0:N_ColumnOfMesh-1,i+1),SandInMesh(N_Total,2) )
	Write(*,"(A40,3(i13,2x))")'N_RowOfMesh,N_ColumnOfMesh,Number:', N_RowOfMesh, N_ColumnOfMesh,i 

End Subroutine Get_Mesh


Subroutine UpdateParticleInformation()
	Implicit None
	Integer( kind=MyIntegerK ) :: 				i, j, M, N
	
	Mesh(:,:,:) = 0
	SandInMesh( :,: ) = 0
	Do i = 1, N_Total!Loop Update Mesh and SandInMesh
		M = Floor( Pt(i)%Y / MeshSize ) + 1
		N = Floor( Pt(i)%X / MeshSize ) + 1
		If( N<1 )then
			N=N_ColumnOfMesh-2
			Pt(i)%x = Pt(i)%x + Length
		End If
		If( N>N_ColumnOfMesh-2 )then
			N=1
			Pt(i)%x = Pt(i)%x - Length
		End If
		If( M>N_RowofMesh-2 )then
			!Print*,'Particle over the uper or lower edge'
			!Print*,'M and i is ', M, i
			M = N_RowofMesh-2
			Pt(i)%vx = Pt(i)%vx
			Pt(i)%vy = -Pt(i)%vy
		End If
		If( M<1 )then
			!Print*,'Particle over the uper or lower edge'
			!Print*,'M and i is ', M, i
			M = 1
			Pt(i)%vx = Pt(i)%vx
			Pt(i)%vy = -Pt(i)%vy
		End If
		Pt(i)%M = M
		Pt(i)%N = N
		SandInMesh( i,1 ) = M
		SandInMesh( i,2 ) = N
		j = Mesh(m,n,1)
		Mesh(m,n,j+2) = Pt(i)%i
		Mesh(m,n,1) = Mesh(m,n,1) + 1
	End Do!Loop Update Mesh and SandInMesh
	Mesh(0,:,1) = 1!Question Mesh(0,:,2) = N_Total + 1 Mesh(N_RowOfMesh-1,:,1) = 1
	Mesh(N_RowOfMesh-1,:,2) = N_Total + 1
	Mesh(:,0,:) = Mesh(:,N_ColumnOfMesh-2,:)
	Mesh(:,N_ColumnOfMesh-1,:) = Mesh(:,1,:)

End Subroutine UpdateParticleInformation


!//Argument of particle electric calculation
Subroutine ArgumentofElectric()
	Implicit None
	Real( Kind=8 ), Parameter ::		Zeta0=3.2E-8, Eps0=8.85E-12 ,Est=1013.25, Tst=373.15, Kb=1.38E-23
	Real( Kind=8 ), Parameter ::		Mol=6.02E23, H=6.63E-34
	Real( Kind=8 ) ::  					Lamda3
	Real( Kind=8 ) ::  					Ps1, Ps, P, T, Xi0, Rh, X, M, Pd11, Pd12, Pd13
   
	!Write(*,"(A40)")'Please input the Temperature '
	!Read( *,* ) Temperature
	Temperature = 25
	T = 273.16 + Temperature
	Call Integral(Pai,T)
	Xi0 = -33.0*1000
	Rh = 20
	Write(*,"(A40)")'Please input the Rh'
	Read( *,* ) Rh
	M = 1.67E-27
	Lamda3 = (2.0*Pi*M*Kb*T/H**2)**(-1.5)
	Ps1 = -7.90298*(Tst/T-1.0) + 5.02808*Log10(Tst/T) - 1.3816*1.0E-7*(10**(11.344*(1.0-T/Tst))-1.0) + &
			8.1328*1.0E-3*(10**(-3.49149*(Tst/T-1.0))-1.0) + Log10(Est)
	Ps = 100*10.0**(Ps1)
	P = Rh*Ps/100.0
	Pd11 = (Kb*T/Lamda3)
	Pd12 = (Pai/P)
	Pd13 = Exp(Xi0/(Kb*T*Mol))
	Pd1 = 1.0/(Pd11*Pd12*Pd13 + 1)
	Write(*,"(A40,3(F13.6,2x))")'Pd Pai Rh', Pd1, Pai, Rh

End Subroutine ArgumentofElectric

!//Collision elecrticity calculation
Real Function DeltaQ( ii,jj,Requiv )
	Implicit None
	Integer( kind=MyIntegerK ) :: 		ii, jj
	Real( Kind=MyRealk ) ::  			A1, A2, Requiv, Ai, Aj
	Real( Kind=MyRealk ), Parameter ::  N=1.0E18, Q=1.6E-19

	If( 2.0*Requiv*ComMax(ii,jj)/Pt(ii)%Radius**2>1.0 )then
		A1 = 2*Pi*Pt(ii)%Radius**2
		Print*,'E maybe too big'
	Else
		A1 = 2*Pi*Pt(ii)%Radius**2*(1.0-Sqrt(1.0 - 2.0*Requiv*ComMax(ii,jj)/Pt(ii)%Radius**2))
	End If
	If( 2.0*Requiv*ComMax(ii,jj)/Pt(jj)%Radius**2>1.0 )then
		A2 = 2*Pi*Pt(jj)%Radius**2
		Print*,'E maybe too big'
	Else
		A2 = 2*Pi*Pt(jj)%Radius**2*(1.0-Sqrt(1.0 - 2.0*Requiv*ComMax(ii,jj)/Pt(jj)%Radius**2))
	End If
	Ai = 4.0*Pi*Pt(ii)%Radius**2
	Aj = 4.0*Pi*Pt(jj)%Radius**2
	DeltaQ = 0.02*( Pt(jj)%Alpha*A2*Pt(jj)%Pd*(1.0-Pt(ii)%Pd) - Pt(ii)%Alpha*A1*Pt(ii)%Pd*(1.0-Pt(jj)%Pd) )
	Pt(ii)%Co = Pt(ii)%Co + 1
	Pt(jj)%Co = Pt(jj)%Co + 1
	Pt(ii)%Alpha = 0.132*0.8**Pt(ii)%Co
	Pt(jj)%Alpha = 0.132*0.8**Pt(jj)%Co
	Pt(ii)%Pd = Pt(ii)%Pd + Pt(ii)%Alpha*( A2*Pt(jj)%Pd*(1.0-Pt(ii)%Pd) - &
										   A1*Pt(ii)%Pd*(1.0-Pt(jj)%Pd) )/Ai
	Pt(jj)%Pd = Pt(jj)%Pd - Pt(ii)%Alpha*( A2*Pt(jj)%Pd*(1.0-Pt(ii)%Pd) - &
										   A1*Pt(ii)%Pd*(1.0-Pt(jj)%Pd) )/Aj

End Function DeltaQ



Subroutine ShowInformation( Pipe )
    Implicit None
    Character :: 	Pipe*15
	Integer( kind=MyIntegerK ) :: 		i
    
    Pipe = trim(adjustl(Pipe))
    Select Case( Pipe )
    Case( 'ReadVariable' )
        Print*
	    Print*,'**********This program claculate collision of sand and sand bed**********'
	    Print*
	    Write(*,*) 'Input time step : DeltaT '
	    Read(*,*) DeltaT
        Print*,'Input lable of module（Spring-1 Hertz-2 ElasticPlastic-3）'
	    Read(*,*) Mx    
        Print*,'Input number of incident sand: N_Incident'
	    Read(*,*) N_Incident
	    Allocate( Incident_R(N_Incident), Incident_V(N_Incident), Incident_A(N_Incident))
        Print*,'Input radius of incident particles((1-5)*E-4) Incident_R'
	    Read(*,*) (Incident_R(i), i=1,N_Incident)
	    Print*,'Input velocity of incident particles: Incident_V'
	    Read(*,*) (Incident_V(i), i=1,N_Incident)	
	    Print*,'Input initial angel of incident particles Incident_A'
	    Read(*,*) (Incident_A(i), i=1,N_Incident)	
    Case( 'ShowResult' )
        Print*
	    Print*,'Finish calculation......'
		Open( 12,File='Variable.txt' )
        Do i = 1,N_Total
            Write(12,"(1x,i5,1x,d20.10,3x,d20.10,3x,d20.10)") i, Pt(i)%Vx, Pt(i)%Vy!, Pt(i)%W  
        End Do 
	    Print*
	    Write(12, "('Step length ',E13.6)") DeltaT
        Write(12, "('Number of incident particles ',i8)")  N_Incident
        Write(12, "('Radius volicity Angle of incident particles ')")
        Write(12, "(8x, 3(F20.10,4x)/)")  (Incident_R(i), Incident_V(i), Incident_A(i), i=1, N_Incident)
		Close( 12 )
    Case( 'Calculation' )
	    Print*
        Print*,'Start calculation......'
    End Select
    
End Subroutine ShowInformation

Subroutine SaveResult( Bh )
	Implicit None
	Integer( kind=MyIntegerK ) :: 			i, Bh, j
	Character :: 							FileName*20

	write( FileName,* ) Bh
	!FileName = Trim(Adjustl(FileName))//'.bin'
	!Open( Bh,File=FileName,Access='Direct',Form='Unformatted',Recl=Int(MyIntegerK/4)+3*Int(MyRealk/4) )
	FileName = Trim(Adjustl(FileName))//'.txt'
	Open( Bh,File=FileName )
	Do i = 1, N_Total
		!write(Bh,Rec=i ) Pt(i)%i, Pt(i)%X, Pt(i)%Y, Pt(i)%Radius
		write(Bh,"(i9,2x,6(E13.6,3x))" ) Pt(i)%i, Pt(i)%X, Pt(i)%Y, Pt(i)%Radius, Pt(i)%Vx, Pt(i)%Vy, Pt(i)%Q
		!Q_Height( Int(Pt(i)%Y/CellOfQHeight)+1 ) = Q_Height( Int(Pt(i)%Y/CellOfQHeight)+1 ) + Pt(i)%Q
	End Do 
	close( Bh )
!	Open( 13,File='Q_Height.txt',Access='Append' )
!	!Open( 13,File='Q_Height.bin',Access='Direct',Form='Unformatted',Recl=4+8+8 )
!	Do i = 1, N_CellofQHeight
!		!Write( 13,Rec=i ) Time, CellOfQHeight*Real(i), Q_Height(i)
!		Write( 13,"(i9,3x,F8.3,3x,E13.6)" ) Time, CellOfQHeight*Real(i), Q_Height(i)
!	End Do 
!	Q_Height(:) = 0.0
!	close( 13 )
	Open( 1234,File='ElectricField',Position='Append')
	Call GetElectricField()
	ElectricField = ElectricField / Length
	Do i = 1, Size(ElectricField)
		Write(1234,'(2(E13.6,2x))') Real(i)*CellOfElectricField,  ElectricField(i)
	End Do 
	close( 1234 )
	Open( 2135,File='Qm',Position='Append')
	forall (i=1:Size(Qm,1))
		Qm(i) = Qm(i)/Real(Step)/DeltaT/Length/0.1
	End forall
	Do i = 1, Size(Qm,1)
		Write(2135,'(2(E13.6,2x))') Real(i)*Delta_Z_Acc, Qm(i)
	End Do 
	Qm = 0
	close( 2135 )

End Subroutine SaveResult

Subroutine Orbit()
	Implicit None
	Integer( kind=MyIntegerK ) :: 	i, i1, j, k, H, L, M, N, Zi, Febot, Fetop
	Real(Kind = Myrealk) :: 		CosAngle, CosAngleT, SinAngle, SinAngleT, ViNt, VjNt, Nnmod, Ntmod
	Real(Kind = Myrealk) :: 		Vi, Vj, Vt, Dl, Vwind, V_Relat, Re, Cd, Requiv, E1, E2, Ubot, Utop, Fl, Fe

	Do i = 1, N_Total!Loop all particle to calculation movement and Velocity
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
		F_DragX = 0.0
		F_DragZ = 0.0
		M = SandInMesh(i,1)
		N = SandInMesh(i,2)
		Do L = N-1, N+1!Loop Coloum of i particle's neighbor mesh grid
			Do H = M-1, M+1!Loop Row of i particle's neighbor mesh grid 	
				k = Mesh( H,L,1 )!Particle number in the mesh (H,L)
				If( k<1 )then
					Cycle
				End If
				!Calculation of compression of i and j
				Do i1 = 2, k+1!Loop particles in neighbor mesh grid 
					j = Mesh( H,L,i1 )
					If( j==i )then
						Cycle
					End If
					If(H==0)then
						If( Pt_Old(i)%Radius>Pt_old(i)%Y )then
							Requiv = Pt(i)%Radius
							E1 = Eee( Pt(i)%OrderofRadius ) * 1.0E6
							E2 = Eee( Pt(j)%OrderofRadius ) * 1.0E6
							E = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
							KnEdge = 4.0*E*Requiv**0.5/3.0
							Fn = KnEdge * ( Pt_Old(i)%Radius-Pt_Old(i)%Y )**1.5
						Else
							Fn = 0.0
						End If
						Fny = Fny + Fn
						Fnx = Fnx 
						Ftx = Ftx
						Fty = Fty
						Cycle
					Elseif(H==N_RowOfMesh-1)then
						If( Height<Pt_old(i)%Y+Pt_old(i)%Radius )then
							Requiv = Pt(i)%Radius
							E1 = Eee( Pt(i)%OrderofRadius ) * 1.0E6
							E2 = Eee( Pt(j)%OrderofRadius ) * 1.0E6
							E = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
							KnEdge = 4.0*E*Requiv**0.5/3.0
							Fn = -KnEdge * ( Pt_old(i)%Y+Pt_old(i)%Radius-Height )
						Else
							Fn = 0.0
						End If
						Fny = Fny + Fn
						Fnx = Fnx 
						Ftx = Ftx
						Fty = Fty
						Cycle
					Else
						If( L==0 )then
							Pt_old(j)%X = Pt_old(j)%X - Length
						Elseif( L==N_ColumnOfMesh-1 )then
							Pt_old(j)%X = Pt_old(j)%X + Length
						End If	
					End If
					Compression = Pt_Old(i)%Radius+Pt_Old(j)%Radius - ( (Pt_Old(i)%X-Pt_Old(j)%X)**2 + (Pt_Old(i)%Y-Pt_Old(j)%Y)**2 )**0.5
					If( Compression>0.0 )then!Particle i and j contact and Compressioned 
						!Calculation of Normal Force
						Requiv = Pt(i)%Radius * Pt(j)%Radius / ( Pt(i)%Radius + Pt(j)%Radius )
						E1 = Eee( Pt(i)%OrderofRadius ) * 1.0E6
						E2 = Eee( Pt(j)%OrderofRadius ) * 1.0E6
!						If( Pt(i)%Y>0.8*Gao )then
!							E1 = E1*10
!						End If
!						If( Pt(j)%Y>0.8*Gao )then
!							E2 = E2*10
!						End If
						E = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
						Kn = 4.0*E*Requiv**0.5/3.0
						Kt = Kn
						Fn = Kn * Compression**1.5
						!Fn = Kn * Compression
						Nnmod = ( (Pt_Old(i)%X-Pt_Old(j)%X)**2 + (Pt_Old(i)%Y-Pt_Old(j)%Y)**2 )**0.5
						CosAngle = ( Pt_Old(i)%X-Pt_Old(j)%X ) / Nnmod
						SinAngle = ( Pt_Old(i)%Y-Pt_Old(j)%Y ) / Nnmod	
						Fnx = Fn * CosAngle + Fnx
						Fny = Fn * SinAngle + Fny
						!Calculation of Tangential Force
						!Tangential vector  Nt = ( Pt_Old(j)%Y-Pt_Old(i)%Y, Pt_Old(i)%x-Pt_Old(j)%x )
						!Velocity vector Vi = ( Pt_Old(i)%Vx, Pt_Old(i)%Vy ) Vj = ( Pt_Old(j)%Vx, Pt_Old(j)%Vy )
						!Decompostion vector(Vi Vj) to tangent dircetion Nt
						Ntmod = Nnmod
						Vi = ( Pt_Old(i)%Vx**2 + Pt_Old(i)%Vy**2 )**0.5
						Vj = ( Pt_Old(j)%Vx**2 + Pt_Old(j)%Vy**2 )**0.5
						ViNt = ( Pt_Old(i)%Vx*(Pt_Old(j)%Y-Pt_Old(i)%Y) + Pt_Old(i)%Vy*(Pt_Old(i)%X-Pt_Old(j)%X) ) / Ntmod
						VjNt = ( Pt_Old(j)%Vx*(Pt_Old(j)%Y-Pt_Old(i)%Y) + Pt_Old(j)%Vy*(Pt_Old(i)%X-Pt_Old(j)%X) ) / Ntmod
						!ViNt = ViNt + Pt(i)%W*Pt(i)%Radius
						!VjNt = VjNt + Pt(j)%W*Pt(j)%Radius
						Vt = VjNt - ViNt
						If( Vt/=0.0 )then
							Ft = Miu * Fn * Vt/Abs(Vt)  !Question
						Else
							Ft = 0.0
						End If
						CosAngleT = (Pt_Old(j)%Y-Pt_Old(i)%Y) / Ntmod
						SinAngleT = (Pt_Old(i)%X-Pt_Old(j)%X) / Ntmod
						Ftx = Ft * CosAngleT + Ftx
						Fty = Ft * SinAngleT + Fty
						!Calculatin of damping force( Tangential and normal direction )
						Bn = ( 4.0*Kn*Pt_Old(i)%Mass*log(Epsl)*log(Epsl)/(Pi*Pi-log(Epsl)*log(Epsl)) )**0.5
						Bt = ( 4.0*Kt*Pt_Old(i)%Mass )**0.5
						Dn = -Bn * Compression / DeltaT 
						Dt = -Bt * Vt!Question Dt is not this form
						Dnx = Dn * CosAngle + Dnx
						Dny = Dn * SinAngle + Dny
						Dtx = Dt * CosAngleT + Dtx
						Dty = Dt * SinAngleT + Dty
						If( Vt/=0.0 )then
							Fw = Fw + Ft !+ Dt
						Else
							Fw = Fw 
						End If
						!Calculatin of electric between i and j
						If( Compression>ComMax(i,j) )then
							ComMax( i,j ) = Compression
						Else 
							Pt(i)%Q = Pt(i)%Q + DeltaQ( i,j,Requiv )
							ComMax( i,j ) = 0.0
						End If
					Else
						ComMax( i,j ) = 0.0
					End If!Particle i and j contact and Compressioned 
					If( L==0 )then
						Pt_old(j)%X = Pt_old(j)%X + Length
					Elseif( L==N_ColumnOfMesh-1 )then
						Pt_old(j)%X = Pt_old(j)%X - Length
					End If	
				End Do!Loop particles in neighbor mesh grid 
			End Do !Loop column of i particle's neighbor mesh grid
		End Do!Loop Row of i particle's neighbor mesh grid
		!Calculation of drag force
		If( Pt_Old(i)%Y-Z0>0.8*Gao )then
			Vwind = Ustar*Log((Pt_Old(i)%Y-GaoW)/Z0)/Kaman
			V_Relat = ( (Pt_old(i)%Vx-Vwind)**2 + (Pt_Old(i)%Vy)**2 )**0.5
			Re = V_Relat*Pt(i)%Radius*2.0/Ga
			Cd = 24.0/Re + 6.0/(1.0+(Re)**0.5) + 0.4
			F_DragX = -(Pi/8.0)*(4.0*Pt(i)%Radius**2)*Rho_Air*Cd*(Pt_Old(i)%Vx-Vwind)*V_Relat
			F_DragZ = -(Pi/8.0)*(4.0*Pt(i)%Radius**2)*Rho_Air*Cd*Pt_Old(i)%Vy*V_Relat
			Febot = Int( (Pt(i)%Y-Gao)/CellOfElectricField ) 
			If( Febot<1 ) Febot = 1
			Fetop = Febot + 1
			If( Fetop>Size( ElectricField,1 ) )Fetop = Febot
			Fe = Pt(i)%Q * 0.5*(ElectricField(Fetop)+ElectricField(Febot))
		End If
		!Calculation of movement
		Fx = Ftx + Fnx + F_DragX! + Dnx + Dtx
		Fy = Fty + Fny + F_DragZ! + Fe !+ Dny + Dty
		Pt(i)%Vx = Pt(i)%Vx + Fx/Pt(i)%Mass * Deltat
		Pt(i)%Vy = Pt(i)%Vy + Fy/Pt(i)%Mass * Deltat + G*Deltat
		Pt(i)%V = ( Pt(i)%vx**2+Pt(i)%vy**2 )**0.5
		Pt(i)%X = Pt(i)%X + Pt(i)%Vx * DeltaT + 0.5*Fx/Pt(i)%Mass*DeltaT**2
		Pt(i)%Y = Pt(i)%Y + Pt(i)%Vy * DeltaT + 0.5*( Fy/Pt(i)%Mass+G )*DeltaT**2
		!Pt(i)%W = Pt(i)%W + Fw*Pt(i)%Radius*DeltaT/Pt(i)%Inertia
		If( Pt(i)%Y<Pt(i)%Radius )then
			Pt(i)%Y = Pt(i)%Radius
			Pt(i)%Vx = 0.9*Pt(i)%Vx
			Pt(i)%Vy = 0.9*Pt(i)%Vy
		End If
		!Calculation of TauP
		If( Pt(i)%Y>Gao )then
			Call Get_TauP( i )
			zi = Int( (Pt(i)%Y-Gao)/Delta_Z_Acc ) + 1
			!Qm( zi ) = Qm( zi ) + Pt(i)%Q/Pt(i)%Mass
			Qm( zi ) = Qm( zi ) + Pt(i)%Mass
		End If

	End Do !Loop all particle to calculation movement and Velocity

End Subroutine Orbit


Subroutine GetElectricField()
	Implicit None
	Integer :: 						i, j, k, l, Row, Col, N
	Real( kind=MyRealk ):: 			Hei, Edis, Etmp

 	ElectricField(:) = 0.0
	Do i = 1, N_Total
		If( Pt(i)%Y<GaoE ) Cycle
		Hei = Height_ElectricField
		k = 0
		Do while( .True. ) 
			Hei = Hei + CellOfElectricField
			K = k + 1
			If( K>Size(ElectricField,1) ) Exit
			If( Pt(i)%Y==Hei ) Hei = Hei + 1.0E-20
			Edis = (Pt(i)%x-0.5*Length)**2 + (Pt(i)%Y-Hei)**2 
			If( Edis<5E-5 ) Cycle
			Etmp = -(Pt(i)%Y-Hei)/Abs(Pt(i)%Y-Hei)*Abs(Pt(i)%Y-Hei)/Sqrt(Edis)*9.0E9*Pt(i)%Q/Edis
			!Etmp = 9.0E9*Pt(i)%Q/Edis
			ElectricField(k) = ElectricField(k) + Etmp
		End Do 
	End Do 

End Subroutine GetElectricField


Subroutine ModifyWind()
    Implicit None
    Real( Kind=MyrealK ) ::                     Z_Surface_Temp(5000), Delta_Z_Surface, Z_Surface_Top, Z_Surface_Down, TauA_Surface, H, Pp
    Real( Kind=Myrealk ), Allocatable ::        Z_Surface(:), U_Surface(:), Dudz_Surface(:), TauP_Surface(:), Dudz(:)
    Integer( Kind=MyIntegerK ) ::               N, i, j, k, Remainder, K1, K2, j1, Ws

    !====================================计算空气剪切力===============================
	ws = N_Z+N_Z_Acc
	Call Interpolation()
	Allocate( Dudz(Ws) )
	TauA(:) = Tau
    Dudz(:) = 0.0
    U_New(:) = 0.0
    !====================================计算风速廓线的一阶导数 Dudz ==================
    Do i = 1, N_Z+N_Z_Acc
        TauA(i) = Tau - TauP(i)
        If( TauA(i)>0.0 )then
            Dudz(i) = ( TauA(i)/(Rho_Air*Kaman**2*Zwind(i)**2) )**0.5
        Else
			TauA(i) = 1.0E-10
            Dudz(i) = ( TauA(i)/(Rho_Air*Kaman**2*Zwind(i)**2) )**0.5
        End If
    End Do
    !====================================将 0-2*Delta_Z_Acc 这个高度内的网格再次细化 ===
    Z_Surface_Temp(:) = 0.0
    Delta_Z_Surface = Delta_Z_Acc
    j = 6
    k = 1
    Z_Surface_Temp(1) = Zwind(2)
	Do While( K<J )
		k = k + 1
		Delta_Z_Surface = Delta_Z_Surface/2.0
		Z_Surface_Temp(k) = Z_Surface_Temp(k-1) - Delta_Z_Surface
	End Do 
    Remainder = Int( (Z_Surface_Temp(j)-Z0)/Delta_Z_Surface )
    If( Remainder<2 )then
        Print*, '细化床面附近的网格时出错，请增大 Delta_Z_Surface 或者 J'
    End If
    Do i = j+1, j+Remainder
        Z_Surface_Temp(i) = Z_Surface_Temp(i-1) - Delta_Z_Surface
    End Do 
    If( Z_Surface_Temp(j+Remainder)<Z0 )then
        Print*, 'Z_Surface_Temp(细化网格的最低点低于Z0）'
    End If
    Allocate( Z_Surface(j+Remainder) )
    Z_Surface(:) = 0.0
    Do i = 1, j+Remainder
        Z_Surface(i) = Z_Surface_Temp(J+Remainder+1-i)
    End Do 
    !====================================以上代码所得网格点处的风速=====================
    Allocate( U_Surface(j+Remainder),Dudz_Surface(j+Remainder),TauP_Surface(j+Remainder) )
    U_Surface(:) = 0.0
    Dudz_Surface(:) = 0.0
    TauP_Surface(:) = 0.0
    Do i = 1, j+Remainder
        If( Z_Surface(i)<=Zwind(1) )then
            Z_Surface_Top = Zwind(1)
            Z_Surface_Down = Z0
            TauP_Surface(i) = ( (Z_Surface(i)-Z_Surface_Down)*TauP(1) + (Z_Surface_Top-Z_Surface(i))*Tau_Z0 ) / ( Z_Surface_Top-Z_Surface_Down ) 
        Else
            Z_Surface_Top = Zwind(2)
            Z_Surface_Down = Zwind(1)
            TauP_Surface(i) = ( (Z_Surface(i)-Z_Surface_Down)*TauP(2) + (Z_Surface_Top-Z_Surface(i))*TauP(1) ) / ( Z_Surface_Top-Z_Surface_Down ) 
        End IF
        TauA_Surface = Tau - TauP_Surface(i)
        Dudz_Surface(i) = TauA_Surface/Abs(TauA_Surface)*( Abs(TauA_Surface)/(Rho_Air*Kaman**2*Z_Surface(i)**2) )**0.5
    End Do 
    U_Surface(1) = Dudz_Surface(1)*( Z_Surface(1)-Z0 )
    U_Surface(2) = U_Surface(1) + 0.5*( Z_Surface(2)-Z_Surface(1) )*( Dudz_Surface(1)+Dudz_Surface(2) )
    U_Surface(3) = U_Surface(2) + (1.0/12.0)*( Z_Surface(3)-Z_Surface(2) )*( -Dudz_Surface(1)+8.0*Dudz_Surface(2)+5.0*Dudz_Surface(3) )
    Do i = 4, j+Remainder
        If( Z_Surface(i)-Z_Surface(i-1)==Delta_Z_Surface )then
            U_Surface(i) = U_Surface(i-1) + (1.0/24.0)*Delta_Z_Surface*( Dudz_Surface(i-3)-5.0*Dudz_Surface(i-2)+19.0*Dudz_Surface(i-1)+9.0*Dudz_Surface(i) )
        Else
            k = Int( Log( (Z_Surface(i)-Z_Surface(i-1))/Delta_Z_Surface )/Log(2.0) )
            U_Surface(i) = U_Surface(i-1) + (1.0/12.0)*( Z_Surface(i)-Z_Surface(i-1) )*( -Dudz_Surface(i-2-k)+8.0*Dudz_Surface(i-1)+5.0*Dudz_Surface(i) )
        End If
    End Do
    Do i = 1, j+Remainder
        If( Abs( Z_Surface(i)-Delta_Z_Acc )<0.1*Delta_Z_Surface )then
            k1 = i
        End If
        If( Abs( Z_Surface(i)-2*Delta_Z_Acc )<0.1*Delta_Z_Surface )then 
            k2 = i
        End If
    End Do 
    !====================================重新计算全局网格点处的风速=====================
    U_New(1) = U_Surface(k1)
    U_New(2) = U_Surface(k2)
    U_New(3) = U_New(2) + (Delta_Z_Acc/12.0)*( -Dudz(1)+8.0*Dudz(2)+5.0*Dudz(3) )
	H = Delta_Z_Acc
	Do i = 4, N_Z_Acc
        U_New(i) = U_New(i-1) + (H/24.0)*( Dudz(i-3)-5.0*Dudz(i-2)+19.0*Dudz(i-1)+9.0*Dudz(i) )
	End Do 
	H = Delta_Z
	i = N_Z_Acc+1
	j = Int( ( Zwind(i)-2*H )/Delta_Z_Acc )
	j1 = Int( ( Zwind(i)-H )/Delta_Z_Acc )
	U_New(i) = U_New(i-1) + (H/12.0)*( -Dudz(j)+8.0*Dudz(j1)+5.0*Dudz(i) )
	i = N_Z_Acc+2
	j = N_Z_Acc
	U_New(i) = U_New(i-1) + (H/12.0)*( -Dudz(j)+8.0*Dudz(i-1)+5.0*Dudz(i) )
	Do i = N_Z_Acc+3, N_Z_Acc+N_Z
        U_New(i) = U_New(i-1) + (H/24.0)*( Dudz(i-3)-5.0*Dudz(i-2)+19.0*Dudz(i-1)+9.0*Dudz(i) )
	End Do 
	Ustar = 0.0
	Do i = 1, ws
		Ustar = Ustar + U_New(i)*Kaman/log(Zwind(i)/Z0)
	End Do 
	Ustar = Ustar/Real(ws)
	Deallocate( Dudz,Z_Surface,U_Surface,Dudz_Surface,TauP_Surface )
	TauP(:) = 0.0
    
End Subroutine ModifyWind

Subroutine Get_TauP( i )
	Implicit None
	Integer( Kind=MyIntegerK ) ::  	i, Zi

	If( Pt(i)%Y-Gao>Delta_Z_Acc*N_Z_Acc )then
		Zi = 20 + Int( (Pt(i)%Y-Gao-Delta_Z_Acc*N_Z_Acc)/Delta_Z ) + 1
	Else
		Zi = Int( (Pt(i)%Y-Gao)/Delta_Z_Acc ) + 1
	End If
	TauP(zi) = TauP(zi) - Pt(i)%vy/Abs(Pt(i)%Vy)*Pt(i)%Mass*Pt(i)%Vx
	

End Subroutine Get_TauP

Subroutine Interpolation()
	Implicit None
	Integer( Kind=MyIntegerK ) :: 		i, j, k, M
	
	TauP(:) = TauP(:)/Edgeright
	If( TauP(1)==0.0 )then
		TauP(1) = 1.0E-10
	End IF
	Do i = N_Z_Acc+N_Z, 1, -1
		If( TauP(i)/=0.0 .and. i>1 )then
			Do j = i-1, 1, -1
				If( TauP(j)==0.0 )then
					Do k = j-1, 1, -1
						If( TauP(k)/=0.0 )then
							Do m = k+1, j
								TauP(m) = ( TauP(k)*(Zwind(i)-Zwind(m))+TauP(i)*(Zwind(m)-Zwind(k)) )/( Zwind(i)-Zwind(k) )
							End Do 
							Exit
						End If
					End Do 
				End If
			End Do 
		End If
	End Do 

End Subroutine Interpolation


End module Qiu_li
	
	

Program Main
    Use Qiu_li
    Implicit none	
	Integer( kind=MyIntegerK ) :: 		Flag, i, j
	Real( kind=MyRealk ) :: 			Start, Finish
	Character :: 						Pipe*15
	
!//---------------------------------------------------------Set Variables---------------------------------------------------------------//
	Call CPU_Time(Start)
	Call SetVariable()
	Call Get_Mesh()
	Call UpdateParticleInformation()
	Pt_old(:) = Pt(:)
!//---------------------------------------------------------Start Calculation---------------------------------------------------------------//
	Pipe = 'Calculation'
	Call ShowInformation( Pipe )
	j = 0
	Open( 789,File='log' )
	Do Flag = 1, N_Time, Step!Loop Step 
		Do Time = Flag, Flag+Step-1!Loop time in step 
			Pt_Old(:) = Pt(:)
			Call Orbit()
			Call UpdateParticleInformation()
		End Do!Loop Time in step
		Call ModifyWind()
		Call SaveResult(Bh)
		Write(*,"(A15,i10,A15,F13.6)") 'Time is:',Time-1, 'New Ustar is:',Ustar
		Write(789,"(A15,i10,A15,F13.6)") 'Time is:',Time-1, 'New Ustar is:',Ustar
		Bh = Bh + 1
	End Do!Loop step 
!//---------------------------------------------------------Result Process---------------------------------------------------------------//
	!Call DataProcess()
	Call CPU_Time(Finish)
	Print*,Finish-Start
	Read(*,*)
    
End program Main
