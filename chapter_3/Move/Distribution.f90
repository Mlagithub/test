Module Distribution
	Implicit None
	Integer,Parameter ::								Myrealk=4, Myintk=4

Contains


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



Subroutine LogNormal_Ratio( Mean,Sigma,N,Ratio )
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

End Subroutine LogNormal_Ratio



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



!Subroutine LogNormalDistribution( Mean,Sigma,N,Radius )
!	Implicit None
!	Integer ::								i, N
!	Real( kind=MyrealK ),Allocatable ::		Radius(:), P(:)
!	Real( Kind=Myrealk ) ::					Mean, Sigma
!
!End Subroutine LogNormalDistribution


End Module Distribution


!Program Main
!	Use Distribution
!	Implicit None
!	Integer ::										N_Total, N_Sorts, i
!	Real( kind=Myrealk ), Allocatable ::			Radius(:), Rad(:), Ratio(:)
!	Real( kind=Myrealk ) ::							Radii
!
!	Open( 11,File='Out.txt' )
!
!	N_Total = 10000
!	N_Sorts = 21
!	Allocate( Radius(N_Sorts),Rad(N_Total),Ratio(N_Sorts) )
!	Radii = 1.0
!	Do i = 1, N_Sorts
!		Radius(i) = Radii
!		Radii = Radii + 0.1
!	End Do 
!	!Radius(:)=(/1.0,2.0,3.0,4.0,5.0,6.0,7.0/)
!
!
!	!Call RandomDistribution( Nt=N_Total,Ns=N_Sorts,Rmax=1.0,Rmin=3.0,Rs=Radius,R=Rad,Ratio=Ratio )
!	Call NormalDistribution( Mean=10.0,Sigma=1.0,Nt=N_Total,Ns=N_Sorts,Rs=Radius,R=Rad,Ratio=Ratio )
!	Do i = 1, N_Total
!		Write( 11,"(E13.6)" ) Rad(i)
!	End DO 
!	Do i = 1, N_Sorts
!		Print*,Ratio(i)
!	End Do 
!
!End Program Main
