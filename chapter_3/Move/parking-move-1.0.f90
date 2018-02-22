Module Mmain
	Implicit None
	Integer,Parameter :: 					Myrealk=4, Myintk=4
	Integer :: 								N_Total, N_Sorts
	Integer(kind=Myintk),Allocatable ::		Neighbor(:,:), Neighborold(:,:), Mesh(:,:,:)
	Real( kind=Myrealk ) :: 				Br, Bt, Bl, Bb, Lij, MaxRadius, Top, Accuracy, MeshSize, MinRadius
	Real( kind=Myrealk ), Allocatable :: 	Radius(:), Ratio(:)
	Real( kind=Myrealk ),Parameter :: 		Pi=3.14159, G=9.8, Rho=2.65E3 
	Type :: Particle
		Integer ::							i, M, N
		Real( kind=Myrealk ) ::				X, Y, Radius, Mass
	End Type
	Type( Particle ), Allocatable :: 		Pt(:)

Contains


Subroutine SetVariable()
	Implicit None
	Integer ::						i, j, Err
	Real( kind=Myrealk ) ::			P
	Character ::					C, Conf

	Print*
	Print*,'If you want to use a configure file please confirm' 
	Print*,'it has been prepared in the current directory'
	Print*,'Input Y/y N/n'
	Read( *,* ) Conf
	If( Conf=='Y' .or. Conf=='y' )then
		Open( 11,File='Configure',Status='Old', Iostat=Err )
		If( Err/=0 )then
			Print*,'Can not find Configure file'
			Read(*,*)
		End If
		Read( 11,*, Iostat = Err ) Accuracy
		If( Err /= 0 )then
			Print*,'Please confirm Configure file is ok'
			Read(*,*)
			Stop
		End If
		Read( 11,* , Iostat = Err ) N_Total
		If( Err /= 0 .or.N_Total<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Total is:', N_Total
			Read(*,*)
			Stop
		End IF
		Allocate( Pt(N_Total) )
		Read( 11,*, Iostat = Err ) N_Sorts
		If( Err /= 0 .or. N_Sorts<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Sorts is:', N_Sorts
			Read(*,*)
			Stop
		End IF
		Read( 11,* ) Bb
		Read( 11,* ) Bt
		Read( 11,* ) Bl
		Read( 11,* ) Br
		Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
		Ratio(:) = 0.0
		Radius(:) = 0.0
		Read( 11,* ) Radius(:)
		Read( 11,* ) Ratio(:)
		Do i = 1, N_Total
			Read( 11,*,Iostat=Err ) Pt(i)%i,Pt(i)%Radius
			If( Err/=0 )then
				Print*,'Configure file maybe wrong'
				Print*,'N_Total of particle in the file is not enough'
				Read(*,*)
				Stop
			End If
		End Do 
		Close( 11 )
	Else
		Print*
		Print*,'Start.....'
		Print*
		Print*,'Input the Accuracy'
		Do 
			Read(*,*, Iostat = Err) Accuracy
			If( Err /= 0 .or. Accuracy>1.0 )then
				Print*,'Accuracy ', Accuracy
				cycle
			Else
				Exit
			End If
		End Do 
		Print*,'Input the total N_Total of sandbed'
		Do 
			Read(*,*, Iostat = Err) N_Total
			If( Err /= 0 )then
				Print*,'Input must be integer'
				cycle
			Else
				Exit
			End If
		End Do 
		Allocate( Pt(N_Total) )
		print*,'Input the sorts of sand: '
		Do 
			Read(*,*, Iostat = Err) N_Sorts
			If( Err /= 0 )then
				Print*,'Input must be Integer'
				cycle
			Else
				Exit
			End IF
		End Do 
		Allocate( Radius(N_Sorts),Ratio(N_Sorts) )
		Radius(:) = 0.0
		Ratio(:) = 0.0
		Print*,'Input the radius of sands(Unit/M)'
		Do i = 1, N_Sorts
			Print*,'Radius',i
			Do 
				Read(*,*,Iostat=Err) Radius(i)
				If( Err /= 0 )then
					Print*,'Try again'
					Cycle
				Else
					Exit
				End If
			End Do 
		End Do
		Bb = 0.0
		Bl = 0.0
		Print*,'Input the edge of sandbed--Br,Bt'
		Read(*,*) Br
		Read(*,*) Bt
		Call Get_Ratio(N_Sorts)
		Call Random_Radius()
	End If 
	MaxRadius = Maxval( Radius )
	MinRadius = Minval( Radius )

End Subroutine SetVariable


Subroutine Get_Ratio(N)
    Implicit None
    Integer :: N, i
    Real( kind=Myrealk ) :: Start, Finish, JianGe, S, Step
    
    Start = -3.0
    Finish = 3.0
    JianGe = (Finish-Start)/Float(N)
    S = Start
    i = 1
    Step = 1.0E-5
    Ratio(:) = 0.0
    Do 
        Ratio(i) = Ratio(i) + Step*1.0/Sqrt(2.0*3.1415)*Exp(-(S)**2*0.5)
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

End Subroutine Get_Ratio



Subroutine Random_Radius()
    Implicit None
    Integer ::                      N(N_Sorts), j, i
    Real( kind=MyrealK ) ::         RadiusList(N_Total), Sj
	Character ::                    a*3, aa*6

	N(:) = 0
    RadiusList(:) = 0.0
    Call RANDOM_SEED
    Do i = 1, N_Sorts    
        N(i) = Int(Real(N_Total)*Ratio(i))    
        If( N(i)==0 )then    
            N(i) = 1    
        End IF    
    End Do     
    N(Int(float(N_Sorts)*0.5)) = N(Int(float(N_Sorts)*0.5)) + N_Total - Sum(N)    
	Do i = 1, N_Sorts	
		Write(a,'(g0)') i	
		Write(aa,*) "R0",Trim(a),"="	
		Write(*,"(a,i8)") aa, N(i)	
	End Do	
    If ( N(Int(float(N_Sorts)*0.5)) < 0 )then    
        Print*,'Particle number wrong'
        Read(*,*)    
    End If    
    j = 1    
    Do i = 1, N_Sorts    
        RadiusList(j:N(i)+j-1) = Radius(i)    
        j = j + N(i)    
    End Do    
    Call UpsetList(RadiusList)          
    Do i = 1, N_Total    
        Pt(i)%i = i    
        Pt(i)%Radius = RadiusList(i)    
    End Do     
ENd Subroutine Random_Radius



Subroutine UpsetList( A )
    Implicit None
    Real( kind=Myrealk ) ::     Temp, R
    Real( kind=Myrealk ) ::     A(:)
    Integer ::                  i, j, Narray
    
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




Subroutine Geti( i )
	Implicit None
	Integer ::						i, j
	Real( kind=Myrealk ) ::			P,Lij

	Pt(i)%i = i
	Call Random_Number(P)
	Pt(i)%x = Pt(i)%Radius + P*(Br-2.0*Pt(i)%Radius)
	Pt(i)%y = Top+2.0*MaxRadius
	If( Pt(i)%Y-Pt(i)%Radius>Bt )then
		Print*, 'Fill number is ',i-1
		Call SaveConfigure()
		Call SaveResult( 11,i-1 )
		Stop
	End If

End Subroutine Geti



Subroutine Get_Mesh()
	Implicit None
	Integer( kind=Myintk ) ::	 			i,j,N_ColumnOfMesh,N_RowOfMesh,M,N

	MeshSize = 2.0*MaxRadius
	N_ColumnOfMesh = Int( Br/MeshSize )
	MeshSize = Br/Real(N_ColumnOfMesh)
	N_RowOfMesh = Int( Bt/MeshSize )
	i = Int( MeshSize**2/Pi/MinRadius**2 ) + 1
	Allocate( Mesh(0:i,0:N_ColumnOfMesh+1,0:N_RowOfMesh+1) )
	Mesh(:,:,:) = 0
	Mesh(1,0,:) = N_Total+1
	Mesh(0,0,:) = 1
	Mesh(1,N_ColumnOfMesh+1,:) = N_Total+3
	Mesh(0,N_ColumnOfMesh+1,:) = 1
	Mesh(1,:,0) = N_Total+2
	Mesh(0,:,0) = 1
	Mesh(1,:,N_RowOfMesh+1) = N_Total+4
	Mesh(0,:,N_RowOfMesh+1) = 1

	Print*,'N_ColumnOfMesh, N_RowOfMesh', N_ColumnOfMesh, N_RowOfMesh, i


End Subroutine Get_Mesh


Subroutine SaveResult(k,N)
	Implicit None
	Integer ::						k, i, N
	Character ::					FileName*30

	Write( FileName,* ) K
	FileName = Trim(Adjustl(FileName))//'.txt'
	Open( k,File=FileName )
	Do i = 1, N
		Write( k,'(i13,2x,3(E13.6,2x))' )  i,Pt(i)%X,Pt(i)%Y,Pt(i)%Radius
	End Do 
	Close( k )

End Subroutine SaveResult


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
		If( Pt(i)%X<Pt(i)%Radius )then
			Pt(i)%X = Pt(i)%Radius
			Njc = N_Total + 1
			Exit
		End If
		M = Floor( Pt(i)%Y/MeshSize ) + 1
		N = Floor( Pt(i)%X/MeshSize ) + 1
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


Subroutine Right( i,Njc )
	Implicit None
	Integer ::							i, j, Jc, Njc
	Real( kind=Myrealk ) ::				Lij, Xstep

	Xstep = 1.0E-3*Pt(i)%Radius
	Jc = 0
	Njc = 0
	Do while ( Jc<1 )
		Pt(i)%X = Pt(i)%X + Xstep
		If( Pt(i)%X+Pt(i)%Radius>Br )then
			Pt(i)%X = Pt(i)%Radius
			Njc = N_Total + 1
			Exit
		End If
		Do j = 1, i!N_Total
			If(j==i) cycle
			Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
			If( Pt(i)%Radius+Pt(j)%Radius>Lij )then
				Jc = Jc + 1
				Njc = j
				Exit 
			End If
		End Do 
	End Do 
	
End Subroutine Right



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
		If( Pt(i)%Y<Pt(i)%Radius )then
			Pt(i)%Y = Pt(i)%Radius
			Njc = N_Total + 2
			Exit
		End If
		M = Floor( Pt(i)%Y/MeshSize ) + 1
		N = Floor( Pt(i)%X/MeshSize ) + 1
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



Subroutine Rotate( i,Tstep,Njc )
	Implicit None
	Integer ::							i, j, k, Njc, l
	Real( kind=Myrealk ) ::				Lij, Tstep, Thetaij, Lik

	l = 0
	j = Njc
	Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
	Thetaij = Acos( (Pt(i)%X-Pt(j)%X)/Lij )
	Do while( .True. )
		If( Tstep>0.0 )then
			If( Pt(i)%X<Pt(j)%X ) Exit
		Else
			If( Pt(i)%X>Pt(j)%X ) Exit
		End If
		Thetaij = Thetaij - Tstep
		Pt(i)%x = Pt(j)%X + Lij*Cos(Thetaij)
		Pt(i)%y = Pt(j)%Y + Lij*sin(Thetaij)
		If( Pt(i)%x<Pt(i)%Radius )then
			Pt(i)%x = Pt(i)%Radius
			Njc = N_Total + 1
			Exit
		End If
		If( Pt(i)%x+Pt(i)%Radius>Br )then
			Pt(i)%x = Br - Pt(i)%Radius
			Njc = N_Total + 3
			Exit
		End If
		If( Pt(i)%Y<Pt(i)%Radius )then
			Pt(i)%Y = Pt(i)%Radius
			Njc = N_Total + 2
			Exit
		End If
		Do k = 1, i!N_Total
			If(k==i .or. k==j) cycle
			Lik = sqrt( (Pt(i)%X-Pt(k)%X)**2 + (Pt(i)%Y-Pt(k)%Y)**2 )
			If( Pt(i)%Radius+Pt(k)%Radius>Lik )then
				Njc = k
				j = Njc
				Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
				Thetaij = Acos( (Pt(i)%X-Pt(j)%X)/Lij )
				Exit 
			End If
		End Do 
	End Do 

End Subroutine Rotate


Subroutine Rotate1( i,Tstep,Njc )
	Implicit None
	Integer ::							i, j, k, Njc
	Integer ::							i1, i2, i3, l, M, N
	Real( kind=Myrealk ) ::				Lij, Tstep, Thetaij, Lik

	l = 0
	j = Njc
	Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
	Thetaij = Acos( (Pt(i)%X-Pt(j)%X)/Lij )
	If( Pt(i)%Y>Pt(j)%Y )then
		Thetaij = Thetaij
	Else
		Thetaij = 2.0*Pi - Thetaij 
	End If
	Do while( .True. )
		If( Tstep<0.0 )then
			If( Pt(i)%X<Pt(j)%X ) Exit
		Else
			If( Pt(i)%X>Pt(j)%X ) Exit
		End If
		Thetaij = Thetaij + Tstep
		Pt(i)%x = Pt(j)%X + Lij*Cos(Thetaij)
		Pt(i)%y = Pt(j)%Y + Lij*sin(Thetaij)
		If( Pt(i)%x<Pt(i)%Radius )then
			Pt(i)%x = Pt(i)%Radius
			Njc = N_Total + 1
			Exit
		End If
		If( Pt(i)%x+Pt(i)%Radius>Br )then
			Pt(i)%x = Br - Pt(i)%Radius
			Njc = N_Total + 3
			Exit
		End If
		If( Pt(i)%Y<Pt(i)%Radius )then
			Pt(i)%Y = Pt(i)%Radius
			Njc = N_Total + 2
			Exit
		End If
		M = Floor( Pt(i)%Y/MeshSize ) + 1
		N = Floor( Pt(i)%X/MeshSize ) + 1
		Loop1:Do i1 = M-1,M
			Do i2 = N-1,N+1
				l = Mesh(0,i2,i1)
				If( l<1 ) Cycle
				Do i3 = 1, l
					k = Mesh(i3,i2,i1)
					If(k==i .or. k==j .or. k>N_Total) cycle
					Lik = sqrt( (Pt(i)%X-Pt(k)%X)**2 + (Pt(i)%Y-Pt(k)%Y)**2 )
					If( Pt(i)%Radius+Pt(k)%Radius>Lik )then
						Njc = k
						j = k
						Lij = sqrt( (Pt(i)%X-Pt(j)%X)**2 + (Pt(i)%Y-Pt(j)%Y)**2 )
						Thetaij = Acos( (Pt(i)%X-Pt(j)%X)/Lij )
						Exit Loop1
					End If
				End Do 
			End Do 
		End DO Loop1
	End Do 

End Subroutine Rotate1


Subroutine ShowInformation()
	Implicit None
	Integer ::								i
	Character(len=100) ::					A, B, C(20)

	A = '//-------------------------------------------------------------------------------------//'
	C(1) = '// '
	C(2) = '// Particle parking simulation. '
	C(3) = '// '
	C(4) = '// Tightly stacked balls in a 2D rectangular. Similar to Y.T.Feng(2003)  '
	C(5) = '//	but they do not consider structural balance and compression of the '
	C(6) = '//	particles. This algorithm use a new method which simulated gravity '
	C(7) = '//	of the structure '
	C(8) = '// '
	C(9) = '// Author:	MU LiAn LanZhou University '
	C(10) = '// Date:	2017/07/06 '
	C(11) = '// Mail:	1791430496@qq.com '
	C(12) = '// Edition:	1.0 '
	Write(*,*)
	Write(*,*)
	Write(*,"(A100)") Adjustl(A)
	Do i = 1, 12
		B = C(i)
		Write(*,"(A100)") Adjustl(B)
	End Do 
	Write(*,"(A100)") Adjustl(A)
	Write(*,*)
	Write(*,*)

End Subroutine ShowInformation


Subroutine MeshUpdate(i)
	Implicit None
	Integer ::			i, M, N, j
	
	M = Floor( Pt(i)%Y/MeshSize ) + 1
	N = Floor( Pt(i)%X/MeshSize ) + 1
	j = Mesh(0,n,m)
	Mesh(j+1,n,m) = Pt(i)%i
	Mesh(0,n,m) = Mesh(0,n,m) + 1
	Pt(i)%M = M
	Pt(i)%N = N

End Subroutine MeshUpdate

Subroutine FillRate( Rate )
	Implicit None
	Integer ::					j, i
	Real( kind=Myrealk ) :: 	V1, V2, Xmax=0.0, Ymax=0.0, Rate

	i = Int( N_Total/2 )
    Do j = 1, i
		If( Pt(j)%X > Xmax )then
			Xmax = Pt(j)%X
		End If
		If( Pt(j)%Y > Ymax )then
			Ymax = Pt(j)%Y
		End If
		V1 = V1 + 3.14*Pt(j)%Radius * Pt(J)%Radius
    End Do
    V2 = Xmax*Ymax
	Rate = V1/V2

End Subroutine FillRate


Subroutine SaveConfigure()
	Implicit None
	Integer ::				j
	Character ::						WFormat*30

	Open( 11, File='Configure' )
	Write( 11,"(E10.3,2x,A30)" ) Accuracy, 'Newton ineration Accuracy'
	Write( 11,"(i10,2x,A30)" ) N_Total, 'Particle N_Total'
	Write( 11,"(i10,2x,A30)" ) N_Sorts, 'Particle sorts'
	Write( 11,"(E10.3,2x,A30)" ) Bb,'Bottom edge'
	Write( 11,'(E10.3,2x,A30)' ) Bt,'Top edge'
	Write( 11,'(E10.3,2x,A30)' ) Bl,'Left edge'
	Write( 11,'(E10.3,2x,A30)' ) Br,'Right edge'
	Write( WFormat,* ) N_Sorts
	WFormat = '('//Trim(Adjustl(WFormat))//'(E10.3),2x,A30)'
	Write( 11,WFormat ) Radius,'Radius'
	Write( 11,WFormat ) Ratio, 'Ratio'
    Do j = 1, N_Total
		Write(11, "(I6, 2(f13.6))") Pt(j)%i, Pt(j)%Radius
    End Do
	Close(11)


End Subroutine SaveConfigure


End Module Mmain






Program Main
	Use Mmain
	Implicit None
	Integer ::							i, j, k, l, Jc, Njc, Time
	Real( kind=Myrealk ) ::				Tstep, Finish, Start, Rate
	
	Call CPU_Time(Start)
	Call ShowInformation()
	Call SetVariable()
	Call Get_Mesh()

	Top = 2.0*MaxRadius
	Pt(1)%X = Pt(1)%Radius
	Pt(1)%Y = Pt(1)%Radius
	Call MeshUpdate(1)
	Do i = 2, N_Total
		Call Geti( i )
		Print*,'Move ',i
		Njc = 0
		J = -1
		! Move particle left and right
		Do while( .True. )
			Call Down( i,Njc ) !First
			Call Left( i,Njc ) !Second
			If( J==Njc ) Exit
			J = Njc
		End Do 
		! Rotate particle along particle Njc
		If( Njc>N_Total)then
			If( Top<Pt(i)%Y+Pt(i)%Radius ) Top = Pt(i)%Y+Pt(i)%Radius
			Call MeshUpdate(i)
			Cycle
		Elseif( Njc>0 )then
			j = Njc
			If( Pt(i)%X>=Pt(j)%X )then
				Do while( .True. )
					Tstep = -0.01
					Call Rotate1( i,Tstep,Njc )
					If( Njc>N_Total ) Exit
					If( Pt(i)%X<Pt(Njc)%x ) Exit
					If( l==Njc ) Exit
					l = j
					j = Njc
				End Do
			Else
				Do while( .True. )
					Tstep = 0.01
					Call Rotate1( i,Tstep,Njc )
					If( Njc>N_Total ) Exit
					If( Pt(i)%X>Pt(Njc)%x ) Exit
					If( l==Njc ) Exit
					l = j
					j = Njc
				End Do
			End If
			If( Top<Pt(i)%Y+Pt(i)%Radius ) Top = Pt(i)%Y+Pt(i)%Radius
			Call MeshUpdate(i)
		Else
			Print*,'Something maybe wrong'
		End If
	End Do 

	Call FillRate( Rate )
	Call CPU_Time(Finish)
	Call SaveConfigure()
	Call SaveResult( 11,N_Total )

	Print*,'Runing time is :', Finish-Start
	Print*,'Fill rate is :', Rate


End Program Main
