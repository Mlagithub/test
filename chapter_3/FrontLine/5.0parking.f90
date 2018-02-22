!//全局变量    
Module Quan_Ju_Bian_Lian
    Implicit None
    Integer,Parameter :: 					Myintk = Selected_Int_kind(8), Myrealk = Selected_Real_Kind(10,10)
    Integer :: 								i, Number, N_Sortofballs
	Integer(kind=Myintk),Allocatable ::		Mesh(:,:,:), SandInMesh(:,:), Neighbor(:,:)
	Real( kind=Myrealk ),Parameter :: 		Pi=3.14159
	Real(Kind = Myrealk) :: 				Br, Rinf , Bl, Bt, FillRate, Bb, Accuracy
    Real(kind = MyrealK), Allocatable :: 	Radius(:), Ratio(:)
    Type :: Particles
        Integer :: Xh
        Real(Kind = Myrealk) :: X, Y, Radius
    End Type Particles
    Type(Particles), Target, Allocatable :: Particle(:)
    Type :: Datalink
        Integer :: Bh
        Type(Particles) :: Item
        Type(Datalink), Pointer :: Prev
        Type(Datalink), Pointer :: Next
    End Type Datalink
End Module Quan_Ju_Bian_Lian
    
!//核心算法  
! Items 存放每次新生成的粒子的信息
! Head  存放所有粒子的信息
! ForwardLine 为活动指针，每次计算新粒子时应当指向当前活动的前进边的始端
Module Fl
    Use Quan_Ju_Bian_Lian
    Implicit None
	Integer :: Hs
    Type(Datalink), Pointer :: ForwardLine, Head, Items
    Real(kind=Myrealk) :: R_work(3), CC(4)
    
Contains

!//cc(1),cc(2)为已知的圆心坐标x1,y1
!//cc(3),cc(4)为已知的圆心坐标x2,y2
!//R_Work(1),R_Work(2),R_Work(3) 分别为三个圆形的半径  
Subroutine dsnse(i)
    Implicit none
    Integer :: l,j,i,Num,N
    Real(kind=Myrealk) :: d,s,yy(2),xx(2), f
    l = 5.0E5
    Num = 2
	N = 0
    XX(1) = Particle(i)%X
    XX(2) = Particle(i)%Y
    Do 
		If( N>5 )then
			Print*,'Newtonian iteration accuracy maybe too small'
			Exit
		End If
        Call fs(yy, xx, f)
        If (f >= Accuracy) Then
            l = l - 1
            If (l==0)then
				N = N + 1
				l = 5.0E5
			End If
            d = 0.0
            Do j = 1, num
                d = d + yy(j)*yy(j)
            End Do
            If (d==0.0) Then
                l = -1
				Print*,'Particle_Get ChuZhi maybe not suitable'
				Read(*,*)
                Return
            End If
            s = f/d
            Do j = 1, num
                xx(j) = xx(j) - s*yy(j)
            End Do
		else
            exit
        End If
    End do 
    Particle(i)%Xh = i
    Particle(i)%X = XX(1)
    Particle(i)%Y = XX(2)
    Return
End Subroutine dsnse 

Subroutine fs(yy, XX, f)
    Implicit none
    Real(kind=Myrealk) :: f, f1, f2, df1, df2, c(6), r1, r2, r3, yy(2), XX(2)
    r1=r_work(1)
    r2=r_work(2)
    r3=r_work(3)
    c(1)=cc(1)
    c(2)=cc(2)
    c(3)=c(1)**2 + c(2)**2 - (r1 + r3)**2
    c(4)=cc(3)
    c(5)=cc(4)
    c(6)=c(4)**2 + c(5)**2 - (r2 + r3)**2
    f1 = xx(1)**2 + xx(2)**2 - 2.0*c(1)*xx(1) - 2.0*c(2)*xx(2) + c(3)
    f2 = xx(1)**2 + xx(2)**2 - 2.0*c(4)*xx(1) - 2.0*c(5)*xx(2) + c(6)
    f = f1*f1 + f2*f2 
    df1 = 2.0*xx(1) - 2.0*c(1)
    df2 = 2.0*xx(1) - 2.0*c(4)
    yy(1) = 2.0*(f1*df1+f2*df2)
    df1 = 2.0*xx(2) - 2.0*c(2)
    df2 = 2.0*xx(2) - 2.0*c(5)
    yy(2) = 2.0*(f1*df1+f2*df2)
    Return
End Subroutine fs    

!//生成颗粒半径    
Subroutine Random_Radius()
    Implicit None
    Integer ::                      N(N_SortofBalls),j
    Real( kind=MyrealK ) ::         RadiusList(Number), Sj
	Character ::                    a*3, aa*6

	N(:) = 0
    RadiusList(:) = 0.0
    Call RANDOM_SEED
	Particle(1)%xh = 1	
	Particle(1)%Radius = Radius(Int(float(N_SortofBalls)*0.5))	
    Do i = 1, N_SortofBalls    
        N(i) = Int(float(Number)*Ratio(i))    
        If( N(i)==0 )then    
            N(i) = 1    
        End IF    
    End Do     
	Do i = 1, N_SortofBalls	
		Write(a,'(g0)') i	
		Write(aa,*) "R0",Trim(a),"="	
		Write(*,"(a,i8)") aa, N(i)	
	End Do	
    N(Int(float(N_SortofBalls)*0.5)) = N(Int(float(N_SortofBalls)*0.5)) + Number - Sum(N)    
    If ( N(Int(float(N_SortofBalls)*0.5)) < 0 )then    
        Print*,'粒子个数计算错误'    
        Read(*,*)    
    End If    
    j = 1    
    Do i = 1, N_SortofBalls    
        RadiusList(j:N(i)+j-1) = Radius(i)    
        j = j + N(i)    
    End Do    
    Call UpsetList(RadiusList)          
    Do i = 2, Number    
        Particle(i)%Xh = i    
        Particle(i)%Radius = RadiusList(i)    
    End Do     
ENd Subroutine Random_Radius
         

!//打乱数组的元素
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

!//得到各种粒子占比例
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

!//得到各种粒子占比例
Subroutine Get_RatioLog(N)
    Implicit None
    Integer :: N, i
    Real( kind=Myrealk ) :: Start, Finish, JianGe, S, Step, Sigma, M
    
    Start = 0.01
    Finish = 0.3
    JianGe = (Finish-Start)/Float(N)
    S = Start
    i = 1
    Step = 1.0E-5
    Ratio(:) = 0.0
	Sigma = 0.279
	m = -2.19
    Do 
        Ratio(i) = Ratio(i) + Step*1.0/(Sigma*s*Sqrt(2.0*3.1415))*Exp(-(log(s)-m)**2/(2.0*Sigma**2))
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


End Subroutine Get_RatioLog


!角点粒子生成（左下角点）
Subroutine CornerParticle
    Implicit None
    Integer :: i
	Real(Kind = Myrealk) :: Theat, R
	Theat = 90
	R = Radius(Int(float(N_SortofBalls)*0.5))
    Particle(1)%X = (R*Sin(theat/180.0*Pi))/(1.0+cos(theat/180.0*Pi))
    Particle(1)%Y = R
	Particle(1)%Radius = R
End Subroutine CornerParticle


!!初始化前进边
Subroutine Initialize_ForwardLine
    Implicit None
    Integer :: 		j, Err, i
	
	
	! Bh=-1 means left edge of the domains
    Allocate(Head)
    Head%Bh = -1
    Head%Item = Particles(-1, -Br, Particle(1)%Radius, Rinf )
    Nullify(Head%Next)
    Nullify(Head%Prev)
    ForwardLine => Head
	
	
	Allocate(ForwardLine%Next, Stat = Err)
	If(Err /= 0) Stop  
	ForwardLine%Next%Bh = 1
	ForwardLine%Next%Item = Particles(Particle(1)%Xh, Particle(1)%X, Particle(1)%Y, Particle(1)%Radius)
	ForwardLine%Next%Prev => ForwardLine
	ForwardLine => ForwardLine%Next
	ForwardLine%Prev%Next => ForwardLine
	Nullify(ForwardLine%Next)
	
	
	! Bh=-3 means bottom edge of the domains
	Allocate( ForwardLine%Next, Stat = Err )
	If(Err /= 0) Stop
	ForwardLine%Next%Bh = -3
	ForwardLine%Next%Item = Particles(-3, Particle(1)%X, -Rinf, Rinf)
	ForwardLine%Next%Prev => ForwardLine
	ForwardLine => ForwardLine%Next
	ForwardLine%Prev%Next => ForwardLine
	Nullify(ForwardLine%Next)
	
	
	! Bh=-2 means right edge of the domains
	Allocate( ForwardLine%Next, Stat = Err )
	If(Err /= 0) Stop
	ForwardLine%Next%Bh = -2
	ForwardLine%Next%Item = Particles(-2, Br + Rinf, Particle(1)%Y, Rinf)
	ForwardLine%Next%Prev => ForwardLine
	ForwardLine => ForwardLine%Next
	ForwardLine%Prev%Next => ForwardLine
	Nullify(ForwardLine%Next)
	
    ForwardLine => Head%Next

End Subroutine Initialize_ForwardLine


!粒子生成----生成第 i 个粒子
Subroutine Particle_Get(i)
    Implicit None
    Integer :: i, jj

    Allocate(Items)
	jj = ForwardLine%Bh
	! If pointer is at right edge, move it to next pointer(left edge)
	If( ForwardLine%Next%Bh == -2 )then
		ForwardLine =>Head
	End If
	If( ForwardLine%Bh == -1 )then
!		If( ForwardLine%Next%Bh==-3 )then
!			Particle(i)%X = Particle(i)%Radius
!			Particle(i)%Y = Particle(i)%Radius
!			Items%Bh = i
!			Items%Item = Particles(Particle(i)%Xh, Particle(i)%X, Particle(i)%Y, Particle(i)%Radius)
!			Return
!		Else
!			ForwardLine =>ForwardLine%Next
!		End If
			ForwardLine =>ForwardLine%Next
	End If
	! If next pointer is bottom edge
	If(ForwardLine%Next%Bh == -3)then
		Particle(i)%Y = Particle(i)%Radius
		Particle(i)%X = ForwardLine%Item%X + &
			& Sqrt((Particle(i)%Radius + ForwardLine%Item%Radius)**2 - (Particle(i)%Y - ForwardLine%Item%Y)**2)		
	! If next pointer is right edge
	Elseif(ForwardLine%Next%Bh == -2)then
		Particle(i)%X = Br - Particle(i)%Radius
		Particle(i)%Y = ForwardLine%Item%Y + &
			& Sqrt((Particle(i)%Radius + ForwardLine%Item%Radius)**2 - (Particle(i)%X - ForwardLine%Item%X)**2)
	! If pointer is left edge
	Elseif(ForwardLine%Bh == -1)then
		Particle(i)%X = Particle(i)%Radius
		Particle(i)%Y = ForwardLine%Next%Item%Y + &
			& Sqrt((Particle(i)%Radius + ForwardLine%Next%Item%Radius)**2 - &
			(Particle(i)%X - ForwardLine%Next%Item%X)**2)
	Else
		CC(1) = ForwardLine%Item%X
		CC(2) = ForwardLine%Item%Y
		CC(3) = ForwardLine%Next%Item%X
		CC(4) = ForwardLine%Next%Item%Y
		R_work(1) = ForwardLine%Item%Radius
		R_work(2) = ForwardLine%Next%Item%Radius
		R_work(3) = Particle(i)%Radius
		Call Chu_Zhi(i)
		Call dsnse(i)
	End IF									
    Items%Bh = i
    Items%Item = Particles(Particle(i)%Xh, Particle(i)%X, Particle(i)%Y, Particle(i)%Radius)

End Subroutine Particle_Get

!为计算新粒子的位置赋初值
Subroutine Chu_Zhi(i)
    Implicit None
    Integer :: i
    Real(Kind = Myrealk) :: Xj, Yj, R, Tmp
	R = Particle(i)%Radius
    Xj = CC(3) - CC(1)    
    Yj = CC(4) - CC(2)    
	If( Yj==0 )then
		If( Xj>0 )then
			Yj = 1.0E10
		Elseif( Xj<0 )then
			Yj = -1.0E10
		Else
			Print*,'Particle Get -- ChuZHi is wrong'
			Read(*,*)
			Stop
		End If
	End If
	If( Xj==0 )then
		Xj = 1.0E-20
	End If
	Tmp = Abs(Yj/Xj)
	If(Xj>0)then
		Particle(i)%y = 1*R + CC(2)
	Else
		Particle(i)%y = -1*R + CC(2)
	End If
	If(Yj>0)then
		Particle(i)%x = -1 * Tmp*R + CC(1)
	Else
		Particle(i)%x = Tmp*R + CC(1)
	End If
End Subroutine Chu_Zhi


!重叠检查
! 情况 1 -- 没有相重叠，生成的粒子合适
! 情况 2 -- 与现在位置之后的粒子 a 相重叠，删除现在位置到 a 之间的粒子（沿正向删除）（不包括 a 和现在位置的粒子）
! 情况 3 -- 与现在位置之前的粒子 b 相重叠，删除现在位置到 b 之间的粒子（沿负向删除）（不包括 b ，但包括现在位置的粒子）
! 情况 4 -- 与现在位置之前的粒（b）子和之后的粒子(a)都有重叠，删除 a b 之间的粒子（不包括 a b）
! 进入此程序时，ForwardLine 应当指向当前活动边的始端
! 之前与之后指链表从头开始，在前的为Q，在后的为H
Subroutine Check_Overlap(i, Qk, Hbh, Qbh, Bj)
    Implicit None
    Integer :: Qk, Hbh, Qbh, Bj, i, Bjq, Bjh, j
    Real(Kind = Myrealk) :: X0, Y0, R0, Dence, Rd, Rc
    Logical ::  AfterBj
	AfterBj = .False.
	Qk = 1
	Hbh = 0
	Qbh = 0
	Bjq = 0
	Bjh = 0
    Rc = 1.0E-7
    Bj = ForwardLine%Bh 
	If( Bj/=-1 )then
		Bjq = ForwardLine%Prev%Bh
	Else
		AfterBj = .True. 
	End IF
	If( Bj/=-2 )then
		Bjh = ForwardLine%Next%Bh
	End IF
	ForwardLine => Head
	If( Particle(i)%Y > 0.4*Bt )then
		Bl = Abs(3.0*Particle(i)%x-Particle(i)%Y)/Sqrt(10.0)
	End If
	If( Particle(i)%X - Particle(i)%Radius < Bl )then! Overlap with left edge
		Qbh = -1
		If( Particle(i)%X + Particle(i)%Radius > Br )then
			Print*,'Overlap with left and right edge, this maybe wrong'
			Print*,'Please enlarge Br'
			Read(*,*)
			stop
		End If
	Else
		If( Particle(i)%X + Particle(i)%Radius > Br )then
			Hbh = -2
		End If
	End If
	Do 
		If( .Not. Associated(ForwardLine%Next) )then
			Exit
		End If
		ForwardLine => ForwardLine%Next
		If( ForwardLine%Bh ==  Bj )then
			AfterBj = .True.
			Cycle
		End If
		If( ForwardLine%Bh==-1 .or. ForwardLine%Bh==-2  )then
			Cycle
		End If
		If( ForwardLine%Bh == Bjh )then
			Cycle
		End If
		X0 = ForwardLine%Item%X
		Y0 = ForwardLine%Item%Y
		R0 = ForwardLine%Item%Radius
		Dence = Sqrt( (X0 - Particle(i)%X)**2 + (Y0 - Particle(i)%Y)**2 )
		Rd = R0 + Particle(i)%Radius
		If( Abs(Dence - Rd) > Rc .And. Dence < Rd ) then
			If( AfterBj )then
				Hbh = ForwardLine%Bh
				Exit
			Else
				Qbh = ForwardLine%Bh
			End If
		End If
	End Do

    If( Qbh==0 )then
		If( Hbh==0 )then
			Qk = 1
		Else
			Qk = 2
		End If
	Else
		If( Hbh==0 )then
			Qk = 3
		Else
			Qk = 4
		End If
	End If
	If( Qk==1 )then
		Do j = 1, i-1
			If( j==Bj .or. j==Bjh )then
				Cycle
			End If
			X0 = Particle(j)%X
			Y0 = Particle(j)%Y
			R0 = Particle(j)%Radius
			Dence = Sqrt( (X0 - Particle(i)%X)**2 + (Y0 - Particle(i)%Y)**2 )
			Rd = R0 + Particle(i)%Radius
			If( Abs(Dence - Rd) > Rc .And. Dence < Rd ) then
				Qk = 5
				Exit
			End If
		End Do 
	End If
	Call Ding_Wei( Bj )

End Subroutine Check_Overlap


!!前进边更新
Subroutine Update_ForwardLine(i, Qk, Hbh, Qbh, Bj) 
    Implicit None
    Integer :: 			Qk, Hbh, Qbh, Bj, i, Tou, Wei, K, K1, Js, Bjq, Bjh


    Bj = ForwardLine%Bh 
	If( Bj/=-1 )then
		Bjq = ForwardLine%Prev%Bh
	End IF
	If( Bj/=-2 )then
		Bjh = ForwardLine%Next%Bh
	End IF
	Select Case( Qk )
	Case( 1 )
        Items%Item = Particles(Particle(i)%Xh, Particle(i)%X, Particle(i)%Y, Particle(i)%Radius)
        Call Ding_Wei(Bj)
        Call Line_Add(ForwardLine, Items)
	Case( 2 )
        Call Ding_Wei(Bj)
		ForwardLine => ForwardLine%Next
		Do 	
			If( .Not. Associated( ForwardLine%Next ) ) Exit
            Call Line_Dele( ForwardLine )    
			If( ForwardLine%Bh == Hbh ) Exit
        End Do     
        Call Ding_Wei( Bj )    
	Case( 3 )
        Call Ding_Wei( Qbh )    
		ForwardLine => ForwardLine%Next
		Do 
			If( .Not. Associated( ForwardLine%Next ) ) Exit
            Call Line_Dele( ForwardLine )    
            If( ForwardLine%Bh == Bjh ) Exit    
        End Do     
        Call Ding_Wei( Qbh )    
	Case( 4 )
		Call Ding_Wei( Qbh )
		ForwardLine => ForwardLine%Next
		Do 
			If( .Not. Associated( ForwardLine%Next ) ) Exit
			Call Line_Dele( ForwardLine )
			If( ForwardLine%Bh == Hbh ) Exit
		End Do 
		Call Ding_Wei( Qbh )
	End Select

End Subroutine Update_ForwardLine   


!前进边增加------调用前应 Allocate(Items), 在 Pos 位置后填加一条数据 
! 调用前后不改变原先指针所指的位置 
Subroutine Line_Add(Pos, Items)
    Implicit None
    Type(Datalink), Pointer :: Pos, Items
    Items%Next => Pos%Next
    Items%Prev => Pos
    If( Associated( Pos%Next ) )then
        Pos%Next%Prev => Items
    End If
    Pos%Next => Items
End Subroutine Line_Add


!前进边删除------将 Items 位置的数据删除 
! 调用后指针指向原 Items 的下一个数据
Subroutine Line_Dele(Items)
    Implicit None
    Type(Datalink), Pointer :: Items, Prev, Next
    Prev => Items%Prev
    Next => Items%Next
    Deallocate(Items)
    If( Associated(Prev) ) Prev%Next => Next
    If( Associated(Next) ) Next%Prev => Prev
    Allocate(Items)
    Items => Next
End Subroutine Line_Dele


Subroutine Ding_Wei(Dw)
    Implicit None
    Integer :: Dw
	ForwardLine => Head
    Do 
        If( ForwardLine%Bh == Dw ) Exit
		If( .Not. Associated(ForwardLine%Next) ) Exit
        ForwardLine => ForwardLine%Next
    End Do
End Subroutine Ding_Wei


Subroutine Get_Mesh()
	Implicit None
	Integer( kind=Myintk ) ::	 			i,j,N_ColumnOfMesh,N_RowOfMesh,M,N
	Real( kind=Myrealk ) :: 					MeshSize

	MeshSize= 2.0*Maxval( Radius )
	N_ColumnOfMesh = Int( Br/MeshSize )+2
	MeshSize = Br/Real(N_ColumnOfMesh-2)
	N_RowOfMesh = Int( Bt/MeshSize )+2
	i = Int( MeshSize**2/Pi/MinVal(Radius)**2 ) + 1
	Allocate( Mesh(0:N_RowOfMesh-1,0:N_ColumnOfMesh-1,i+1),SandInMesh(Number,2), Neighbor(Number,i) )
	Neighbor(:,:) = 0
	Mesh(:,:,:) = 0
	SandInMesh( :,: ) = 0
	Do i = 1, Number
		M = Floor( Particle(i)%Y / MeshSize ) + 1
		N = Floor( Particle(i)%X / MeshSize ) + 1
		SandInMesh( i,1 ) = M
		SandInMesh( i,2 ) = N
		j = Mesh(m,n,1)
		Mesh(m,n,j+2) = Particle(i)%xh
		Mesh(m,n,1) = Mesh(m,n,1) + 1
	End Do

End Subroutine Get_Mesh

Subroutine GetNeighbor()
	Implicit None
	Integer :: 					i, j, k, N, Col, Row, l, Njc
	Real( kind=Myrealk ) :: 	SmallValue, Lij, Cd
	Character :: 				Ch*40

	Open( 111,File='Neighbor' )
	SmallValue = 1.0E-1*Minval(Radius)
	Njc = 0
	Do i = 1, Number
		Do Col = SandInMesh(i,2)-1, SandInMesh(i,2)+1
			Do Row = SandInMesh(i,1)-1, SandInMesh(i,1)+1
				k = Mesh(Row,Col,1)
				If( k>0 )then
					Do l = 1, k
						j = Mesh(Row,Col,l+1)
						Lij = Sqrt( (Particle(i)%x-Particle(j)%x)**2 + (Particle(i)%y-Particle(j)%y)**2 )
						If( Abs(Lij-Particle(i)%Radius-Particle(j)%Radius)<SmallValue )then
							n = Neighbor(i,1)
							Neighbor(i,n+2) = j
							Neighbor(i,1) = Neighbor(i,1)+1
						End If
						If( i/=j .and. (Particle(i)%Radius+Particle(j)%Radius)-Lij>SmallValue )then
							Njc = Njc + 1
							Print*, Njc,i,j
						End If
					End DO 
				End If
			End Do 
		End Do 
		N = Neighbor(i,1)
		Write( Ch,* ) N+2
		Ch = '('//Trim(Adjustl(Ch))//'(i9))'
		Write( 111,ch ) i,N,Neighbor(i,2:N+1)
		Ch = ''
	End Do
	Close( 111 )

End Subroutine GetNeighbor

Subroutine FillingRate()
	Implicit None
	Integer ::					j 
	Real( kind=Myrealk ) :: 	V1, V2, Xmax=0.0, Ymax=0.0

    Do j = 1, Number
		If(Particle(j)%Y + Particle(j)%Radius <= Bt)then
            If( Particle(j)%X > Xmax )then
                Xmax = Particle(j)%X
            End If
            If( Particle(j)%Y > Ymax )then
                Ymax = Particle(j)%Y
            End If
            V1 = V1 + 3.14*Particle(j)%Radius * Particle(J)%Radius
		End If
    End Do
    V2 = Xmax*Ymax
	FillRate = V1/V2

End Subroutine FillingRate

Subroutine Sorting()
	Implicit None
	Real( Kind=Myrealk ) ::		Temp
	Integer ::					i, j

	Do i = 1, Size(Radius)
		Do j = i+1, Size(Radius)
			If( Radius(i)>Radius(j) )then
				Radius(i) = Temp
				Radius(i) = Radius(j)
				Radius(j) = Temp
			End If
		End Do 
	End Do 

End Subroutine Sorting

Subroutine Information()
	Implicit None

	Print*,'//----------------------------------------------------------------------//'
	Print*,'//Instruction of this Algorithm.'
	Print*,'//  1.This Algorithm is a unequal particle parking program.'
	Print*,'//  2.Main content please read the paper--Filling domains with disks'
	Print*,'//  3.The sorts and number of particles you want to use can be use'
	Print*,'//	  configure file or interactivly.'
	Print*,'//  4.Detail Information of configure file please read README file'
	Print*,'//----------------------------------------------------------------------//'

End Subroutine Information

End Module Fl
    
    
Program Open_Form
    Use Quan_Ju_Bian_Lian
    Use Fl
    Implicit None
    Integer :: 					Qk = 1, Hbh, Qbh, Bj, j, Fid, Err, Ntmp=1, iii
	Real(Kind = Myrealk) :: 	Start, Finish
	Logical ::					Lg
	Character :: 				Conf='N', C, WFormat*60

!//---------------------------------------------参数设置-----------------------------------------//
	Call Information()
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
		If( Err /= 0 .or. Accuracy>1.0 )then
			Print*,Accuracy
			Print*,'Input a real value'
			Read(*,*)
			Stop
		End If
		Read( 11,* , Iostat = Err ) Number
		If( Err /= 0 .or. Number<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'Number is:', Number
			Read(*,*)
			Stop
		End IF
		Allocate( Particle(Number) )
		Read( 11,*, Iostat = Err ) N_Sortofballs
		If( Err /= 0 .or. N_Sortofballs<0 )then
			Print*,'Configure file maybe wrong'
			Print*,'N_Sortofballs is:', N_Sortofballs
			Read(*,*)
			Stop
		End IF
		Read( 11,* ) Bt
		Read( 11,* ) Bt
		Read( 11,* ) Bl
		Read( 11,* ) Br
		Allocate( Radius(N_SortofBalls),Ratio(N_SortofBalls) )
		Radius(:) = 0.0
		Ratio(:) = 0.0
		Read( 11,* ) Radius(:)
		Read( 11,* ) Ratio(:)
		Rinf = 1.0E5
		Accuracy = Accuracy*Radius(1)
		Call Sorting()
		Do i = 1, Number
			Read( 11,*,Iostat=Err ) Particle(i)%Xh,Particle(i)%Radius
			If( Err/=0 )then
				Print*,'Configure file maybe wrong'
				Print*,'Number of particle in the file is not enough'
				Read(*,*)
				Stop
			End If
		End Do 
		Close( 11 )
	Else
		Call CPU_TIME(Start)
		Print*
		Print*,'Start.....'
		Print*
		Print*,'Input the Newton ineration Accuracy '
		Do 
			Read(*,*, Iostat = Err) Accuracy
			If( Err /= 0 .or. Accuracy>1.0 )then
				Print*,'Input a real value'
				cycle
			Else
				Exit
			End If
		End Do 
		Print*,'Input the total number of sandbed'
		Do 
			Read(*,*, Iostat = Err) Number
			If( Err /= 0 )then
				Print*,'Input must be integer'
				cycle
			Else
				Exit
			End If
		End Do 
		Allocate( Particle(Number) )
		print*,'Input the sorts of sand: '
		Do 
			Read(*,*, Iostat = Err) N_Sortofballs
			If( Err /= 0 )then
				Print*,'Input must be Integer'
				cycle
			Else
				Exit
			End IF
		End Do 
		Allocate( Radius(N_SortofBalls),Ratio(N_SortofBalls) )
		Radius(:) = 0.0
		Ratio(:) = 0.0
		Print*,'Input the radius of sands(Unit/M)'
		Do i = 1, N_SortofBalls
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
		!	Ratio = (/0.3438,0.5337,0.1138,0.0087/)
		!	Radius = (/0.03E-3,0.05E-3,0.1E-3,0.15E-3/)
		Bb = 0.0
		Bl = 0.0
		Rinf = 1.0E5
		Print*,'Input the edge of sandbed--Br,Bt'
		Read(*,*) Br
		Read(*,*) Bt
		!Call Get_RatioLog(N_SortofBalls)
		Call Get_Ratio(N_SortofBalls)
		Call Random_Radius()
		Call Sorting()
	End If !//---------------------------------------------开始计算-----------------------------------------//
	Print*,Radius
	Read(*,*)
	Print*
	Print*,'Start to filling domains'
    Call CornerParticle
	Print*
	Print*,'Corner particle ok'
    Call Initialize_ForwardLine
	Print*
	Print*,'Initialization of forwardline ok'
    i = 1
    Do 
        i = i + 1
        If( i > Number ) Exit
        Call Particle_Get(i)
        Call Check_Overlap(i, Qk, Hbh, Qbh, Bj)
		If( Qk==5 )then
			i = i - 1
			ForwardLine => ForwardLine%Next
			Cycle
		End IF
        Call Update_ForwardLine(i, Qk, Hbh, Qbh, Bj)
        If( Qk /= 1 )then
!			Do iii = 1, Size(Radius)
!				If( Particle(i)%Radius>Radius(iii) )then
!					Particle(i)%Radius = Radius(iii)
!				End If
!			End Do 
			i = i - 1
			Cycle
		Else 
			If(ForwardLine%Next%Next%Bh == -3)then
				ForwardLine => ForwardLine%Next	
			ElseIf(ForwardLine%Next%Next%Bh == -2)then
				ForwardLine => Head
			Else
				ForwardLine => ForwardLine%Next%Next	
			End If
        End If
		If( Particle(i)%Y<Bt )then
			Ntmp = Ntmp + 1
		End If
	End Do 
	Print*
	Print*,'Filling ok'
	Open( 11, File='Configure' )
	Open( 12, File='Particle' )
	Write( 11,"(G10.3,2x,A30)" ) Accuracy, 'Newton ineration Accuracy'
	Write( 11,"(i10,2x,A30)" ) Ntmp, 'Particle Number'
	Write( 11,"(i10,2x,A30)" ) N_Sortofballs, 'Particle sorts'
	Write( 11,"(G10.3,2x,A30)" ) Bb,'Bottom edge'
	Write( 11,'(G10.3,2x,A30)' ) Bt,'Top edge'
	Write( 11,'(G10.3,2x,A30)' ) Bl,'Left edge'
	Write( 11,'(G10.3,2x,A30)' ) Br,'Right edge'
	Write( WFormat,* ) N_Sortofballs
	WFormat = '('//Trim(Adjustl(WFormat))//'(G10.3),2x,A30)'
	Write( 11,WFormat ) Radius,'Radius'
	Write( 11,WFormat ) Ratio, 'Ratio'
    Do j = 1, Number
		If( Particle(j)%Y<Bt )then
			Write(11, "(I6, 2(f13.6))") Particle(j)%Xh, Particle(j)%Radius
			Write(12, "(I6, 3(f13.6))") Particle(j)%Xh, Particle(j)%X, Particle(j)%Y, Particle(j)%Radius
		End If
    End Do
	Close(11)
	Close(12)
	Call Get_Mesh()
	Call GetNeighbor()
	Call FillingRate()
	Write(*,"(A30,i13)") 'Number in the area is:', Ntmp
	Write(*,"(A30,F13.6)") 'Filling rate is:', FillRate
	Call CPU_TIME(Finish)
	Write(*,"(A30,F13.6)") 'Runtime is(second):', Finish-Start
    !Call system('"/usr/local/bin/matlab/bin/matlab" -r gif')
	Print*
	Print*,'Finish'
	Read(*, *)
End Program Open_Form    
