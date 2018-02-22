!******************************************************************************************************************************
!                                   程序说明                           
!   1. 该程序比较了基于两颗粒正碰的三种碰撞模型。分别为：弹簧模型 赫兹模型 弹塑性模型。
!   2. 模型介绍：
!           颗粒大小为 CM 级别。
!           碰撞速度为 m/s 
!           下球固定，上球下落与之碰撞。
!           上球为球1，下球为球2
!   3.注意： 模型单位用 m。
!           
!******************************************************************************************************************************
    
    
Module Mla
    Implicit None
    Integer, Parameter ::  MyRealk = Selected_Real_Kind(20,10)
    Integer :: Mx
    Real(kind=Myrealk), Parameter :: G=9.8, Den=2200, Ym=0.5E9, P=0.4, Pai=3.1415926, K_Spring=0.05E8
    Real(Kind=Myrealk) :: Deltamax=0.0, R_Equ=0.0, V1=0.0, V2=0.0, K12=0.0, V_Impact=0.0, R1, R2, R_Impact, &
                            Epsilon, Deltt, Deltaq, Deltamax_Zhou, Deltaxing=0.0, Delta_P=0.0
Contains
    !******************************************以下代码计算两球碰撞的最大压缩量**********************************************
Subroutine Ya_Suo_Lian
    Implicit None    
    Real(Kind=Myrealk) :: X1, X2, Gama=0.0, Gama1=0.0, Fn1=0.0, Fn2=0.0, Lamda, Delta, M1, M2, M, Deltg, Deltaa, Y1, Y2    
    Real(Kind=Myrealk) :: Const=0.0, Bn=0.0, X=0.0, Sigma_C, Sigma_S, Sigma_Y, Ra, A_Y, A_P, V_Y, A_P1, A_P2
    Real(kind = Myrealk), Parameter :: A = 0.8**2.5, B = 1.0 - 0.8**2.5    
    Logical :: Flag=.True. , Flag1 = .False.   
    Integer :: i, jj=0              
   
    i = 0     
    R2 = 10.0E-3
    R_Equ = R1*R2/(R1+R2)                           !等效粒径    
    M1 = Den*4.0*Pai*R1**3/3.0                      !上球的质量    
    M2 = Den*4.0*Pai*R2**3/3.0                      !下球的质量    
    M =  M1*M2 / (M1 + M2)    
    x1 = 0.0                                        !x1， y1， x2， y2 为两球的位置。    
    x2 = 0.0    
    Y1 = R1 + 2*R2    
    Y2 = R2                                         !两球的碰撞速度    
    K12 = 4.0*Ym*R_Equ**0.5/3.0                     !两球碰撞的弹性模量    
    Flag = .True.    
    jj = 0
    Gama = 0.0
    Gama1 = 0.0
    V1 = V_Impact
    Sigma_C = 0.0
    !Deltamax_Zhou = (1.25*M*V_Impact**2/K12)**0.4
    Deltamax_Zhou = ((R1+R2)/R1)*(1.25*M1*V_Impact**2/K12)**0.4
    
    Delta_P = 0.0
    
    Do     
        i = i + 1
        Gama = (Y1-Y2) - (R1 + R2)    
        Delta = Abs(Gama)    
        If( Gama<=0.0 )then    
            Deltg= Abs(Gama)-Gama1    
            Deltaa = Deltg/Deltt    
            If( Deltaa<0.0 .And. Flag==.True. )then    
                Deltamax = Abs(Gama1)    
                jj = jj + 1    
                Flag = .False.    
            End If    
            Select Case( Mx )    
            Case( 1 )    
	            Fn1 = k12*Abs(Gama)**1.5	
            Case( 2 )    
                Const = (log(Epsilon))**2                
                Bn = Sqrt((4.0*M1*k_Spring*Const) / (Pai*Pai + Const))    
                Fn1 = k_Spring*Abs(Gama) + Bn*Deltaa    
            Case( 3 )
                x = 1.0 - 0.8*Epsilon*Epsilon - (0.2*Epsilon**8*A)/(1.0 - (1.0 - A)*Epsilon) + Epsilon**16*(1.0-Epsilon)
                DeltaXing = Deltamax_Zhou * (1.0-X) 
     
                If( Deltaa >= 0.0 )then    
                    Fn1 = K12*Delta**1.5    
                Elseif( Deltaa < 0.0 )then    
                    Lamda = K12*Deltamax_Zhou**0.5*(1.25/Epsilon/Epsilon*(1.0 + Epsilon**6*0.25*A/(1.0 - B*Epsilon) - 1.25*Epsilon**14*(1.0-Epsilon))**(-1) &    
                            - 0.8**(0.5)*(1.0 + (Epsilon**6*0.25*A)/(1.0 - B*Epsilon) - 1.25*Epsilon**14*(1.0 - Epsilon) )**(0.5))    
                    If( Delta < DeltaXing )then    
                        Fn1 = 0.0  
                        Flag1 = .True.
                    Elseif( Delta >= DeltaXing )then    
                        Fn1 = Lamda*(Delta - DeltaXing) + Epsilon**3*K12*(Delta - DeltaXing)**1.5    
                    End If    
                End if 
            Case( 4 )
                Sigma_Y = 23.0E6
                V_Y = (0.5*Pai/Ym)**2*(8.0*Pai*R_Equ**3/15.0/M)**0.5*Sigma_Y**2.5
                A_Y = (15.0*R_Equ**2*M*V_Y**2/(16.0*Ym))**0.2                
                If( Deltaa >= 0.0 )then
                    If(Sigma_C <= Sigma_Y)then
                        Fn1 = 4.0*Ym*Sqrt(R_Equ)*Delta**1.5/3.0
                        A_P1 = Gama1
                    Else
                        Fn1 = -(Sigma_Y*Pai)**3*R_Equ**2/(12.0*Ym**2) + Pai*R_Equ*Sigma_Y*Delta
                        A_P2 = Gama1
                    End If
                    If(Delta /= 0.0)then
                        Ra = Sqrt(Delta*R_Equ)
                        Sigma_C = 3.0*Fn1/(2.0*Pai*Ra**2)
                    Else
                        Ra = 0.0
                        Sigma_C = 0.0
                    End If
                Else
                    !A_P = Sqrt(Delta**2 - A_Y**2)
                    A_P = A_P2 - A_P1
                    Delta_P = A_P**2/R_Equ
                    If(Delta > Delta_P)then
                        Fn1 = 4.0*Ym*Sqrt(R_Equ)*(Delta-Delta_P)**1.5/3.0
                    Else
                        Fn1 = 0.0
                        Flag1 = .True.
                    End If
                End IF
                DeltaXing = Delta_P 
            End Select     
        Else    
            Flag = .True.    
        End If
        V1 = V1 + (Fn1/M1-G)*Deltt    
        Y1 = Y1 + V1*Deltt + 0.5*(Fn1/M1-G)*Deltt**2
        
        !If( Mx == 4 ) Write(25, "(i10, 2x, E20.11E4)") i, Y1
             
        Fn1 = 0.0  
        Gama1 = Abs(Gama)    
        If ( (Gama>0.0 .And. V1>0.0) .or. (Flag1==.True. .and. Flag==.False.) )then
            Exit 
        End If
    End Do 
     
    
    
End Subroutine Ya_Suo_Lian
    
!//函数积分求解
Subroutine Integral(Inte)
	Implicit None
	Real(Kind = Myrealk) :: X0, X1, X, Bc, Inte, Pai1
	Real(kind = Myrealk), Parameter :: Zeta0 = 3.2E-8, Eps0 = 8.85E-12 ,Kb = 1.38E-23, T = 273.16 + 22.0, Q = 1.6E-19
	Inte = 0.0
	X0 = 1.0
	X1 = 1.0E5
	Bc = 1.0E-2
	X = X0
	Do 
		X = X + Bc
		Pai1 = Exp((-Q*Q*X)/(4*Pai*Zeta0*Eps0*Kb*T))/X**2
		If(X >= X1) Exit
		Inte = Inte + Pai1*Bc
	End Do
End Subroutine Integral



!//碰撞电的计算
Subroutine Electric_Calculation(Pai_D)
	Implicit None
	Real(Kind = Myrealk), Parameter ::  N = 1.0E18, Q = 1.6E-19, Kb = 1.38E-23, Den = 2.2E3,&
										Zeta0 = 3.2E-8, Eps0 = 8.85E-12 ,Est = 1013.25, Tst = 373.15,&
										Mol = 6.02E23, H = 6.63E-34, Alpha = 0.5, Rho=0.006
	Real(Kind = Myrealk) :: Pd1, Pd2, A1, A2, Lamda3, Ps1, Ps, P, T, Xi0, Rh, X, M, Pd11, Pd12, Pd13, Pai_D, Delta_P	
    
	T = 273.16 + 22.0
	Xi0 = -30.0*1000
	Rh = 20
	M = 1.67E-27
	Lamda3 =(2.0*Pai*M*Kb*T/H**2)**(-1.5)
	Ps1 = -7.90298*(Tst/T-1.0) + 5.02808*Log10(Tst/T) - 1.3816*1.0E-7*(10**(11.344*(1.0-T/Tst))-1.0) + &
			8.1328*1.0E-3*(10**(-3.49149*(Tst/T-1.0))-1.0) + Log10(Est)
	Ps = 100*10.0**(Ps1)
	P = Rh*Ps/100.0
	Pd11 = (Kb*T/Lamda3)
	Pd12 = (Pai_D/P)
	Pd13 = Exp(Xi0/(Kb*T*Mol))
	Pd1 = 1.0/(Pd11*Pd12*Pd13 + 1)
	Pd2 = Pd1
    Delta_P = 0.0
    
    If(Mx == 3 .or. Mx == 4)then
        Delta_P = DeltaXing
        If(Deltamax > DeltaXing)then
            Delta_P = DeltaXing
            A1 = 2*Pai*R1**2*(Sqrt(1.0 - R_Equ*Delta_P/R1**2) - Sqrt(1.0 - R_Equ*Deltamax/R1**2))
            A2 = 2*Pai*R2**2*(Sqrt(1.0 - R_Equ*Delta_P/R2**2) - Sqrt(1.0 - R_Equ*Deltamax/R2**2))
        Else
            A1 = 2*Pai*R1**2*(1.0-Sqrt(1.0 - R_Equ*Deltamax/R1**2))
            A2 = 2*Pai*R2**2*(1.0-Sqrt(1.0 - R_Equ*Deltamax/R2**2))
        End If
    !Elseif(Mx == 4)then
    !    Delta_P = DeltaXing
    !    A1 = 2*Pai*R1**2*(Sqrt(1.0 - R_Equ*Delta_P/R1**2) - Sqrt(1.0 - R_Equ*Deltamax/R1**2))
    !    A2 = 2*Pai*R2**2*(Sqrt(1.0 - R_Equ*Delta_P/R2**2) - Sqrt(1.0 - R_Equ*Deltamax/R2**2))
    Else
	    A1 = 2*Pai*R1**2*(1.0-Sqrt(1.0 - R_Equ*Deltamax/R1**2))
	    A2 = 2*Pai*R2**2*(1.0-Sqrt(1.0 - R_Equ*Deltamax/R2**2))
    End If
    
	Deltaq = Alpha*Rho*Pd1*(1-Pd1)*(A2-A1)
    
    Select Case( Mx )
    Case( 1 )
        Write(21,"(F13.6,5x,E13.6)")  R1/R2, -Deltaq    !R1/R2, -Deltaq
    Case( 2 )
        Write(22,"(F13.6,5x,E13.6)")  R1/R2, -Deltaq    !R1/R2, -Deltaq
    Case( 3 )
        Write(23,"(F13.6,5x,E13.6)")  R1/R2, -Deltaq    !R1/R2, -Deltaq
    Case( 4 )
        Write(24,"(F13.6,5x,E13.6)")  R1/R2, -Deltaq    !R1/R2, -Deltaq
    End Select

    
End Subroutine Electric_Calculation
    
End Module Mla
    
    
Program Main
    Use Mla
    Implicit none	
	Integer :: i
	Real(kind = MyRealk) :: Start, Finish, Pai_D
	 
	Print*
	Print*,'**********此程序计算的是非等粒径沙床碰撞过程**********'
	Print*
	Open(21, File = '10-1.txt', Access = 'Append')
    Open(22, File = '10-2.txt', Access = 'Append')
    Open(23, File = '10-3.txt', Access = 'Append')
    Open(24, File = '10-4.txt', Access = 'Append')
    !Open(25, File = 'Location.txt', Access = 'Append')
    
    Call Integral(Pai_D)
    
    Mx = 0
    Do Mx = 1, 4 
        R_Impact = 0.0E-3
        Epsilon = 0.4
        R1 = R_Impact
        V_Impact = -3.96
        Print*, Mx
        Print* 
         
        DO while(R1 <= 9.9E-3)
            R1 = R1 + 0.1E-3
            V_Impact = -3.96
            !V_Impact = V_Impact - 0.1
            Deltt = 1.0E-7 
            Call Ya_Suo_Lian
            Call Electric_Calculation(Pai_D)
        End do 

        Print*, Deltamax
        Deltamax = 0.0
    End Do 

    Close(21)
	Close(22)
    Close(23)
    Close(24)
    

End program Main