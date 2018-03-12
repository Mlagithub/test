Module Wind
	Use GlobalVariable
	Implicit None

    Real( kind=MyRk ), Parameter :: 		Rho_Air=1.29 !// Density of air  Ga=1.5E-5, Kaman=0.4, Z0=9.3E-4
    Real( kind=MyRk ), Parameter :: 		Ga=1.5E-5 !// Air viscous coefficient
    Real( kind=MyRk ), Parameter :: 		Kaman=0.4 !// Kaman constant
    Real( kind=MyRk ), Parameter :: 		Z0=9.3E-4 !// 下垫面粗糙度
	Real( kind=MyRk ), Allocatable :: 		TauA(:) !// Air shear stress
	Real( kind=MyRk ), Allocatable :: 		TauP(:) !// Sand shear stress
	Real( kind=MyRk ), Allocatable :: 		Zwind(:) !// Wind velocity at different height
	Real( kind=MyRk ), Allocatable :: 		U_New(:) !// Frictional wind speed at different height 
	Real( kind=MyRk ), Allocatable :: 		Qm(:) !// Q(electricity)/M(mass)
	Real( kind=MyRk ), Allocatable :: 		QmSum(:) !// Q(electricity)/M(mass)
	Real( kind=MyRk ) ::					H_Qm !// Frictional wind speed
	Real( kind=MyRk ) ::					CellOfQm !// Step length of calculaion Q/M
	Real( kind=MyRk ) ::					Ustar, Ustaro !// Frictional wind speed
	Real( kind=MyRk ) ::					UstarT !// Frictional wind speed
	Real( kind=MyRk ) ::					Tau !// Total shear stress 
	Real( kind=MyRk ) ::					TauT !// Critical shear stress
	Real( kind=MyRk ) ::					Tau_Z0
	Real( kind=MyRk ) ::					Delta_Z=0.01
	Real( kind=MyRk ) ::					Delta_Z_Acc=0.001
	Real( kind=MyRk ) ::					Height_Z_Acc=0.02
	Integer( kind=MyIk ) ::					N_Z
	Integer( kind=MyIk ) ::					N_Z_Acc=20


Contains


Subroutine ModifyWind()
    Implicit None
    Real( Kind=MyRk ) ::                    Z_Surface_Temp(5000), Delta_Z_Surface, Z_Surface_Top, Z_Surface_Down, TauA_Surface, H, Pp
    Real( Kind=MyRk ), Allocatable ::       Z_Surface(:), U_Surface(:), Dudz_Surface(:), TauP_Surface(:), Dudz(:)
    Integer( Kind=MyIk ) ::					N, i, j, k, Remainder, K1, K2, j1, Ws

    !====================================计算空气剪切力===============================
	ws = N_Z+N_Z_Acc
	TauP(:) = TauP(:)/Edgeright
	Call Interpolation()
	Allocate( Dudz(Ws) )
	TauA(:) = Tau
    Dudz(:) = 0.0
    !====================================计算风速廓线的一阶导数 Dudz ==================
	Forall( i=1:N_Z+N_Z_Acc )
		TauA(i) = Max( 1.0E-10,Tau-TauP(i) )
        Dudz(i) = ( TauA(i)/(Rho_Air*Kaman**2*Zwind(i)**2) )**0.5
	End Forall
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
	Do i = 1, 28!ws
		Ustar = Ustar + U_New(i)*Kaman/log(Zwind(i)/Z0)
	End Do 
	Ustar = Ustar/28.0 !Real(ws)
	If( Ustar>Ustaro ) Ustar = 0.95*Ustaro
	Deallocate( Dudz,Z_Surface,U_Surface,Dudz_Surface,TauP_Surface )
	TauP(:) = 0.0
    
End Subroutine ModifyWind

Subroutine Get_TauP( i )
	Implicit None
	Integer( Kind=MyIk ) ::  	i, Zi

	If( Pt(i)%Y-H_SandBed>Delta_Z_Acc*N_Z_Acc )then
		Zi = N_Z_Acc + Int( (Pt(i)%Y-H_SandBed-Delta_Z_Acc*N_Z_Acc)/Delta_Z ) + 1
	Else
		Zi = Int( (Pt(i)%Y-H_SandBed)/Delta_Z_Acc ) + 1
	End If
	TauP(zi) = TauP(zi) - Pt(i)%vy/Abs(Pt(i)%Vy)*Pt(i)%Mass*Pt(i)%Vx
	

End Subroutine Get_TauP

Subroutine Interpolation()
	Implicit None
	Integer( Kind=MyIk ) :: 		i, j, k, M
	
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



Subroutine Initialize_Wind()
	Implicit None
	Integer( kind=MyIk ) :: i

	N_Z = Int( (Height-N_Z_Acc*Delta_Z_Acc-H_SandBed)/Delta_Z )
	Allocate( Zwind(N_Z+N_Z_Acc),TauA(N_Z_Acc+N_Z),TauP(N_Z+N_Z_Acc),U_New(N_Z+N_Z_Acc) )
	TauA(:) = 0.0
	TauP(:) = 0.0
	U_New(:) = 0.0
	Do i = 1, N_Z_Acc
		Zwind(i) = Real( Delta_Z_Acc*i )
	End Do 
	Do i = 1+N_Z_Acc, N_Z+N_Z_Acc
		Zwind(i) = Real( Delta_Z*Float(i-N_Z_Acc) ) + Zwind(N_Z_Acc)
	End Do 

End Subroutine Initialize_Wind


!// Drag Force 
Real Function DragX( i )
	Implicit None
	Integer, Intent(In) ::	i
	Real( kind=MyRk ) ::	Vwind, Cd, Re, V_Relat

	V_Relat = ( (Pt_old(i)%Vx-Vwind)**2 + (Pt_Old(i)%Vy)**2 )**0.5
	Re = V_Relat*Pt(i)%Radius*2.0/Ga
	Cd = 24.0/Re + 6.0/(1.0+(Re)**0.5) + 0.4
	Vwind = Ustar*Log((Pt_Old(i)%Y-H_SandBed)/Z0)/Kaman
	DragX = -(Pi/8.0)*(4.0*Pt(i)%Radius**2)*Rho_Air*Cd*(Pt_Old(i)%Vx-Vwind)*V_Relat

End Function DragX
	

Real Function DragZ( i )
	Implicit None
	Integer, Intent(In) ::	i
	Real( kind=MyRk ) ::	Vwind, Cd, Re, V_Relat

	V_Relat = ( (Pt_old(i)%Vx-Vwind)**2 + (Pt_Old(i)%Vy)**2 )**0.5
	Re = V_Relat*Pt(i)%Radius*2.0/Ga
	Cd = 24.0/Re + 6.0/(1.0+(Re)**0.5) + 0.4
	Vwind = Ustar*Log((Pt_Old(i)%Y-H_SandBed)/Z0)/Kaman
	DragZ = -(Pi/8.0)*(4.0*Pt(i)%Radius**2)*Rho_Air*Cd*Pt_Old(i)%Vy*V_Relat

End Function DragZ



End Module Wind
