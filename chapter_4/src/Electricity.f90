Module Electricity
	Use GlobalVariable
	Implicit None
	Real( kind=MyRk ), Allocatable ::                       ComMax(:,:)
	Real( Kind=MyRk ) ::					Temperature !// Global temperature
	Real( Kind=MyRk ) ::					Rh !// Relativy humidity
	Real( Kind=MyRk ) ::					Xfn !// Adsorption energy
	Real( Kind=MyRk ) ::					Pdd



Contains


!//Solve function integral
Subroutine Integral( T,Inte )
	Implicit None
	Real( kind=8 ), Parameter ::		Zeta0=3.2E-8, Eps0=8.85E-12 ,Kb=1.38E-23, Q=1.6E-19
	Real( kind=MyRk ), Intent(In) ::  	T
	Real( kind=8 ), Intent(out) ::		Inte
	Real( kind=8 ) ::					X0, X1, X, Bc, Pai1

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



!// Calculation of electron transfer probability
Real Function Pd( Ttemp,Rh,Xfn )
	Real( kind=8 ), Parameter ::		Zeta0=3.2E-8, Eps0=8.85E-12 ,Est=1013.25
	Real( kind=8 ), Parameter ::		Mol=6.02E23, H=6.63E-34, Tst=373.15, Kb=1.38E-23
	Real( kind=8 ) ::					Ps1, Ps, P, Xi0, X, M, Pd11, Pd12, Pd13, Lamda3
	Real( kind=8 ) ::					Pai
	Real( kind=MyRk ) ::				T !// Temperature
	Real( kind=MyRk ), Intent(In) ::	Ttemp !// Temperature
	Real( kind=MyRk ), Intent(In) ::  	Rh !// Relativy humidity
	Real( kind=MyRk ), Intent(In) ::  	Xfn !// Adsorption energy
   
	T = 273.16 + Ttemp
	Call Integral( T=T,Inte=Pai )
	Pai = Real( Pai,MyRk )
	Xi0 = Xfn*1000
	M = 1.67E-27
	Lamda3 = (2.0*Pi*M*Kb*T/H**2)**(-1.5)
	Ps1 = -7.90298*(Tst/T-1.0) + 5.02808*Log10(Tst/T) - 1.3816*1.0E-7*(10**(11.344*(1.0-T/Tst))-1.0) + &
			8.1328*1.0E-3*(10**(-3.49149*(Tst/T-1.0))-1.0) + Log10(Est)
	Ps = 100*10.0**(Ps1)
	P = Rh*Ps/100.0
	Pd11 = (Kb*T/Lamda3)
	Pd12 = (Pai/P)
	Pd13 = Exp(Xi0/(Kb*T*Mol))
	Pd = 1.0/(Pd11*Pd12*Pd13 + 1)

End Function Pd


!// Calculation of collision charge of two particles
Real Function DeltaQ( i,j )
	Implicit None
	Real( Kind=MyRk ), Parameter ::  			N=1.0E18, Q=1.6E-19
	Real( Kind=MyRk ) ::						A1, A2, Ai, Aj, Requiv
	Integer( kind=MyIk ), Intent(In) :: 		i, j


	Requiv = Pt(i)%Radius * Pt(j)%Radius / ( Pt(i)%Radius + Pt(j)%Radius )

	If( 2.0*Requiv*ComMax(i,j)/Pt(i)%Radius**2>1.0 )then
		A1 = 2*Pi*Pt(i)%Radius**2
		Print*,'E maybe too big'
	Else
		A1 = 2*Pi*Pt(i)%Radius**2*(1.0-Sqrt(1.0 - 2.0*Requiv*ComMax(i,j)/Pt(i)%Radius**2))
	End If
	If( 2.0*Requiv*ComMax(i,j)/Pt(j)%Radius**2>1.0 )then
		A2 = 2*Pi*Pt(j)%Radius**2
		Print*,'E maybe too big'
	Else
		A2 = 2*Pi*Pt(j)%Radius**2*(1.0-Sqrt(1.0 - 2.0*Requiv*ComMax(i,j)/Pt(j)%Radius**2))
	End If

	Ai = 4.0*Pi*Pt(i)%Radius**2
	Aj = 4.0*Pi*Pt(j)%Radius**2
	DeltaQ = 0.002*( Pt(j)%Alpha*A2*Pt(j)%Pd*(1.0-Pt(i)%Pd) - Pt(i)%Alpha*A1*Pt(i)%Pd*(1.0-Pt(j)%Pd) )

	Pt(i)%Co = Pt(i)%Co + 1
	Pt(j)%Co = Pt(j)%Co + 1

	Pt(i)%Alpha = 0.132*0.8**Pt(i)%Co
	Pt(j)%Alpha = 0.132*0.8**Pt(j)%Co

	Pt(i)%Pd = Pt(i)%Pd + Pt(i)%Alpha*&
		( A2*Pt(j)%Pd*(1.0-Pt(i)%Pd) - A1*Pt(i)%Pd*(1.0-Pt(j)%Pd) ) / Ai
	Pt(j)%Pd = Pt(j)%Pd - Pt(i)%Alpha*&
		( A2*Pt(j)%Pd*(1.0-Pt(i)%Pd) - A1*Pt(i)%Pd*(1.0-Pt(j)%Pd) ) / Aj

End Function DeltaQ




End Module Electricity
