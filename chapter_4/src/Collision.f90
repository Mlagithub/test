Module Collision
	Use GlobalVariable
	Use MeshGrid
	Use Wind
	Use Electricity
	Implicit None

Contains



!// Two-dimensional normal contact spring force
!// F(orce)n(ormal)S(pring)2(two)d(imension)_H(ert)z
Real Function FnS2d_Hz( i,j,DeltaN,E1,E2,Miv1,Miv2  )
	Implicit None
	Integer(kind=MyIk) ::		i !// Particle 1
	Integer(kind=MyIk) ::		j !// Particle 2
	Real( kind=MyRk ) ::		R1 !// Radius of particle 1
	Real( kind=MyRk ) ::		R2 !// Radius of particle 2
	Real( kind=MyRk ) ::		E1 !// Young's modulus of particle 1
	Real( kind=MyRk ) ::		E2 !// Young's modulus of particle 2
	Real( kind=MyRk ) ::		Miv1 !// Poisson ratio of particle 1
	Real( kind=MyRk ) ::		Miv2 !// Poisson ratio of particle 2
	Real( kind=MyRk ) ::		DeltaN !// Normal compression of two particles
	Real( kind=MyRk ) ::		Requiv !// Equivalent radius
	Real( kind=MyRk ) ::		Eequiv !// Equivalent Young's modulus
	Real( kind=MyRk ) ::		K !// Hertz contact constant
	

	R1 = Pt(i)%Radius
	R2 = Pt(j)%Radius
	Requiv = R1*R2/(R1+R2)
	Eequiv = ( (1.0-Miv1*Miv1)/E1 + (1.0-Miv2*Miv2)/E2 )**(-1)
	K = 4.0*Eequiv*Requiv**0.5/3.0
	FnS2d_Hz = K * DeltaN**1.5

End Function FnS2d_Hz



!// Two-dimensional normal contact spring force
!// F(orce)n(ormal)S(pring)2(two)d(imension)_Sp(ring)
Real Function FnS2d_Sp( K,DeltaN )
	Implicit None
	Real( kind=MyRk ) :: DeltaN !// Normal compression of two particles
	Real( kind=MyRk ) :: K !// Spring stiffness

	FnS2d_Sp = K * DeltaN

End Function FnS2d_Sp



!// Two-dimensional normal contact damping force
!// F(orce)n(ormal)D(amping)2(two)d(imension)_M(in)d(lin)
Real Function FnD2d_Md( i,j,E1,E2,Miv1,Miv2,DeltaN,VnRel )
	Implicit None
	!// Damping constant
	Real( kind=MyRk ),Parameter ::		Beta=0.3033156 !Beta = log(2.71829) / sqrt( (log(2.71828))**2+Pi**2 )
	Integer(kind=MyIk) ::		i !// Particle 1
	Integer(kind=MyIk) ::		j !// Particle 2
	Real( kind=MyRk ) ::		VnRel !// Relative velocity
	Real( kind=MyRk ) ::		S !// Tangential stiffness
	Real( kind=MyRk ) ::		M !// Equivalent mass
	Real( kind=MyRk ) ::		R1 !// Radius of particle 1
	Real( kind=MyRk ) ::		R2 !// Radius of particle 2
	Real( kind=MyRk ) ::		M1 !// Mass of particle 1
	Real( kind=MyRk ) ::		M2 !// Mass of particle 2
	Real( kind=MyRk ) ::		E1 !// Young's modulus of particle 1
	Real( kind=MyRk ) ::		E2 !// Young's modulus of particle 2
	Real( kind=MyRk ) ::		Miv1 !// Poisson ratio of particle 1
	Real( kind=MyRk ) ::		Miv2 !// Poisson ratio of particle 2
	Real( kind=MyRk ) ::		Requiv !// Equivalent shear modulus
	Real( kind=MyRk ) ::		Eequiv !// Equivalent Young's modulus
	Real( kind=MyRk ) ::		DeltaN !// Normal compression of two particles

	R1 = Pt(i)%Radius
	R2 = Pt(j)%Radius
	M1 = Pt(i)%Mass
	M2 = Pt(j)%Mass
	Requiv = R1*R2/(R1+R2)
	M = M1*M2/(M1+M2)
	Eequiv = ( (1.0-Miv*Miv)/E1 + (1.0-Miv*Miv)/E2 )**(-1)
	S = 2.0 * Eequiv * sqrt( Requiv*DeltaN )
	FnD2d_Md = -2.0 * sqrt(5.0/6.0) * Beta * sqrt(S*M)*Abs(VnRel)

End Function FnD2d_Md



!// Two-dimensional tangential contact spring force
!// F(orce)t(angential)S(pring)2(two)d(imension)_M(in)d(lin)
Real Function FtS2d_Md( i,j,E1,E2,Miv1,Miv2,DeltaN,DeltaT )
	Implicit None
	Integer(kind=MyIk) ::		i !// Particle 1
	Integer(kind=MyIk) ::		j !// Particle 2
	Real( kind=MyRk ) ::		S !// Tangential stiffness
	Real( kind=MyRk ) ::		R1 !// Radius of particle 1
	Real( kind=MyRk ) ::		R2 !// Radius of particle 2
	Real( kind=MyRk ) ::		G1 !// Shear modulus of particle 1
	Real( kind=MyRk ) ::		G2 !// Shear modulus of particle 2
	Real( kind=MyRk ) ::		E1 !// Young modulus of particle 1
	Real( kind=MyRk ) ::		E2 !// Young modulus of particle 2
	Real( kind=MyRk ) ::		Miv1 !// Poisson ratio of particle 1
	Real( kind=MyRk ) ::		Miv2 !// Poisson ratio of particle 2
	Real( kind=MyRk ) ::		Gequiv !// Equivalent shear modulus
	Real( kind=MyRk ) ::		Requiv !// Equivalent Radius
	Real( kind=MyRk ) ::		DeltaN !// Normal compression of two particles
	Real( kind=MyRk ) ::		DeltaT !// Tangential compression of two particles

	R1 = Pt(i)%Radius
	R2 = Pt(j)%Radius
	Requiv = R1*R2/(R1+R2)
	G1 = E1/2.0/(1.0+Miv1)
	G2 = E2/2.0/(1.0+Miv2)
	Gequiv = (2.0-Miv1)/G1 + (2.0-Miv2)/G2 
	S = 8.0 * Gequiv * sqrt( Requiv*DeltaN )
	FtS2d_Md = S * DeltaT


End Function FtS2d_Md



!// Two-dimensional tangential contact damping force
!// F(orce)t(tangential)D(amping)2(two)d(imension)_M(in)d(lin)
Real Function FtD2d_Md( i,j,G1,G2,Miv1,Miv2,VtRel,DeltaN )
	Implicit None
	!// Damping constant
	Real( kind=MyRk ),Parameter ::		Beta=0.3033156 !Beta = log(2.71829) / sqrt( (log(2.71828))**2+Pi**2 )
	Integer(kind=MyIk) ::		i !// Particle 1
	Integer(kind=MyIk) ::		j !// Particle 2
	Real( kind=MyRk ) ::		VtRel !// Relative velocity
	Real( kind=MyRk ) ::		S !// Tangential stiffness
	Real( kind=MyRk ) ::		M !// Equivalent mass
	Real( kind=MyRk ) ::		M1 !// Mass of particle 1
	Real( kind=MyRk ) ::		M2 !// Mass of particle 2
	Real( kind=MyRk ) ::		G1 !// Shear modulus of particle 1
	Real( kind=MyRk ) ::		G2 !// Shear modulus of particle 2
	Real( kind=MyRk ) ::		R1 !// Radius of particle 1
	Real( kind=MyRk ) ::		R2 !// Radius of particle 2
	Real( kind=MyRk ) ::		Miv1 !// Poisson ratio of particle 1
	Real( kind=MyRk ) ::		Miv2 !// Poisson ratio of particle 2
	Real( kind=MyRk ) ::		Requiv !// Equivalent Radius
	Real( kind=MyRk ) ::		Gequiv !// Equivalent shear modulus
	Real( kind=MyRk ) ::		DeltaN !// Normal compression of two particles

	R1 = Pt(i)%Radius
	R2 = Pt(j)%Radius
	M1 = Pt(i)%Mass
	M2 = Pt(j)%Mass
	M = M1*M2/(M1+M2)
	Requiv = R1*R2/(R1+R2)
	Gequiv = (2.0-Miv1)/G1 + (2.0-Miv2)/G2 
	S = 8.0 * Gequiv * sqrt( Requiv*DeltaN )
	FtD2d_Md = -2.0*sqrt(5.0/6.0)*Beta*sqrt(S*M)*Abs(VtRel)

End Function FtD2d_Md



!// Orbit of sand 
Subroutine Orbit( Sw_W,Sw_E )
	Implicit None
	Logical, Intent(In) ::		Sw_W !// Switch of weather couping wind and sand
	Logical, Intent(In) ::		Sw_E !// Switch of weather consider electric field force of sand movement
	Integer ::					i, j
	Real( kind=MyRk ) ::		Lij !// Distance of two sands
	Real( kind=MyRk ) ::		Rt !// Tangental compression of two sands
	Real( kind=MyRk ) ::		VnRel !// Relative speed of two particles at normal direction
	Real( kind=MyRk ) ::		VtRel !// Relative speed of two particles at tangental direction


	Integer( kind=MyIk ) :: 	i1, k, H, L, M, N, Zi
	Real( kind=MyRk ) :: 		Fx, Fy, Fn, Ft, Fnx, Fny, Ftx, Fty, Dn, Dt, Dnx, Dny, Dtx, Dty, Fw
	Real( kind=MyRk ) :: 		F_DragX, F_DragZ
	Real( kind=MyRk ) :: 		CosAngle, SinAngle, ViNt, VjNt, ViNn, VjNn
	Real( kind=MyRk ) ::		Requiv, E1, E2, KnEdge, Compression




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
		Do H = Pt_old(i)%M-1, Pt_old(i)%M+1!Loop Coloum of i particle's neighbor mesh grid
			Do L = Pt_old(i)%N-1, Pt_old(i)%N+1!Loop Row of i particle's neighbor mesh grid 	
				k = Mesh( 0,L,H )!Particle number in the mesh (H,L)
				If( k<1 ) Cycle
				!Calculation of compression of i and j
				Do i1 = 1, k!Loop particles in neighbor mesh grid 
					j = Mesh( i1,L,H )
					If( j==i .or. j>N_Total) Cycle
					If(H==0)then
						If( Pt_Old(i)%Radius>Pt_old(i)%Y )then
							Requiv = Pt_old(i)%Radius
							KnEdge = 4.0*Ym*Requiv**0.5/3.0
							Fny = Fny + KnEdge * ( Pt_Old(i)%Radius-Pt_Old(i)%Y )**1.5
						End If
						Cycle
					Elseif(H==N_RowOfMesh)then
						If( Height<Pt_old(i)%Y+Pt_old(i)%Radius )then
							!Requiv = Pt_old(i)%Radius
							!KnEdge = 4.0*Ym*Requiv**0.5/3.0
							!Fny = Fny - KnEdge * ( Pt_old(i)%Y+Pt_old(i)%Radius-Height )
							Fny = 0.0 + Fny
						End If
						Cycle
					Else
						If( L==0 )then
							Pt_old(j)%X = Pt_old(j)%X - Length
						Elseif( L==N_ColumnOfMesh+1 )then
							Pt_old(j)%X = Pt_old(j)%X + Length
						End If	
					End If
					Lij = Distance( x1=Pt_old(i)%X,x2=Pt_old(j)%X, &
									y1=Pt_old(i)%Y,y2=Pt_old(j)%Y )
					Compression = Pt_Old(i)%Radius+Pt_Old(j)%Radius - Lij
					If( Compression>0.0 )then!Particle i and j contact and Compressioned 
						!Calculation of Normal Force
						Fn = FnS2d_Hz( i=i,j=j,DeltaN=Compression,E1=Ym,E2=Ym,&
										Miv1=Miv,Miv2=Miv )
						CosAngle = ( Pt_Old(i)%X-Pt_Old(j)%X ) / Lij
						SinAngle = ( Pt_Old(i)%Y-Pt_Old(j)%Y ) / Lij
						Fnx = Fn * CosAngle + Fnx
						Fny = Fn * SinAngle + Fny
					
						!Calculation of Tangential Force
						!Normal vector Nn = ( Pt_Old(i)%X-Pt_Old(j)%X, Pt_Old(i)%Y-Pt_Old(j)%Y )
						!Tangential vector  Nt = ( Pt_Old(j)%Y-Pt_Old(i)%Y, Pt_Old(i)%x-Pt_Old(j)%x )
						!Velocity vector Vi = ( Pt_Old(i)%Vx, Pt_Old(i)%Vy ) Vj = ( Pt_Old(j)%Vx, Pt_Old(j)%Vy )
						ViNt = Vector_Projection( x1=Pt_old(i)%Vx, y1=Pt_old(i)%Vy, &
												  x2=Pt_old(j)%Y-Pt_old(i)%Y, y2=Pt_old(i)%x-Pt_old(j)%x )
						VjNt = Vector_Projection( x1=Pt_old(j)%Vx, y1=Pt_old(j)%Vy, &
												  x2=Pt_old(j)%Y-Pt_old(i)%Y, y2=Pt_old(i)%x-Pt_old(j)%x )
						!ViNt = ViNt + Pt(i)%W*Pt(i)%Radius
						!VjNt = VjNt + Pt(j)%W*Pt(j)%Radius
						VtRel = ( VjNt+Pt(j)%W*Pt(j)%Radius ) - ( ViNt+Pt(i)%W*Pt(i)%Radius ) 
						Rt = Compression*0.5
						Ft = Fts2d_Md( i=i,j=j,Miv1=Miv,Miv2=Miv,E1=Ym,E2=Ym,&
										DeltaN=Compression,DeltaT=Rt )
						If( VtRel/=0.0 )then
							Ft = Abs(VtRel)/VtRel*Ft 
							Ftx = Ftx + Vector_Projection( &
								x1 = Ft*(Pt_old(j)%Y-Pt_old(i)%Y)/Lij,&
								y1 = Ft*(Pt_old(i)%x-Pt_old(j)%X)/Lij,&
								x2 = 1.0,&
								y2 = 0.0 )
							Fty = Fty + Vector_Projection( &
								x1 = Ft*(Pt_old(j)%Y-Pt_old(i)%Y)/Lij,&
								y1 = Ft*(Pt_old(i)%x-Pt_old(j)%X)/Lij,&
								x2 = 0.0,&
								y2 = 1.0 )
						Else
							Ft = 0.0
						End If
						!Calculatin of damping force( Tangential and normal direction )
						ViNn = Vector_Projection( x1=Pt_old(i)%Vx, y1=Pt_old(i)%Vy, &
												  x2=Pt_old(i)%X-Pt_old(j)%X, y2=Pt_old(i)%Y-Pt_old(j)%Y )
						VjNn = Vector_Projection( x1=Pt_old(j)%Vx, y1=Pt_old(j)%Vy, &
												  x2=Pt_old(i)%X-Pt_old(j)%X, y2=Pt_old(i)%Y-Pt_old(j)%Y )
						VnRel = VjNn -ViNn
						Dn = FnD2d_Md( i=i,j=j,E1=Ym,E2=Ym,Miv1=Miv,&
							Miv2=Miv,DeltaN=Compression,VnRel=VnRel )
						Dt = FtD2d_Md( i=i,j=j,Miv1=Miv,Miv2=Miv,G1=Ym/2.0/(1.0+Miv),&
							G2=Ym/2.0/(1.0+Miv),DeltaN=Compression,VtRel=VtRel )
						Dnx = Dn * CosAngle + Dnx
						Dny = Dn * SinAngle + Dny
						If( VtRel/=0.0 )then
							Dtx = Dtx + Vector_Projection( &
								x1 = Dt*(Pt_old(j)%Y-Pt_old(i)%Y)/Lij,&
								y1 = Dt*(Pt_old(i)%x-Pt_old(j)%X)/Lij,&
								x2 = 1.0,&
								y2 = 0.0 )
							Dty = Dty + Vector_Projection( &
								x1 = Dt*(Pt_old(j)%Y-Pt_old(i)%Y)/Lij,&
								y1 = Dt*(Pt_old(i)%x-Pt_old(j)%X)/Lij,&
								x2 = 0.0,&
								y2 = 1.0 )
						End If
						!// Transshipment angular velocity
						Fw = Fw + Ft*Pt(i)%Radius
						!Calculatin of electric between i and j
						If( Compression>ComMax(i,j) )then
							ComMax( i,j ) = Compression
						Else 
							Pt(i)%Q = Pt(i)%Q + DeltaQ( i,j )
							ComMax( i,j ) = 0.0
						End If
					Else
						ComMax( i,j ) = 0.0
					End If!Particle i and j contact and Compressioned 
					If( L==0 )then
						Pt_old(j)%X = Pt_old(j)%X + Length
					Elseif( L==N_ColumnOfMesh+1 )then
						Pt_old(j)%X = Pt_old(j)%X - Length
					End If	
				End Do!Loop particles in neighbor mesh grid 
			End Do !Loop column of i particle's neighbor mesh grid
		End Do!Loop Row of i particle's neighbor mesh grid
		!Calculation of drag force
		If( Sw_W .and. Pt_Old(i)%Y-Z0>H_SandBed )then
			F_DragX = DragX( i )
			F_DragZ = DragZ( i )
		Else
			F_DragX = 0.0
			F_DragZ = 0.0 
		End If

		!Calculation of movement
		Fx = Ftx + Fnx  + Dnx + Dtx + F_DragX
		Fy = Fty + Fny  + Dny + Dty + F_DragZ
		Pt(i)%Vx = Pt(i)%Vx + Fx/Pt(i)%Mass * Deltat
		Pt(i)%Vy = Pt(i)%Vy + Fy/Pt(i)%Mass * Deltat + G*Deltat
		Pt(i)%V = Distance( x1=Pt(i)%Vx, y1=Pt(i)%Vy, x2=0.0, y2=0.0 )
		Pt(i)%X = Pt(i)%X + Pt(i)%Vx * DeltaT + 0.5*Fx/Pt(i)%Mass*DeltaT**2
		Pt(i)%y = Pt(i)%Y + Pt(i)%Vy * DeltaT + 0.5*( Fy/Pt(i)%Mass+G )*DeltaT**2
		Pt(i)%W = Pt(i)%W + Fw*DeltaT/Pt(i)%Inertia
		If( Pt(i)%Y<Pt(i)%Radius )then
			Pt(i)%Y = Pt(i)%Radius
			Pt(i)%Vx = 0.9*Pt(i)%Vx
			Pt(i)%Vy = 0.9*Pt(i)%Vy
			Pt(i)%V = Distance( x1=Pt(i)%Vx , y1=Pt(i)%Vy, x2=0.0, y2=0.0 )
		End If
		!Calculation of TauP
		If( Sw_W .and. Pt(i)%Y>H_SandBed )then
			Call Get_TauP( i )
			!If( Pt(i)%Y+Pt(i)%Radius<H_Qm )then
			!	zi = Int( (Pt(i)%Y-H_SandBed)/CellOfQm ) + 1
			!	Qm( zi ) = Qm( zi ) + Pt(i)%Mass!/Abs(Pt(i)%Vy)/DeltaT
			!End If
		End If
		!// Control of sand charge 
		If( Pt(i)%Y+Pt(i)%Radius<H_SandBed )then
			Pt(i)%Q = 0.0
			Pt(i)%Co = 0
			Pt(i)%Alpha = 0.132
			Pt(i)%Pd = Pdd
		End If

	End DO 


End Subroutine Orbit


!// Distance of two particles
Real Function Distance( x1,y1,x2,y2 )
	Implicit None
	Real( kind=MyRk ), Intent(In) ::		x1,y1,x2,y2

	Distance = ( (x1-x2)**2 + (y1-y2)**2 )**0.5

End Function Distance


!// Projection of vector (x1,y1) to vector (x2,y2)
Real Function Vector_Projection( x1,y1,x2,y2 )
	Implicit None
	Real( kind=MyRk ), Intent(In) ::		x1,y1,x2,y2

	Vector_Projection = ( x1*x2+y1*y2 )/sqrt( x2**2+y2**2 )

End Function Vector_Projection


!// Angle of vector (x1,y1) and vector (x2,y2)
Real Function Vector_angle( x1,y1,x2,y2 )
	Implicit None
	Real( kind=MyRk ), Intent(In) ::		x1,y1,x2,y2

	Vector_angle = ( x1*x2+y1*y2 )/&
		sqrt( x2**2+y2**2 ) / sqrt( x1**2+y1**2 )

End Function Vector_angle



End Module Collision
