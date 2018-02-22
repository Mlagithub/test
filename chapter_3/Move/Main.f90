Program Main
	Use Accumulation
	Use Sedimentation
	Use FileOperate
	Implicit None
	Character( len=100 ) ::				A(100)
	Integer ::							i, Err



	Call Description()

	A(1) = '|    Please select the label of the following module'
	A(2) = '|'
	A(3) = '|        1>    Particle accumulation'
	A(4) = '|        2>    Particle sedimentation'
	A(5) = '|        3>    Particle accumulation and sedimentation'
	Call ShowInformation( A(1:5),5,3 )
	Do
		Read( *,*,Iostat=Err ) i
		If( Err/=0 .or. i<0 .or. i>3 )then
			Print*,'Try again'
		Else
			Exit
		End If
	End Do 


	Select Case( i )
	Case( 1 )
		Call ParticleFilling()
	Case( 2 )
		Call ParticleSettlement()
	Case( 3 )
		Call ParticleFilling()
		Call ParticleSettlement()
	End Select 



End Program Main

