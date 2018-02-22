Module ScheduleBar
	Implicit None
	private
	Logical, parameter, public ::	CMD_PROGRESS_ABSOLUTE = .True.
	Type, public ::					CLS_CMD_Progress
		Integer , private ::		N, lens, i
		Character ::				M = "*", O = "."
		Character(len=64) ::		Prefix
	Contains
		Procedure ::				Set
		Procedure ::				Put
	End Type CLS_CMD_Progress

contains

Subroutine Set( this, N, L )
	Class( CLS_CMD_Progress ) ::	this
	Integer , Intent( IN )	::		N , L
	this % N    = N
	this % lens = L
	this % i = 0
	this % Prefix = " Progress: " !//
End Subroutine Set


Subroutine Put( this , k, bAbsol )
	Class( CLS_CMD_Progress ) ::		this
	Integer , Intent( IN ) ::			k
	Logical , optional ::				bAbsol
	Character( len=1 ) ::				br
	Integer ::							jm

	this%i = this%i + k 
	If( present(bAbsol) )then
		if( bAbsol ) this%i = k
	End If
	If( this%i>this%n ) this%i = this%n
	jm = Nint( Real(this%i*this%lens)/Real(this%n) )
	If( jm<1 ) jm=1
	If( this%i<this%n )then
		br = char(13)
	Else
		br = char(10)
	End If
	Write( *,"(5A,F6.2,2A\)" ) Trim(this%Prefix),'[', &
		Repeat( this%m,jm ), Repeat( this%o,this%lens-jm ), ']', &
		this%i*100.0/this%n, "%", br
End Subroutine Put


End Module ScheduleBar


