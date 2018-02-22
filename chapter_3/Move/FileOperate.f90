!// This module contains all functions of file operate and screen output 
Module FileOperate
	Use GlobalVariable
	Implicit None

Contains 

!// Description of the project  
Subroutine Description()
	Implicit None
	Integer ::								i
	Character(len=100) ::					A, B, C(20)

	A = '|-------------------------------------------------------------------------------------|'
	C(1) = '| '
	C(2) = '|  Particle parking simulation. '
	C(3) = '| '
	C(4) = '|  Tightly stacked balls in a 2D rectangular. Similar to Y.T.Feng(2003)  '
	C(5) = '|  but they do not consider structural balance and compression of the '
	C(6) = '|  particles. This algorithm use a new method which simulated gravity '
	C(7) = '|  of the structure '
	C(8) = '| '
	C(9) = '|  Author:	MU LiAn LanZhou University '
	C(10) = '|  Date:	2017/07/06 '
	C(11) = '|  Mail:	1791430496@qq.com '
	C(12) = '|  Edition:	1.0 '
	Write(*,*)
	Write(*,*)
	Write(*,"(A100)") Adjustl(A)
	Do i = 2, 12
		B = C(i)
		Write(*,"(A100)") Adjustl(B)
	End Do 
	Write(*,"(A100)") Adjustl(A)
	Write(*,*)
	Write(*,*)

End Subroutine Description


!// Terminal output
Subroutine ShowInformation( C,N,FormatLabel )
	Implicit None
	Integer ::								i, N, FormatLabel
	Character(len=100) ::					A, B
	Character(len=100) ::					C(:)

	A = '|-------------------------------------------------------------------------------------|'
	Select Case( FormatLabel )
	Case( 0 )
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
	Case( 1 )
		Write(*,"(A100)") Adjustl(A)
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
	Case( 2 )
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
		Write(*,"(A100)") Adjustl(A)
	Case( 3 )
		Write(*,"(A100)") Adjustl(A)
		Do i = 1, N
			Write(*,"(A100)") Adjustl(C(i))
		End Do 
		Write(*,"(A100)") Adjustl(A)
	End Select 

End Subroutine ShowInformation

	
!// Write result to specific file named FileName
Subroutine SaveResult( FileName,N )
	Implicit None
	Integer ::						i, N
	Character ::					FileName*512

	FileName = Trim(Adjustl(FileName))
	Open( 99,File=FileName )
	Do i = 1, N
		Write( 99,'(i13,2x,3(E13.6,2x))' )  i,Pt(i)%X,Pt(i)%Y,Pt(i)%Radius
	End Do 
	Close( 99 )

End Subroutine SaveResult


!// Close specific file
Subroutine CloseFile( Filepointer )
	Integer ::			Filepointer

	Close( Filepointer )

End Subroutine


!// Delete specific file
Subroutine DeleteFile( Filepointer )
	Integer ::			Filepointer

	Close( Filepointer,Status='Delete' )

End Subroutine DeleteFile


!// Write program parameter into configuration file
Subroutine SaveConfigure()
	Implicit None
	Integer ::				j
	Character ::						WFormat*30

	Open( 11, File='Configure' )
	Write( 11,"(E10.3,2x,A30)" ) Accuracy, 'Newton ineration Accuracy'
	Write( 11,"(i10,2x,A30)" ) N_Total, 'Particle N_Total'
	Write( 11,"(i10,2x,A30)" ) N_Sorts, 'Particle sorts'
	Write( 11,'(E10.3,2x,A30)' ) RHO,'Density of particle'
	Write( 11,'(E10.3,2x,A30)' ) Ym,'Young Modulus'
	Write( 11,'(E10.3,2x,A30)' ) Miv,'Poisson ratio'
	Write( 11,"(E10.3,2x,A30)" ) Bb,'Bottom edge'
	Write( 11,'(E10.3,2x,A30)' ) Bt,'Top edge'
	Write( 11,'(E10.3,2x,A30)' ) Bl,'Left edge'
	Write( 11,'(E10.3,2x,A30)' ) Br,'Right edge'
	If( N_Sorts==N_Total )then
		Write( 11,* ) 'Too many Radius'
		Write( 11,* ) 'To many Ratio'
	Else
		Write( WFormat,* ) N_Sorts
		WFormat = '('//Trim(Adjustl(WFormat))//'(E10.3),2x,A30)'
		Write( 11,WFormat ) Radius,'Radius'
		Write( 11,WFormat ) Ratio, 'Ratio'
	End If
    Do j = 1, N_Total
		Write(11, "(I9,3x 2(E13.6))") Pt(j)%i, Pt(j)%Radius
    End Do
	Close(11)

End Subroutine SaveConfigure


End Module FileOperate
