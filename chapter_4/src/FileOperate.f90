!// This module contains all functions of file operate and screen output 
Module FileOperate
	Use GlobalVariable
	Use Electricity
	Implicit None

Contains 

!// Description of the project  
Subroutine Description()
	Implicit None
	Integer ::								i
	Character(len=100) ::					A, B, C(20)

	A = '|-------------------------------------------------------------------------------------|'
	C(1) = '| '
	C(2) = '|  Wind sand salation simulaion. '
	C(3) = '| '
	C(4) = '|  Tightly stacked balls in a 2D rectangular. Similar to Y.T.Feng(2003)  '
	C(5) = '|  but they do not consider structural balance and compression of the '
	C(6) = '|  particles. This algorithm use a new method which simulated gravity '
	C(7) = '|  of the structure '
	C(8) = '| '
	C(9) = '|  Author:	MU LiAn LanZhou University '
	C(10) = '|  Date:	2017/07/26 '
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
	Character(len=512) ::					A, B
	Character(len=*) ::						C(:)

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
Subroutine SaveResult( FileFlag,FileName,N )
	Implicit None
	Integer ::						i, N
	Integer, optional ::			FileFlag
	Character(len=*), optional ::	FileName
	Character(len=512) ::			FN

	If( present(FileName) )then
		FN = FileName
		FN = Trim(Adjustl(FN))
	Else
		If( present(FileFlag) )then
			Write( FN,* ) FileFlag
			FN = Trim(Adjustl(FN))//'.txt'
		Else
			Print*,'Input at least one of FileFlag or FileName'
			Print*,'Please input FileName'
			Read( *,* ) FN
			FN = Trim(Adjustl(FN))
		End If
	End If

	Open( 99,File=FN )
	Do i = 1, N
		Write( 99,'(i13,2x,6(E13.6,2x))' ) i,Pt(i)%X,Pt(i)%Y,Pt(i)%Radius,Pt(i)%Vx,Pt(i)%Vy,Pt(i)%Q
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
	Character ::			WFormat*30

	Open( 11,File='Configure' )
	Write( 11,"(i10,2x,A30)" ) N_Total, 'Total number of sands'
	Write( 11,"(i10,2x,A30)" ) N_Incident, 'Number of incident sand'
	Write( 11,"(i10,2x,A30)" ) N_Sorts, 'Number of radius'
	Write( 11,'(E10.3,2x,A30)' ) RH,'Relativy humidity'
	Write( 11,'(E10.3,2x,A30)' ) Temperature,' Global Temperature'
	Write( 11,'(E10.3,2x,A30)' ) Xfn,'Xfn'
	Write( 11,"(E10.3,2x,A30)" ) EdgeRight, 'Right boder of sand bed'
	Write( 11,"(E10.3,2x,A30)" ) EdgeTop, 'Upper boder of sand bed'
	If( N_Sorts==N_Total )then
		Write( WFormat,* ) 2
		WFormat = '('//Trim(Adjustl(WFormat))//'(E10.3),2x,A30)'
		Write( 11,WFormat ) MaxRadius,MinRadius,'Max and Min Radius'
	Else
		Write( WFormat,* ) N_Sorts
		WFormat = '('//Trim(Adjustl(WFormat))//'(E10.3),2x,A30)'
		Write( 11,WFormat ) Radius,'Radius'
	End If
	Close(11)

End Subroutine SaveConfigure



!// Read program parameter form configuration file
Subroutine ReadConfigure()
	Implicit None
	Integer ::		Err

	Open( 11,File='Configure',Status='old',Iostat=Err  )
	If( Err/=0 )then
		Print*,'Open Confiugre file fail'
		Stop
	End If
	
	Read( 11,*,Iostat=Err ) N_Total 
	If( Err/=0 .or. N_Total<=0 )then
		Print*,'Read Configure fail'
		Print*,'Total number of sand is:', N_Total
		Stop
	End If

	Read( 11,*,Iostat=Err ) N_Incident
	If( Err/=0 .or. N_Incident<0 .or. N_Incident>N_Total )then
		Print*,'Read Configure fail'
		Print*,'Number of inciden sand is:', N_Incident
		Stop
	End If

	Read( 11,*,Iostat=Err ) N_Sorts
	If( Err/=0 .or. N_Sorts<=0 )then
		Print*,'Read Configure fail'
		Print*,'Number of radius is:', N_Sorts
		Stop
	End If

	Read( 11,*,Iostat=Err ) Rh
	If( Err/=0 .or. Rh<=0.0 )then
		Print*,'Read Configure fail'
		Print*,'Relative humidity is:', Rh
		Stop
	End If

	Read( 11,*,Iostat=Err ) Temperature
	If( Err/=0 .or. Temperature<=0.0 )then
		Print*,'Read Configure fail'
		Print*,'Temperature is:', Temperature
		Stop
	End If

	Read( 11,*,Iostat=Err ) Xfn
	If( Err/=0 )then
		Print*,'Read Configure fail'
		Print*,'Xfn is:', Xfn
		Stop
	End If

	Read( 11,*,Iostat=Err ) EdgeRight
	If( Err/=0 .or. EdgeRight<0.0 )then
		Print*,'Read Configure fail'
		Print*,'Right boder of sand bed is:',EdgeRight
		Stop
	End If

	Read( 11,*,Iostat=Err ) EdgeTop
	If( Err/=0 .or. EdgeTop<0.0 )then
		Print*,'Read Configure fail'
		Print*,'Top boder of sand bed is:',EdgeTop
		Stop
	End If

	If( N_Sorts/=N_Total )then
		Allocate( Radius(N_Sorts) )
		Read( 11,*,Iostat=Err ) Radius(:)
		If( Err/=0 )then
			Print*,'Read Configure fail'
			Print*,'Radius', Radius
			Stop
		End If
	Else
		Allocate( Radius(2) )
		Read( 11,*,Iostat=Err ) Radius(:)
		If( Err/=0 )then
			Print*,'Read Configure fail'
			Print*,'Radius', Radius
			Stop
		End If
	End If

	Close(11)

End Subroutine ReadConfigure



!//Get rows of file
Subroutine Getfilen( FileName,Hs )
	Implicit None
	Integer( kind=MyIK ) , Intent(Out) :: 	Hs
	Integer( kind=MyIK ) ::					Err
	Character(*) ::							FileName
	Character ::							A*1

	!FileName = Trim(Adjustl(FileName))
    Open( 11,File=FileName,Status='old',Iostat=Err )
	If( Err/=0 )then
		Print*,'Can not find file'
		Read(*,*)
		Stop
	End If
	Hs = 0
	Rewind(11)
	Do 
		Read( 11, *, Iostat=Err ) A
		If(Err /=  0) Exit
		Hs = Hs + 1
	End Do 
	Rewind( 11 )
	Close( 11 )

End Subroutine Getfilen



End Module FileOperate
