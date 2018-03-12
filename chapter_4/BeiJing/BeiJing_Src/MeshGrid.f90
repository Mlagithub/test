Module MeshGrid
	Implicit None
	Integer :: 							N_ColumnOfMesh ! Column of mesh 
	Integer :: 							N_RowOfMesh ! Row of mesh
	Integer,Allocatable ::				Mesh(:,:,:) ! Mesh information 
	Real( kind=4 ) ::					MeshSize ! Mesh size


Contains


!// Mesh the two-dimensional area by the given size
Subroutine Get_Mesh_2D( N,Boder_Left,Boder_Right,Boder_lower,Boder_upper,Rmin,MeshSize )
	Implicit None
	Integer ::									i, j
	Integer, Intent(In) ::						N
	Real( kind=4 ), Parameter ::				PI=3.14159
	Real( kind=4 ), Intent(In) ::				Rmin
	Real( kind=4 ), Intent(Inout)  ::			MeshSize ! Mesh size
	Real( kind=4 ), Intent(In)  :: 				Boder_Left ! The left border of the area
	Real( kind=4 ), Intent(In)  :: 				Boder_Right ! The right border of the area
	Real( kind=4 ), Intent(In)  :: 				Boder_lower ! The lower border of the area
	Real( kind=4 ), Intent(In)  :: 				Boder_upper ! The upper border of the area

	If( MeshSize<0.0 )then
		Print*,'MeshSize is not corrent'
		Print*,'MeshSize: ', MeshSize
		Stop
	End If
	!// Row and column number of mesh grid
	N_ColumnOfMesh = Int( (Boder_Right-Boder_Left)/MeshSize )
	MeshSize = (Boder_Right-Boder_Left)/Real(N_ColumnOfMesh)
	N_RowOfMesh = Int( (Boder_Upper-Boder_lower)/MeshSize )
	i = Int( MeshSize**2/Pi/Rmin**2 ) + 1
	Allocate( Mesh(0:i,0:N_ColumnOfMesh+1,0:N_RowOfMesh+1) )

	!// Initialize the mesh array 
	!// Set the out grid of area to a specific particle
	Mesh(:,:,:) = 0
	!// Left grid 
	Mesh(1,0,:) = N+1
	Mesh(0,0,:) = 1
	!// Right grid 
	Mesh(1,N_ColumnOfMesh+1,:) = N+3
	Mesh(0,N_ColumnOfMesh+1,:) = 1
	!// Lower grid 
	Mesh(1,:,0) = N+2
	Mesh(0,:,0) = 1
	!// Upper grid 
	Mesh(1,:,N_RowOfMesh+1) = N+4
	Mesh(0,:,N_RowOfMesh+1) = 1

!	Print*,'Number of column and row of grid :', N_ColumnOfMesh, N_RowofMesh

End Subroutine Get_Mesh_2D


!// Update particle mesh information
Subroutine MeshUpdate2D( Boder_Left,Boder_lower )
	Use GlobalVariable
	Implicit None
	Integer ::						i, M, N, j
	Real( kind=4 ) ::				X, Y
	Real( kind=4 ) :: 				Boder_Left ! The left border of the area
	Real( kind=4 ) :: 				Boder_lower ! The lower border of the area
	
	Mesh(:,:,:) = 0
	Do i = 1, N_Total
		X = Pt(i)%X
		Y = Pt(i)%Y
		M = Floor((Y-Boder_lower)/MeshSize) + 1
		N = Floor((X-Boder_left)/MeshSize) + 1
		If( N<1 )then
			N = N_ColumnOfMesh
			Pt(i)%x = Pt(i)%x + Length
		End If
		If( N>N_ColumnOfMesh )then
			N = 1
			Pt(i)%x = Pt(i)%x - Length
		End If
		If( M>N_RowofMesh )then
			M = N_RowofMesh
			Pt(i)%vx = 1.0E-3!Pt(i)%vx
			Pt(i)%vy = -1.0E-3!-Pt(i)%vy
		End If
		If( M<1 )then
			M = 1
			Pt(i)%Y = Pt(i)%Radius
			Pt(i)%vx = Pt(i)%Vx
			Pt(i)%vy = -Pt(i)%Vy
		End If
		Pt(i)%M = M
		Pt(i)%N = N
		j = Mesh(0,n,m)
		Mesh(j+1,n,m) = Pt(i)%i
		Mesh(0,n,m) = Mesh(0,n,m) + 1
	End Do 

End Subroutine MeshUpdate2D


!// Set periodic boundary condition
Subroutine PeriodicBoundary2D( Boder )
	Implicit None
	Integer ::							i
	Logical , Intent(In) ::				Boder(2)

	If( Boder(1) )then
		Mesh(:,0,:) = Mesh(:,N_ColumnOfMesh,:)
		Mesh(:,N_ColumnOfMesh+1,:) = Mesh(:,1,:)
	End If
	If( Boder(2) )then
		Mesh(:,:,0) = Mesh(:,:,N_RowOfMesh)
		Mesh(:,:,N_RowOfMesh+1) = Mesh(:,:,0)
	End If

End Subroutine PeriodicBoundary2D




End Module MeshGrid
