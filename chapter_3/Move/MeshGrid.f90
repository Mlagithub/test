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
	Integer ::						i, j, N
	Real( kind=4 ),Parameter ::		PI=3.14159
	Real( kind=4 ) ::				Rmin
	Real( kind=4 ) ::				MeshSize ! Mesh size
	Real( kind=4 ) :: 				Boder_Left ! The left border of the area
	Real( kind=4 ) :: 				Boder_Right ! The right border of the area
	Real( kind=4 ) :: 				Boder_lower ! The lower border of the area
	Real( kind=4 ) :: 				Boder_upper ! The upper border of the area

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

End Subroutine Get_Mesh_2D


!// Update particle mesh information
Subroutine MeshUpdate( i,X,Y,M,N,Boder_Left,Boder_lower )
	Implicit None
	Integer ::						i, M, N, j
	Real( kind=4 ) ::				X, Y
	Real( kind=4 ) :: 				Boder_Left ! The left border of the area
	Real( kind=4 ) :: 				Boder_lower ! The lower border of the area
	
	M = Floor( (Y-Boder_lower)/MeshSize ) + 1
	N = Floor( (X-Boder_Left)/MeshSize ) + 1
	If( M>N_RowOfMesh )then
		Print*,'MeshGrid:',i,M
		stop
	End If
	If( N>N_ColumnOfMesh )then
		Print*,'MeshGrid:',i,N
		stop
	End If
	j = Mesh(0,n,m)
	Mesh(j+1,n,m) = i
	Mesh(0,n,m) = Mesh(0,n,m) + 1

End Subroutine MeshUpdate


End Module MeshGrid
