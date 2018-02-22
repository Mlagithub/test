Program Main
    Implicit None
    Integer :: i, N
    Real(kind=8) :: Y, X
    
    Open(11, File="Log.txt", Access="Append")
    N = 1E6
    
    Do i = 1, N
        X = Real(i)
        Y = - Log(X) / Log(Real(N)) + 1
        Write(11,"(i10, 2x, F13.6, 2x, F13.6)") i, X, Y
    End Do 
    
End Program Main
        