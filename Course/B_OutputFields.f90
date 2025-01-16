subroutine B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradP_Error, V, divV, DivV_Exact, DivV_Error,&
    LaplP,LaplP_Exact,LaplP_Error, RotV, RotV_Exact, RotV_Error)

    real, dimension(NI,NJ):: X,Y
    real, dimension(0:NI,0:NJ):: P, divV, divV_Exact, DivV_Error
    real, dimension(0:NI,0:NJ,2) :: GradP, GradP_Error, V
    real, dimension(0:NI, 0:NJ) :: LaplP,LaplP_Exact,LaplP_Error
    real, dimension(0:NI,0:NJ):: rotV, rotV_Exact, RotV_Error




    Write(IO,*) 'VARIABLES = "X", "Y", "P", "GradP_X", "GradP_Y", "GradP_X_Err", "GradP_Y_Err",&
    "Vx", "Vy", "divV", "divV_Exact", "divV_Err", "LaplP", "LaplP_Exact", "LaplP_Err",&
    "rotV", "rotV_exact", "rotV_Err"'
    Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
    Write(IO,'(100F14.7)') X(1:NI,1:NJ)
    Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
    Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
    Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)
    Write(IO,'(100F14.7)') GradP_Error(1:NI-1,1:NJ-1,1)
    Write(IO,'(100F14.7)') GradP_Error(1:NI-1,1:NJ-1,2)
    Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
    Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)
    Write(IO,'(100F14.7)') DivV(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') DivV_Exact(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') DivV_Error(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') LaplP(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') LaplP_Exact(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') LaplP_Error(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') RotV_Exact(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') RotV_Error(1:NI-1,1:NJ-1)

End Subroutine