Program Main
    USE OMP_LIB
    implicit none
    character(*), parameter:: InputFile='input.txt',OutputFile='data.plt' ! Имя input и output файла
    character MeshFile*30 ! Имя файла с сеткой
    integer, parameter:: IO = 12 ! input-output номер
    integer :: i, j, Ni, Nj, grad_method, div_sheme
    real,allocatable,dimension(:,:):: X,Y,P,CellVolume ! Скалярные массивы
    real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! Массив геометрии

    real,allocatable,dimension(:,:,:):: GradP, GradP_Exact, GradP_Error

    real,allocatable,dimension(:,:):: DivV, DivV_Exact, DivV_Error ! Дивергенция скорости
    real,allocatable,dimension(:,:,:):: V ! Поле скорости


    real,allocatable,dimension(:,:):: LaplP, LaplP_Exact, LaplP_Error ! Лапласиан скаляра

    real,allocatable,dimension(:,:):: RotV, RotV_Exact, RotV_Error ! 2D ротор скорости


    real init_p, calc_divV_exact, calc_laplP_exact, calc_rotV_exact


    ! === For FLOS ====

    real Re, Pr, CFL, VNM, Vs, Ls, a_diff, t2, t1
    integer niter, scheme_res, iter
    real,allocatable,dimension(:,:,:):: VFlos, GradT ! V и T из Flos
    real,allocatable,dimension(:,:):: T, ResT, TFlos, T_Err, dtau ! Поля T









    !=== READ INPUT FILE ===
    WRITE(*,*) 'Чтение входного файла: ', InputFile
    OPEN(IO,FILE=InputFile)
    READ(IO,*) MeshFile ! сетка
    read(IO, *) grad_method ! Метод вычисления градиента
    read(IO, *) div_sheme ! Схема интерполяция
    CLOSE(IO)

    !=== READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
    WRITE(*,*) 'Чтение размеров сетки из файла: ', MeshFile
    OPEN(IO,FILE = MeshFile)
    READ(IO,*) NI,NJ
    WRITE(*,*) 'NI, NJ = ',NI,NJ

    !=== ALLOCATE ALL ARRAYS ===
    WRITE(*,*) 'Выделение памяти под массивы'
    allocate(X(NI,NJ)) ! Координаты X узлов сетки
    allocate(Y(NI,NJ)) ! Координаты Y узлов сетки
    allocate(P(0:NI,0:NJ)) ! Скалярная величина
    allocate(CellVolume(NI-1,NJ-1)) ! Объем ячеек
    allocate(CellCenter(0:NI,0:NJ,2)) ! Центры ячеек
    allocate(IFaceCenter( NI,NJ-1,2)) ! Центры граней для I-граней
    allocate(IFaceVector( NI,NJ-1,2)) ! Вектора нормали для I-граней
    allocate(JFaceCenter( NI-1,NJ,2)) ! Центры граней для J-граней
    allocate(JFaceVector( NI-1,NJ,2)) ! Вектора нормали для J-граней


    allocate(GradP(0:NI,0:NJ,2)) ! Градиент скаляра
    allocate(GradP_Exact(0:NI,0:NJ,2)) ! Точное значение GradP
    allocate(GradP_Error(0:NI,0:NJ,2)) ! Ошибка в вычислении GradP

    allocate(V(0:NI,0:NJ,2)) ! Поле скорости
    allocate(DivV(0:NI,0:NJ)) ! divV
    allocate(DivV_Exact(0:NI,0:NJ)) ! Точное значение divV
    allocate(DivV_Error(0:NI,0:NJ)) ! Ошибка в вычислении divV

    allocate(LaplP(0:NI,0:NJ)) ! Лапласиан скаляра
    allocate(LaplP_Exact(0:NI,0:NJ)) ! Точное значение LaplP
    allocate(LaplP_Error(0:NI,0:NJ)) ! Ошибка в вычислении LaplP


    allocate(RotV(0:NI,0:NJ)) ! Ротор скаляра
    allocate(RotV_Exact(0:NI,0:NJ)) ! Точное значение rotV
    allocate(RotV_Error(0:NI,0:NJ)) ! Ошибка в вычислении rotV


    !=== READ GRID ===
    WRITE(*,*) 'Чтение сетки из файла: ', MeshFile
    READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
    CLOSE(IO)
    P = 0.0

    !=== CALCULATE METRIC ===
    WRITE(*,*) 'Вычисление метрик'
    Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)


    !=== Initial Field ===
    write(*,*) 'Инициализация полей'
    do i = 0, NI
        do j = 0, NJ
            P(i, j) = init_p(CellCenter(i,j,1), CellCenter(i,j,2))
            call calc_gradP_exact(CellCenter(i,j,1), CellCenter(i,j,2), GradP_Exact(i,j,:))
            call init_V(CellCenter(i,j,1), CellCenter(i,j,2), V(i,j,:))
            DivV_Exact(i,j) = calc_divV_exact(CellCenter(i,j,1), CellCenter(i,j,2))
            LaplP_Exact(i,j) = calc_laplP_exact(CellCenter(i,j,1), CellCenter(i,j,2)) 
            RotV_Exact(I,J) = calc_rotV_exact(CellCenter(I,J,1),CellCenter(i,j,2))
        end do
    end do

    gradP = 0.0


    !=== Calculate Gradient ====
    write(*,*) 'Вычисление градиента скаляра'
    Call B_CalcGradient(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
    grad_method,P,GradP)

    GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)

    write(*,*) 'Max error GradP:', maxval(GradP_Error(1:NI-1,1:NJ-1,:))


    !=== Calculate Divergence ===

    Call B_CalcDiv(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
    V,DivV,P,GradP,div_sheme)
    DivV_Error = ABS((DivV_Exact-DivV)/DivV_Exact)

    write(*,*) 'Max error divV:', maxval(DivV_Error(1:NI-1,1:NJ-1))

    !=== Calculate Laplacian ===

    Call B_CalcLapl(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
    P,GradP,LaplP)
    LaplP_Error = ABS((LaplP_Exact-LaplP)/LaplP_Exact)
    write(*,*) 'Max error LaplP:', maxval(LaplP_Error(1:NI-1,1:NJ-1))

    
    !=== Calculate Rotor ===

    Call B_CalcRotor(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
    V,RotV)
    RotV_Error = ABS((RotV_Exact-RotV)/RotV_Exact)
    write(*,*) 'Max error rotV:', maxval(RotV_Error(1:NI-1,1:NJ-1))





    !=== OUTPUT FIELDS ===
    WRITE(*,*) 'Запись полей в файл: ', OutputFile
    Open(IO,FILE=OutputFile)
    Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP,GradP_Error, V, divV, DivV_Exact, DivV_Error,&
    LaplP,LaplP_Exact,LaplP_Error, RotV, RotV_Exact, RotV_Error)
    Close(IO)

    deallocate(X, Y)


    

    WRITE(*,*) 'Чтение входного файла для задачи к-д переноса '
    OPEN(IO,FILE='inputT.txt')
    READ(IO,*) MeshFile ! сетка
    READ(IO,*) scheme_res ! схема для resT
    READ(IO,*) Vs ! масштаб скорости
    READ(IO,*) Ls ! масштаб длины
    READ(IO,*) Re ! Re
    READ(IO,*) Pr ! Pr
    READ(IO,*) CFL ! CFL
    READ(IO,*) VNM ! VNM
    READ(IO,*) niter ! число итераций
    CLOSE(IO)

    
    allocate(X(NI,NJ)) ! Координаты X узлов сетки
    allocate(Y(NI,NJ)) ! Координаты Y узлов сетки
    allocate(GradT(0:NI,0:NJ,2)) ! GradT
    allocate(VFlos(0:NI,0:NJ,2)) ! Поле скорости из Flos
    allocate(T(0:NI,0:NJ)) ! T
    allocate(ResT(0:NI,0:NJ)) ! невязка T
    allocate(TFlos(0:NI,0:NJ)) ! T из Flos
    allocate(T_Err(0:NI,0:NJ)) ! ошибка T
    allocate(dtau(0:NI,0:NJ)) ! шаг по времени

    !=== READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
    WRITE(*,*) 'Чтение размеров сетки из файла: ', MeshFile
    OPEN(IO,FILE = MeshFile)
    READ(IO,*) NI,NJ
    WRITE(*,*) 'NI, NJ = ',NI,NJ
    WRITE(*,*) 'Чтение сетки из файла: ', MeshFile
    READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
    CLOSE(IO)

    !=== READ FLOS VELOCITY FIELD ===
    OPEN(IO,FILE = 'VelocityField.txt')
    READ(IO,*) ((VFlos(I,J,1),VFlos(I,J,2),I=0,NI),J=0,NJ)
    CLOSE(IO)
    !=== READ FLOS TEMPERATURE FIELD ===
    OPEN(IO,FILE = 'Temperature.txt')
    READ(IO,*) ((TFlos(I,J),I=0,NI),J=0,NJ)
    CLOSE(IO)









    !=== CALCULATE METRIC ===
    WRITE(*,*) 'Вычисление метрик'
    Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)

    !=== CALCULATE RESIDUAL T ===
    ! инициализация поля
    T(:,:) = 1.0
    ! стенки при постоянной T
    T(0,:) = 2.0 ! левая hot
    T(NI,:) = 1.0 ! правая cold
    dtau = 0.001
    a_diff = Vs*Ls/(Re*Pr) !коэффициент температуропроводности



 


    OPEN(IO,FILE = 'Residual.plt')
    write(io,*) 'Variables = "iterations", "Res T"'
    t1 = OMP_GET_WTIME()
    !$OMP PARALLEL
    PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()
    do iter = 1,niter
        ! вычисление GradT
        Call B_CalcGradient(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,&
        JFaceCenter,JFaceVector,grad_method,T,GradT)
        ! вычисление невязки для GradT
        Call B_CalcResT(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,&
        JFaceCenter,JFaceVector,Re,Pr,VFlos,T,GradT,ResT,scheme_res,&
        CFL,VNM,a_diff,dtau)

        !$OMP single
        write(io,*) iter, maxval(abs(ResT(1:NI-1,1:NJ-1)))
        !$OMP end single

        ! итерационно обновляем T c учетом Res и dtau
        !$OMP DO PRIVATE(I,J)
        do j = 1, NJ-1
            do i = 1, NI-1
                T(i,j) = T(i,j) - ResT(i,j)*dtau(i,j)
            end do
        end do
        !$OMP END DO
    end do
    !$OMP END PARALLEL
    t2 = OMP_GET_WTIME()
    close(io)

    !вычисление ошибки T
    T_Err = ABS((TFlos-T)/TFlos)
    print*, maxval( ABS((TFlos-T)/TFlos))
    print*, 'Время выполнения:', t2-t1, ' s'


    open(IO, file='T.plt')
    Write(IO,*) 'VARIABLES = "X", "Y", "T_Flos", "V_x", "V_y", "T", "Error(T)", "GradT_x", "GradT_y"'
    Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-10]=CELLCENTERED)'
    Write(IO,'(100F14.7)') X(1:NI,1:NJ)
    Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
    Write(IO,'(100F14.7)') TFlos(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') VFlos(1:NI-1,1:NJ-1,1)
    Write(IO,'(100F14.7)') VFlos(1:NI-1,1:NJ-1,2)
    Write(IO,'(100F14.7)') T(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') T_Err(1:NI-1,1:NJ-1)
    Write(IO,'(100F14.7)') GradT(1:NI-1,1:NJ-1,1)
    Write(IO,'(100F14.7)') GradT(1:NI-1,1:NJ-1,2)


    close(IO)







    END PROGRAM Main