Subroutine B_CalcGradient(NI,NJ,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
    grad_method,P,GradP)
    implicit none
    integer Ni, Nj, k
    real RLinearInterp, init_p
    REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),& ! центры и объемы ячеек
     IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! центры и нормали граней для I-граней
     JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали граней для J-граней
     P(0:NI,0:NJ),& ! поле давления
     GradP(0:NI,0:NJ,2),GradP_tmp(0:NI,0:NJ,2) ! GradP
    integer grad_method !метод расчета
    real Sf(4,2),rf(4,2),Pf(4),GPE(2),rE(2)
    integer NeighCell(4,2),GrG_iter
    real d,dn, p_dum, PE
    integer i,j,inc,jnc,iface

    select case (grad_method)

    !метод Грина-Гаусса
    case(1)
        !$OMP DO private(i,j,Sf,rf,NeighCell,iface,inc,jnc,d,dn,Pf)
        do i = 1,NI-1
            do j = 1,NJ-1
            ! Вектор площади грани для каждой ячейки
            Sf(1,:) = -IFaceVector(i,j,:)
            Sf(2,:) = IFaceVector(i+1,j,:)
            Sf(3,:) = -JFaceVector(i,j,:)
            Sf(4,:) = JFaceVector(i,j+1,:)
            ! Координаты центров граней
            rf(1,:) = IFaceCenter(i,j,:)
            rf(2,:) = IFaceCenter(i+1,j,:)
            rf(3,:) = JFaceCenter(i,j,:)
            rf(4,:) = JFaceCenter(i,j+1,:)
            ! Индексы соседних ячеек
            NeighCell(1,:) = [i-1,j]
            NeighCell(2,:) = [i+1,j]
            NeighCell(3,:) = [i,j-1]
            NeighCell(4,:) = [i,j+1]
            ! Инициализация градиента
            GradP(i,j,:) = 0
            ! Прохождение по всем граням
            do iface = 1,4
                inc = NeighCell(iface,1)
                jnc = NeighCell(iface,2)
                d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! Расстояние от центра яч до грани
                dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! Расстояние до центра соседней яч
                Pf(iface) = RLinearInterp(d,dn,P(i,j),P(inc,jnc)) ! Интерполяция P на грань
                ! Расчет в приграничной ячейке
                if (dn < 1e-7) then
                    ! Значение в центре заграничной ячейке
                    p_dum = init_p(2*rf(iface,1) - CellCenter(i,j,1),2*rf(iface,2) - CellCenter(i,j,2))
                    Pf(iface) = 0.5*(p_dum + P(i,j))
                end if
                ! Обновляем GradP
                GradP(i,j,:) = GradP(i,j,:) + Pf(iface)*Sf(iface,:)
            end do
            GradP(i,j,:) = GradP(i,j,:)/CellVolume(i,j)
            end do
        end do
        !$OMP END DO



    
    !метод Грина-Гаусса с итерациями
    case(2)

        GrG_iter = 10 !кол-во итераций
        do k = 1,GrG_iter
            !$OMP DO private(i,j,Sf,rf,NeighCell,iface,inc,jnc,d,dn,PE,rE,GPE,Pf)
            do i = 1,NI-1
                do j = 1,NJ-1
                    ! Инициализация GradP
                    GradP_tmp(i,j,:) = 0

                    ! Вектор площади грани с учетом внешней нормали
                    Sf(1,:) = -IFaceVector(i,j,:)
                    Sf(2,:) = IFaceVector(i+1,j,:)
                    Sf(3,:) = -JFaceVector(i,j,:)
                    Sf(4,:) = JFaceVector(i,j+1,:)

                    ! Координаты центров граней
                    rf(1,:) = IFaceCenter(i,j,:)
                    rf(2,:) = IFaceCenter(i+1,j,:)
                    rf(3,:) = JFaceCenter(i,j,:)
                    rf(4,:) = JFaceCenter(i,j+1,:)

                    ! Индексы соседних ячеек
                    NeighCell(1,:) = [i-1,j]
                    NeighCell(2,:) = [i+1,j]
                    NeighCell(3,:) = [i,j-1]
                    NeighCell(4,:) = [i,j+1]

                    do iface = 1,4
                        inc = NeighCell(iface,1)
                        jnc = NeighCell(iface,2)
                        d = norm2(rf(iface,:) - CellCenter(i,j,:)) ! расстояние от центра яч до грани
                        dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:)) ! расстояние до соседней ячейки
                        PE = RLinearInterp(d,dn,P(i,j),P(inc,jnc)) ! Интерполяция P для получ знач в (.) E
                        ! координаты точки E
                        rE(1) = RLinearInterp(d,dn,CellCenter(i,j,1),CellCenter(inc,jnc,1))
                        rE(2) = RLinearInterp(d,dn,CellCenter(i,j,2),CellCenter(inc,jnc,2))
                        ! компоненты градиента в точке E
                        GPE(1) = RLinearInterp(d,dn,GradP(i,j,1),GradP(inc,jnc,1))
                        GPE(2) = RLinearInterp(d,dn,GradP(i,j,2),GradP(inc,jnc,2))
                        ! вычисление P на грани c поправки
                        Pf(iface) = PE + dot_product(GPE(:),rf(iface,:) - rE(:))
                        ! Обновление GradP с учетом итерации
                        GradP_tmp(i,j,:) = GradP_tmp(i,j,:) + Pf(iface)*Sf(iface,:)
                    end do
                    GradP(i,j,:) = GradP_tmp(i,j,:)/CellVolume(i,j)

                end do
            end do
            !$OMP END DO
        end do
    end select
End Subroutine