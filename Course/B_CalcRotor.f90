Subroutine B_CalcRotor(NI,NJ, CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
    V,RotV)
   real,dimension(0:NI,0:NJ,2):: V
   INTEGER :: NeighCell(4,2)
   REAL SF(4,2), rF(4,2)
   REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   !центры и объемы ячеек
         IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !центры и нормали граней для I-граней
         JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),&   ! центры и нормали граней для J-граней
         RotV(0:NI,0:NJ)
   REAL rE(2), GPE(2), Vf(2)
   REAL  d, dn
  
  
      DO I=1, NI-1
        DO J=1, NJ-1
  
          SF(1,:)=-IFaceVector(I,J,:) ! Подготавливаем вектор грани с учетом внешней нормали
          SF(2,:)=IFaceVector(I+1,J,:)
          SF(3,:)=-JFaceVector(I,J,:)
          SF(4,:)=JFaceVector(I,J+1,:)
  
          rF(1,:)=IFaceCenter(I,J,:) ! координаты  центров граней
          rF(2,:)=IFaceCenter(I+1,J,:)
          rF(3,:)=JFaceCenter(I,J,:)
          rF(4,:)=JFaceCenter(I,J+1,:)
  
          NeighCell(1,:)=[i-1,j]!подготавливаем массив соседних ячеек
          NeighCell(2,:)=[I+1,j]
          NeighCell(3,:)=[i,j-1]
          NeighCell(4,:)=[i,j+1]
  
      RotV(I,J)=0 !инициализация
  
      DO IFACE=1,4 !цикл по всем граням
      IN=NeighCell(iface,1)
      JN=NeighCell(iface,2)
  
      d=Norm2(rF(iface,:)-CellCenter(I,j,:))!расстояние от центра ячейки до центра грани
      dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:))  !расстояние до центра соседней ячейки
  
      Vf(1)=rLinearInterp(d,dn,V(I,J,1),V(IN,JN,1)) !интерполяция скорости на грани
      Vf(2)=rLinearInterp(d,dn,V(I,J,2),V(IN,JN,2))
  
                RotV(i,j) = RotV(i,j) + (Sf(iface,1)*Vf(2) - Sf(iface,2)*Vf(1)) !обновляем ротор
  
              end do
  
              RotV(i,j) = RotV(i,j)/CellVolume(i,j)
  
            end do
          end do
  
  End Subroutine