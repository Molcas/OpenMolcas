!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Translation(ifg,jfgrd,jfhss,tr,jndgrd,jndhss,coorm,nirrep,indgrd,indhss)

implicit integer(a-z)
logical jfhss(4,3,4,3), jfgrd(3,4), tr(4), ifg(4), eq, alike
integer jndgrd(3,4,0:nirrep-1), jndhss(4,3,4,3,0:nirrep-1)
integer indgrd(3,4,0:nirrep-1), indhss(4,3,4,3,0:nirrep-1)
real*8 Coorm(3,4)

Alike = .false.
if (IfG(1) .and. IfG(2) .and. IfG(3) .and. IfG(4)) then
  do iCent=1,3
    if (.not. Alike) then
      do jCent=iCent+1,4
        if (EQ(CoorM(1,iCent),CoorM(1,jCent))) then
          do kCent=1,4
            iMax = max(jCent,kCent)
            iMin = min(jCent,kCent)
            do kCar=1,3
              do lCar=1,3
                do mIrrep=0,nIrrep-1
                  JndHss(iMax,kCar,iMin,lCar,mIrrep) = 0
                end do
                JfHss(iMax,kCar,iMin,lCar) = .false.
              end do
            end do
          end do
          do jCar=1,3
            do mIrrep=0,nIrrep-1
              JndGrd(jCar,jCent,mIrrep) = 0
            end do
            jfGrd(jCar,jCent) = .false.
          end do
          IfG(jCent) = .false.
          if (.not. Alike) then
            IfG(iCent) = .false.
            Tr(iCent) = .true.
            Alike = .true.
            do kCent=1,4
              if ((.not. EQ(CoorM(1,iCent),CoorM(1,kCent))) .or. (kCent == iCent)) then
                iMax = max(iCent,kCent)
                iMin = min(iCent,kCent)
                do kCar=1,3
                  if (iMax == iMin) then
                    iStop = kCar
                  else
                    iStop = 3
                  end if
                  do lCar=1,iStop
                    do mIrrep=0,nIrrep-1
                      JndHss(iMax,kCar,iMin,lCar,mIrrep) = -IndHss(iMax,kCar,iMin,lCar,mIrrep)
                    end do
                    JfHss(iMax,kCar,iMin,lCar) = .false.

                    ! Set the derivatives that are needed for the translation
                    ! invariance calculations.

                    do iiCent=1,4
                      if (.not. EQ(CoorM(1,iiCent),CoorM(1,iCent))) then
                        do jjCent=1,iiCent
                          if (.not. EQ(CoorM(1,jjCent),CoorM(1,iCent))) then
                            do kkCar=1,3
                              if (iiCent == jjCent) then
                                iStop = kkCar
                              else
                                iStop = 3
                              end if
                              do llCar=1,iStop
                                JfHss(iiCent,kkCar,jjCent,llCar) = .true.
                              end do
                            end do
                          end if
                        end do ! icent
                        do kkCar=1,3
                          JfGrd(kkCar,iiCent) = .true.
                        end do
                      end if
                    end do
                  end do
                end do
              end if
            end do
            do jCar=1,3
              do mIrrep=0,nIrrep-1
                JndGrd(jCar,iCent,mIrrep) = -IndGrd(jCar,iCent,mIrrep)
              end do
              JfGrd(jCar,iCent) = .false.
            end do
          end if
        end if
      end do
    end if
  end do
end if

! If all centers are different delete the fourth center

if (.not. Alike) then
  IfG(4) = .false.
  Tr(4) = .true.
  call lCopy(12,[.true.],0,JFGRD,1)
  !call lCopy(144,[.true.],0,JFHss,1)
  do iiC=1,4
    do jjC=1,iiC
      do iiCar=1,3
        iStop = 3
        if (iic == jjc) iStop = iiCar
        do jjCar=1,iStop
          JfHss(iiC,iiCar,jjc,jjCar) = .true.
        end do
      end do
    end do
  end do

  do kCar=1,3
    do lCent=1,4
      iStop = 3
      if (lCent == 4) iStop = kCar
      do lCar=1,iStop
        do mIrrep=0,nirrep-1
          JndHss(4,kCar,lCent,lCar,mIrrep) = -IndHss(4,kCar,lCent,lCar,mIrrep)
        end do
        JfHss(4,kCar,lCent,lCar) = .false.
      end do
    end do
  end do

  do lCar=1,3
    do mIrrep=0,nirrep-1
      Jndgrd(lCar,4,mIrrep) = -IndGrd(lCar,4,mIrrep)
    end do
    jfgrd(lCar,4) = .false.
  end do

end if

return

end subroutine Translation
