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

use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(inout) :: ifg(4), jfgrd(3,4), jfhss(4,3,4,3), tr(4)
integer(kind=iwp), intent(in) :: nirrep, indgrd(3,4,0:nirrep-1), indhss(4,3,4,3,0:nirrep-1)
integer(kind=iwp), intent(inout) :: jndgrd(3,4,0:nirrep-1), jndhss(4,3,4,3,0:nirrep-1)
real(kind=wp), intent(in) :: Coorm(3,4)
integer(kind=iwp) :: iCent, iiC, iiCar, iiCent, iMax, iMin, iStop, jCent, jjC, jjCent, kCar, kCent, kkCar, lCar, lCent
logical(kind=iwp) :: alike
logical(kind=iwp), external :: EQ

Alike = .false.
if (IfG(1) .and. IfG(2) .and. IfG(3) .and. IfG(4)) then
  do iCent=1,3
    if (.not. Alike) then
      do jCent=iCent+1,4
        if (EQ(CoorM(1,iCent),CoorM(1,jCent))) then
          do kCent=1,4
            iMax = max(jCent,kCent)
            iMin = min(jCent,kCent)
            JndHss(iMax,:,iMin,:,0:nIrrep-1) = 0
            JfHss(iMax,:,iMin,:) = .false.
          end do
          JndGrd(:,jCent,0:nIrrep-1) = 0
          JfGrd(:,jCent) = .false.
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
                  JndHss(iMax,kCar,iMin,1:iStop,0:nIrrep-1) = -IndHss(iMax,kCar,iMin,1:iStop,0:nIrrep-1)
                  JfHss(iMax,kCar,iMin,1:iStop) = .false.
                  do lCar=1,iStop

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
                              JfHss(iiCent,kkCar,jjCent,1:iStop) = .true.
                            end do
                          end if
                        end do ! icent
                        JfGrd(:,iiCent) = .true.
                      end if
                    end do
                  end do
                end do
              end if
            end do
            JndGrd(:,iCent,0:nIrrep-1) = -IndGrd(:,iCent,0:nIrrep-1)
            JfGrd(:,iCent) = .false.
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
  JfGrd(:,:) = .true.
  !JfHss(:,:,:,:) = .true.
  do iiC=1,4
    do jjC=1,iiC
      do iiCar=1,3
        if (iic == jjc) then
          iStop = iiCar
        else
          iStop = 3
        end if
        JfHss(iiC,iiCar,jjc,1:iStop) = .true.
      end do
    end do
  end do

  do kCar=1,3
    do lCent=1,4
      if (lCent == 4) then
        iStop = kCar
      else
        iStop = 3
      end if
      JndHss(4,kCar,lCent,1:iStop,nIrrep-1) = -IndHss(4,kCar,lCent,1:iStop,nIrrep-1)
      JfHss(4,kCar,lCent,1:iStop) = .false.
    end do
  end do

  Jndgrd(:,4,nIrrep-1) = -IndGrd(:,4,nIrrep-1)
  JfGrd(:,4) = .false.

end if

return

end subroutine Translation
