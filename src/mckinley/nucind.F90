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

subroutine NucInd(coor,kdc,ifgrd,ifhss,indgrd,indhss,jfgrd,jfhss,jndgrd,jndhss,tr,ifg)

use McKinley_global, only: sIrrep
use Index_Functions, only: iTri
use Symmetry_Info, only: nIrrep
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Coor(3,4)
integer(kind=iwp) :: kdc, IndGrd(0:2,0:1,0:(nIrrep-1)), IndHss(0:1,0:2,0:1,0:2,0:(nIrrep-1)), JndGrd(0:2,0:3,0:(nIrrep-1)), &
                     JndHss(0:3,0:2,0:3,0:2,0:(nIrrep-1))
logical(kind=iwp) :: IfGrd(0:2,0:1), IfHss(0:1,0:2,0:1,0:2), JfGrd(0:2,0:3), JfHss(0:3,0:2,0:3,0:2), Tr(0:3), IfG(0:3)
#include "Molcas.fh"
#include "disp.fh"
integer(kind=iwp) :: iAtom, iCar, iCent, iComp, iIrrep, iStop, jAtom, jCar, kCar, kCent, Maxi, Mini, nDisp, nnIrrep
logical(kind=iwp), external :: EQ, TF

!                                                                      *
!***********************************************************************
!                                                                      *
call ICopy(nIrrep*16*9,[0],0,JndHss,1)
call iCopy(nIrrep*4*3,[0],0,JndGrd,1)
call LCopy(144,[.false.],0,jfHss,1)
call LCopy(4,[.false.],0,Tr,1)
call LCopy(12,[.false.],0,jfGrd,1)

! COPY CNTLR MATRIXES

do iAtom=0,1
  do iCar=0,2
    JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
    do iIrrep=0,nIrrep-1
      JndGrd(iCar,iAtom,iIrrep) = IndGrd(iCar,iAtom,iIrrep)
    end do
    do jAtom=0,1
      do jCar=0,2
        JfHss(iAtom,iCar,jAtom,jCar) = IfHss(iAtom,iCar,jAtom,jCar)
        do iIrrep=0,nIrrep-1
          JndHss(iAtom,iCar,jAtom,jCar,iIrrep) = IndHss(iAtom,iCar,jAtom,jCar,iIrrep)
        end do ! iirrep
      end do !jCar
    end do ! jAtom
  end do !iCar
end do !iAtom

! Derivatives with respect to the operator is computed via the translational invariance.

nnIrrep = nIrrep
if (sIrrep) nnIrrep = 1
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(kdc,iIrrep)
  do iCar=0,2
    iComp = 2**iCar
    if (TF(kdc,iIrrep,iComp)) then
      nDisp = nDisp+1

      ! Reset flags for the basis set centers so that we
      ! will explicitly compute the derivatives with
      ! respect to those centers. Activate flag for the
      ! third center so that its derivative will be computed
      ! by the translational invariance.

      JndGrd(iCar,0,iIrrep) = abs(JndGrd(iCar,0,iIrrep))
      JndGrd(iCar,1,iIrrep) = abs(JndGrd(iCar,1,iIrrep))
      JndGrd(iCar,2,iIrrep) = -nDisp
      JfGrd(iCar,0) = .true.
      JfGrd(iCar,1) = .true.
      JfGrd(iCar,2) = .false.
    else
      JndGrd(iCar,2,iIrrep) = 0
    end if
  end do
end do

! The third center is calculated by translation invariance
! This requires the 2nd derivatives on the other centers.


do iCar=0,2
  do jAtom=0,2
    if (jAtom == 2) then
      iStop = iCar
    else
      iStop = 2
    end if
    do jCar=0,iStop
      do iIrrep=0,nIrrep-1
        if ((JndGrd(iCar,2,iIrrep) /= 0) .and. (JndGrd(jCar,jAtom,iIrrep) /= 0)) then
          JndHss(2,iCar,jAtom,jCar,iIrrep) = -iTri(abs(JndGrd(iCar,2,iIrrep)),abs(JndGrd(jCar,jAtom,iIrrep)))

          Tr(2) = .true.
          if (jAtom == 2) then
            Maxi = max(iCar,jCar)
            Mini = min(iCar,jCar)
            jfHss(0,Maxi,0,Mini) = .true.
            jfHss(1,Maxi,1,Mini) = .true.
            jfHss(1,iCar,0,jCar) = .true.
            jfHss(1,jCar,0,iCar) = .true.
          else
            Maxi = max(iCar,jCar)
            Mini = min(iCar,jCar)
            jfHss(jAtom,Maxi,jAtom,Mini) = .true.
            jfHss(1,iCar,0,jCar) = .true.
            jfHss(1,jCar,0,iCar) = .true.
          end if ! jAtom == 2
        end if ! if indgrd
      end do  ! iirrep
    end do ! jCar
  end do ! jAtom
end do ! iCar

IfG(0) = .true.
IfG(1) = .true.
IfG(2) = .false.
IfG(3) = .false.
do iCent=0,1
  if (EQ(Coor(1,iCent+1),Coor(1,3))) then
    IfG(iCent) = .false.
    do iCar=0,2
      jfGrd(iCar,iCent) = .false.
      do kCar=0,2
        do kCent=0,3
          jfHss(iCent,iCar,kCent,kCar) = .false.
          jfHss(kCent,kCar,iCent,iCar) = .false.
          do iIrrep=0,nIrrep-1
            jndHss(iCent,iCar,kCent,kCar,iIrrep) = 0
            jndHss(kCent,kCar,iCent,iCar,iIrrep) = 0
          end do !iIrrep
        end do ! kcent
      end do !kCar
      do iIrrep=0,nIrrep-1
        jndGrd(iCar,iCent,iIrrep) = 0
      end do !iIrrep
    end do !ICat
  end if ! uf eq
end do !icent

return

end subroutine NucInd
