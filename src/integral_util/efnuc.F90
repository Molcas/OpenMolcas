!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine EFNuc(CoOP,Chrg,Coor,nAtm,ESIT,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the electricstatic interaction tensor contribution*
!         from the nuclei. In the case that the test charge coincide   *
!         with a nucleau we will remove that center.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, April '95.                                      *
!***********************************************************************

use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
integer nAtm, nOrdOp
real*8 Chrg(nAtm), Coor(3,nAtm), ESIT((nOrdOp+1)*(nOrdOp+2)/2)
integer, allocatable :: C_ESIT(:)
real*8 CoOp(3)
integer nTot, iPowr, iAtom, ix, iy, iz
real*8 Fact, x, y, z, r2, r, eix, eiy, eiz, temp
real*8, parameter :: Thr = 1.0e-12_wp
#ifdef _DEBUGPRINT_
integer n, nElem
! Statement function
nElem(n) = (n+1)*(n+2)/2
#endif

! Compute the nuclear contribution to the electrostatic interation
! tensor, ESIT.

ESIT(:) = Zero

nTot = (nOrdOp+1)**6
call mma_allocate(C_ESIT,nTot,Label='ESIT')
call InitIA(C_ESIT,nOrdOp)

iPowR = 2*nOrdOp+1
Fact = One
if (nOrdOp >= 1) Fact = -One
do iAtom=1,nAtm
  x = CoOp(1)-Coor(1,iAtom)
  y = CoOp(2)-Coor(2,iAtom)
  z = CoOp(3)-Coor(3,iAtom)
  r2 = x**2+y**2+z**2
  if (r2 > Thr) then
    r = Chrg(iAtom)/sqrt(r2)**iPowR
    do ix=nOrdOp,0,-1
      do iy=nOrdOp-ix,0,-1
        iz = nOrdOp-ix-iy
        if (ix == 0) then
          EIx = One
        else
          EIx = x**ix
        end if
        if (iy == 0) then
          EIy = One
        else
          EIy = y**iy
        end if
        if (iz == 0) then
          EIz = One
        else
          EIz = z**iz
        end if
        temp = Fact*EIx*EIy*EIz*r

        call ContEI(C_ESIT,nOrdOp,ESIT,ix,iy,iz,temp)

      end do
    end do       ! End loop over cartesian combinations
  end if
end do           ! End loop over atoms

call mma_deallocate(C_ESIT)

#ifdef _DEBUGPRINT_
call RecPrt(' The Electrostatic Interaction Tensor',' ',ESIT,nElem(nOrdOp),1)
#endif

return

end subroutine EFNuc
