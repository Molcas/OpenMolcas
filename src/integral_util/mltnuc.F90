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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine MltNuc(CoOP,Chrg,Coor,nAtm,rNucMm,ir)
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments for the nuclei.             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtm, ir
real(kind=wp), intent(in) :: CoOp(3), Chrg(nAtm), Coor(3,nAtm)
real(kind=wp), intent(out) :: rNucMm(nTri_Elem1(ir))
integer(kind=iwp) :: iAtom, ip, ix, iy, iz
real(kind=wp) :: CCoMx, CCoMy, CCoMz, Temp

#ifdef _DEBUGPRINT_
call RecPrt(' In MltNuc:Coor',' ',Coor,3,nAtm)
call RecPrt(' In MltNuc:Chrg',' ',Chrg,nAtm,1)
call RecPrt(' In MltNuc:CoOp',' ',CoOp,1,3)
#endif

! Compute the nuclear contribution to the multipole moments

ip = 0
do ix=ir,0,-1
  do iy=ir-ix,0,-1
    ip = ip+1
    iz = ir-ix-iy
    temp = Zero
    !write(u6,*) ' ix,iy,iz=',ix,iy,iz
    do iAtom=1,nAtm
      if (ix == 0) then
        CCoMx = One
      else
        CCoMx = (Coor(1,iAtom)-CoOp(1))**ix
      end if
      if (iy == 0) then
        CCoMy = One
      else
        CCoMy = (Coor(2,iAtom)-CoOp(2))**iy
      end if
      if (iz == 0) then
        CCoMz = One
      else
        CCoMz = (Coor(3,iAtom)-CoOp(3))**iz
      end if
      !write(u6,*) CCoMx,CCoMy,CCoMz,temp
      temp = temp+Chrg(iAtom)*CCoMx*CCoMy*CCoMz
    end do
    rNucMm(ip) = temp
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' Nuclear Multipole Moments',' ',rNucMm,ip,1)
#endif

end subroutine MltNuc
