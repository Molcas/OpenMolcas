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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

! Subroutine relating to generalized 1-e density matrix (GD) called
! in CMSNewton
!     GD^KL_tu = sum_{MN}{U^KM * U^LN * GD^MN_tu}
subroutine RotGD(GD,R,nGD,lRoots,NAC2)

use CMS, only: RGD
use Constants, only: Zero, One

implicit none
integer nGD, lRoots, NAC2
real*8 GD(nGD), R(lRoots**2)
integer iNAC2, iLoc, lRoots2
!real*8 RGD(lRoots**2), RGDR(lRoots**2)

lRoots2 = lRoots**2

!write(u6,*) 'rotation matrix in RotGD'
!call RecPrt(' ',' ',R,lRoots,lRoots)

!write(u6,*) 'GD matrix after rotation'
!call RecPrt(' ',' ',GD,lRoots2,NAC2)

do iNAC2=1,NAC2
  iLoc = (iNAC2-1)*lRoots2+1
  call DGEMM_('T','N',lRoots,lRoots,lRoots,One,R,lRoots,GD(iLoc),lRoots,Zero,RGD,lRoots)
  call DGEMM_('N','N',lRoots,lRoots,lRoots,One,RGD,lRoots,R,lRoots,Zero,GD(iLoc),lRoots)
end do

!write(u6,*) 'GD matrix after rotation'
!call RecPrt(' ',' ',GD,lRoots2,NAC2)

return

end subroutine RotGD
