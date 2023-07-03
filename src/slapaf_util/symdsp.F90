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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

logical function SymDsp(iBsFnc)
!***********************************************************************
!                                                                      *
! Object: to establish if a translation or a rotation belongs to the   *
!         total symmetric irreducible representation.                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             April '92                                                *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iOper

implicit real*8(A-H,O-Z)
integer jPrmt(0:7)
#include "real.fh"
data jPrmt/1,-1,-1,1,-1,1,1,-1/
! Statement function
iPrmt(i,j) = jPrmt(iand(i,j))

!write(6,*) ' iBsFnc=',iBsFnc
SymDsp = .true.
mask = 0
do i=0,nIrrep-1
  do j=1,3
    if (iand(iOper(i),2**(j-1)) /= 0) mask = ior(mask,2**(j-1))
  end do
end do
jBsFnc = iand(mask,iBsFnc)
!write(6,*) ' jBsFnc=',jBsFnc

! Loop over operators

iAcc = 0
do i=0,nIrrep-1
  iAcc = iAcc+iPrmt(iOper(i),jBsFnc)
end do
if (iAcc == 0) SymDsp = .false.

return

end function SymDsp
