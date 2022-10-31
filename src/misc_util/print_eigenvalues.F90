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

subroutine Print_EigenValues(H,nH)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nH
real(kind=wp) :: H(nTri_Elem(nH))
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ipEVal, ipEVec

call GetMem('EVal','Allo','Real',ipEVal,nH*(nH+1)/2)
call GetMem('EVec','Allo','Real',ipEVec,nH*nH)

! Copy elements for H

call dcopy_(nH*(nH+1)/2,H,1,Work(ipEVal),1)

! Set up a unit matrix

call dcopy_(nH*nH,[Zero],0,Work(ipEVec),1)
call dcopy_(nH,[One],0,Work(ipEVec),nH+1)

! Compute eigenvalues and eigenvectors

call Jacob(Work(ipEVal),Work(ipEVec),nH,nH)
call Jacord(Work(ipEVal),Work(ipEVec),nH,nH)

! Print out the result

write(u6,*)
write(u6,*) 'Eigenvalues of the matrix'
write(u6,*)
write(u6,'(10F15.8)') (Work(i*(i+1)/2+ipEVal-1),i=1,nH)

call GetMem('EVec','Free','Real',ipEVec,nH*nH)
call GetMem('EVal','Free','Real',ipEVal,nH*(nH+1)/2)

return

end subroutine Print_EigenValues
