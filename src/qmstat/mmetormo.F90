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

subroutine MMEtoRMO(nAObas,nMObas,ipAvRed,iMME)

use qmstat_global, only: MxMltp
use Index_functions, only: nTri3_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAObas, nMObas, ipAvred, iMME(nTri3_Elem(MxMltp))
#include "WrkSpc.fh"
integer(kind=iwp) :: iMlt, nUniqueM
real(kind=wp), allocatable :: MmeMO(:,:), Sq(:,:), TEMP(:,:)

! First all multipoles are transformed to MO-basis...

call mma_allocate(Sq,nAObas,nAObas,label='Squared')
call mma_allocate(TEMP,nMObas,nAObas,label='TEMP')
call mma_allocate(MmeMO,nMObas,nMObas,label='Final')
nUniqueM = 1+3+6
do iMlt=1,nUniqueM
  call Square(Work(iMME(iMlt)),Sq,1,nAObas,nAObas)
  call Dgemm_('T','N',nMObas,nAObas,nAObas,One,Work(ipAvRed),nAObas,Sq,nAObas,Zero,TEMP,nMObas)
  call Dgemm_('N','N',nMObas,nMObas,nAObas,One,TEMP,nMObas,Work(ipAvRed),nAObas,Zero,MmeMO,nMObas)
  call SqToTri_Q(MmeMO,Work(iMME(iMlt)),nMObas)
end do
call mma_deallocate(Sq)
call mma_deallocate(TEMP)
call mma_deallocate(MmeMO)

return

end subroutine MMEtoRMO
