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
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nAObas, nMObas, ipAvred, iMME(nTri3_Elem(MxMltp))
#include "WrkSpc.fh"
integer(kind=iwp) :: iMlt, iMmeMO, iSq, iTEMP, nUniqueM

! First all multipoles are transformed to MO-basis...

call GetMem('Squared','Allo','Real',iSq,nAObas**2)
call GetMem('TEMP','Allo','Real',iTEMP,nAObas*nMObas)
call GetMem('Final','Allo','Real',iMmeMO,nMObas**2)
nUniqueM = 1+3+6
do iMlt=1,nUniqueM
  call Square(Work(iMME(iMlt)),Work(iSq),1,nAObas,nAObas)
  call Dgemm_('T','N',nMObas,nAObas,nAObas,One,Work(ipAvRed),nAObas,Work(iSq),nAObas,Zero,Work(iTEMP),nMObas)
  call Dgemm_('N','N',nMObas,nMObas,nAObas,One,Work(iTEMP),nMObas,Work(ipAvRed),nAObas,Zero,Work(iMmeMO),nMObas)
  call SqToTri_Q(Work(iMmeMO),Work(iMME(iMlt)),nMObas)
end do
call GetMem('Squared','Free','Real',iSq,nAObas**2)
call GetMem('TEMP','Free','Real',iTEMP,nAObas*nMObas)
call GetMem('Final','Free','Real',iMmeMO,nMObas**2)

return

end subroutine MMEtoRMO
