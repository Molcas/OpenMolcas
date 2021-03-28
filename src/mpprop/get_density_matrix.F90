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

subroutine Get_Density_Matrix_mpprop(ip_D,nDens,nBas,nSym)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
integer nBas(nSym)

if (nSym == 1) then
  nDens = nBas(1)*(nBas(1)+1)/2
  call GetMem('D1ao','Allo','Real',ip_D,nDens)
  call Get_D1ao(Work(ip_D),nDens)
# ifdef _DEBUGPRINT_
  call RecPrt('D',' ',Work(ip_D),1,nDens)
# endif
else
  write(6,*) 'MpProp cannot handle symmetry'
  call Abend()
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(nBas)

end subroutine Get_Density_Matrix_mpprop
