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

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ip_D, nDens
integer(kind=iwp), intent(in) :: nSym, nBas(nSym)
#include "WrkSpc.fh"

if (nSym == 1) then
  nDens = nBas(1)*(nBas(1)+1)/2
  call GetMem('D1ao','Allo','Real',ip_D,nDens)
  call Get_D1ao(Work(ip_D),nDens)
# ifdef _DEBUGPRINT_
  call RecPrt('D',' ',Work(ip_D),1,nDens)
# endif
else
  write(u6,*) 'MpProp cannot handle symmetry'
  call Abend()
end if

return

end subroutine Get_Density_Matrix_mpprop
