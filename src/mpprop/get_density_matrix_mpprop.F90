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

subroutine Get_Density_Matrix_mpprop(D,nBas,nSym)

#include "intent.fh"

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym)
real(kind=wp), intent(_OUT_) :: D(*)
integer(kind=iwp) :: nDens

if (nSym == 1) then
  nDens = nBas(1)*(nBas(1)+1)/2
  call Get_dArray_chk('D1ao',D,nDens)
# ifdef _DEBUGPRINT_
  call RecPrt('D',' ',D,1,nDens)
# endif
else
  write(u6,*) 'MpProp cannot handle symmetry'
  call Abend()
end if

return

end subroutine Get_Density_Matrix_mpprop
