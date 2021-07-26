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

subroutine V_EF_PCM(nAt,nTs,DoPot,DoFld,AtmC,Tessera,V,EF_n,EF_e)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs
logical(kind=iwp), intent(in) :: DoPot, DoFld
real(kind=wp), intent(in) :: AtmC(3,nAt), Tessera(4,nTs)
real(kind=wp), intent(out) :: V(nTs), EF_n(3,nTs), EF_e(3,nTs)
integer(kind=iwp) :: nOrdOp

! Compute potential on tesserae

if (DoPot) then
  V(:) = Zero
  nOrdOp = 0
  call Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
end if

! Compute electric field on tesserae

if (DoFld) then
  EF_n(:,:) = Zero
  EF_e(:,:) = Zero
  nOrdOp = 1
  call Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
end if

return

end subroutine V_EF_PCM
