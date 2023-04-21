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

subroutine ExtT1(H,T1,dima,adda)
! this routine does:
! Extract H(a',i) <- T1(a,i) for given aGrp
!
! parameter description:
! H    - Output file (O)
! T1   - T1 amplitudes (I)
! dima - dimension of given Group (I)
! adda - shift of a' in full a set (I)
!
! N.B. Kvajt odflaknute

use chcc_global, only: no, nv
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, adda
real(kind=wp), intent(out) :: H(dima,no)
real(kind=wp), intent(in) :: T1(nv,no)

H(:,:) = T1(adda+1:adda+dima,:)

return

end subroutine ExtT1
