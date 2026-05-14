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

subroutine GetChVHlp4(H,T1,dima,adda)
! this routine does:
! Extract H(i,a') <- T1(a,i) for given aGrp
!
! parameter description:
! H       - Output file (O)
! T1      - T1 amplitudes (I)
! dima    - dimension of given Group (I)
! adda    - shift of a' in full a set (I)
!
! N.B. Kvajt odflaknute, je to ExtT1 len transponovane

use chcc_global, only: no, nv
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, adda
real(kind=wp), intent(out) :: H(no,dima)
real(kind=wp), intent(in) :: T1(nv,no)
integer(kind=iwp) :: a

do a=1,dima
  H(:,a) = T1(adda+a,:)
end do

return

end subroutine GetChVHlp4
