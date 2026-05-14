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

subroutine mo_transp(cmo,cmo_t,no,nv,ndel,nbas)
! CMO(p,alpha) <- CMO_t(alpha,p+del),  p=o+v

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no, nv, ndel, nbas
real(kind=wp), intent(out) :: cmo(no+nv,nbas)
real(kind=wp), intent(in) :: cmo_t(nbas,no+nv+ndel)
integer(kind=iwp) :: i

do i=1,nbas
  cmo(:,i) = cmo_t(i,1:no+nv)
end do

return

end subroutine mo_transp
