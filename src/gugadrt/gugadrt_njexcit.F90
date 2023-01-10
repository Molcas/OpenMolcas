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

subroutine gugadrt_njexcit(indjk,ljk,iextbit,nextbit,ivalid,jstep,kttmp,k0)

use gugadrt_global, only: iref_occ, n_ref
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ljk, iextbit, nextbit, jstep, k0
integer(kind=iwp), intent(inout) :: indjk(ljk), kttmp
integer(kind=iwp), intent(out) :: ivalid
integer(kind=iwp) :: i, idxref, inm, itexcit(n_ref), ival, kp

kp = k0
do idxref=1,n_ref
  call upacknod(indjk,idxref,ival,nextbit,iextbit,ljk)
  if (jstep == 1 .or. jstep == 2) then
    if (iref_occ(kp+1,idxref) == 0) ival = ival+1
  end if
  if (jstep == 3) then
    if (iref_occ(kp+1,idxref) == 0) ival = ival+2
    if (iref_occ(kp+1,idxref) == 1) ival = ival+1
  end if
  if (ival > 2) ival = 3
  itexcit(idxref) = ival
end do
inm = minval(itexcit)

if (inm > 2) then
  ivalid = 0
else
  kttmp = inm
  ivalid = 1
  if (jstep /= 0) then
    do i=1,n_ref
      ival = itexcit(i)
      call packnod(indjk,i,ival,nextbit,iextbit,ljk)
    end do
  end if
end if

return

end subroutine gugadrt_njexcit
