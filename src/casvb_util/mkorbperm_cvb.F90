!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mkorbperm_cvb()

use casvb_global, only: cvb, cvbdet, iorbprm, ipr, norb, orbs, owrk2
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iorb, jorb
real(kind=wp) :: sgn

if (ipr(1) >= 1) then
  write(u6,'(/,a)') ' Permuting orbitals :'
  write(u6,'(1x,30i4)') iorbprm(1:norb)
end if
do iorb=1,norb
  jorb = abs(iorbprm(iorb))
  sgn = real(sign(1,iorbprm(iorb)),kind=wp)
  owrk2(:,iorb) = orbs(:,jorb)
  owrk2(:,iorb) = sgn*owrk2(:,iorb)
end do
orbs(:,:) = owrk2(:,:)
call str2vbc_cvb(cvb,cvbdet)
call permvb_cvb(cvbdet,iorbprm)
call vb2strc_cvb(cvbdet,cvb)

return

end subroutine mkorbperm_cvb
