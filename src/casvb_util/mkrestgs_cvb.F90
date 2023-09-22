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

subroutine mkrestgs_cvb(orbsao,irdorbs,cvb,cvbdet,iapr,ixapr,iabind,cvbdet1)

use casvb_global, only: nbas_mo
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: orbsao(nbas_mo,norb), cvb(nvb), cvbdet(ndetvb), cvbdet1(*)
integer(kind=iwp) :: irdorbs(norb), iapr(ndetvb), ixapr(nda+1), iabind(*)
#include "files_cvb.fh"
integer(kind=iwp) :: ia, ib, idetvb1, idum(1), ioffs, iorb, ixa, nalf1, nbet1, ndetvb1, norb1

ioffs = 0
call rdis_cvb(idum,1,recn_tmp04,ioffs)
ndetvb1 = idum(1)
call rdis_cvb(idum,1,recn_tmp04,ioffs)
norb1 = idum(1)
call rdis_cvb(idum,1,recn_tmp04,ioffs)
nalf1 = idum(1)
call rdis_cvb(idum,1,recn_tmp04,ioffs)
nbet1 = idum(1)
if ((norb1 /= norb) .or. (nalf1 /= nalf) .or. (nbet1 /= nbet)) then
  write(u6,'(a)') ' Inconsistency between previous and current VB wavefunction definitions.'
  write(u6,*) ' NORB now ',norb,' before ',norb1
  write(u6,*) ' NALF now ',nalf,' before ',nalf1
  write(u6,*) ' NBET now ',nbet,' before ',nbet1
  call abend_cvb()
end if
do iorb=1,norb
  irdorbs(iorb) = 1
  call rdrs_cvb(orbsao(1,iorb),norb,recn_tmp04,ioffs)
end do
call rdis_cvb(iabind,ndetvb1,recn_tmp04,ioffs)
call rdrs_cvb(cvbdet1,ndetvb1,recn_tmp04,ioffs)

call fzero(cvbdet,ndetvb)
do idetvb1=1,ndetvb1
  ! NDA & string definitions assumed the same:
  ib = (iabind(idetvb1)-1)/nda+1
  ia = iabind(idetvb1)-(ib-1)*nda
  do ixa=ixapr(ia),ixapr(ia+1)-1
    if (ib == iapr(ixa)) cvbdet(ixa) = cvbdet1(idetvb1)
  end do
end do
kbasiscvb = kbasis
call vb2strc_cvb(cvbdet,cvb)

return

end subroutine mkrestgs_cvb
