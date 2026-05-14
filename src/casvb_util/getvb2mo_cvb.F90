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

subroutine getvb2mo_cvb(orbs)

use casvb_global, only: ifvb, recn_vbwfn
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: orbs(*)
integer(kind=iwp) :: ierr, ioff, ioffs_cvb, ioffs_orbs, ioffs_orbsao, ioffs_orbslao, iorb, kbasiscvb1, nbas_mo1, norb1, nvb1

if (ifvb == 1) call cvbinit_cvb()
call rdheader_cvb(recn_vbwfn,norb1,nbas_mo1,nvb1,kbasiscvb1,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
ioff = 1
do iorb=1,norb1
  call rdgspr_cvb(recn_vbwfn,orbs(ioff),iorb,norb1,1,ierr)
  if (ierr /= 0) then
    write(u6,*) ' Error in VB orbital read :',ierr
    call abend()
  end if
  ioff = ioff+norb1
end do

return

end subroutine getvb2mo_cvb
