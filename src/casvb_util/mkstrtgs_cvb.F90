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

subroutine mkstrtgs_cvb(orbsao,irdorbs,cvb,recn,kbasis1)

use casvb_global, only: nbas_mo

implicit real*8(a-h,o-z)
logical use_ao, ifmos_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
dimension orbsao(nbas_mo,norb), irdorbs(norb), cvb(*)

call rdheader_cvb(recn,norb1,nbas_mo1,nvb1,kbasiscvb,ioffs_orbs,ioffs_cvb,ioffs_orbsao,ioffs_orbslao)
use_ao = ifmos_cvb() .and. ((.not. variat) .or. (variat .and. (nmcscf == 1))) .and. (nbas_mo == nbas_mo1) .and. (ioffs_orbsao > 0)
do iorb=1,norb
  if (.not. use_ao) then
    irdorbs(iorb) = 1
    call rdgspr_cvb(recn,orbsao(1,iorb),iorb,norb,1,ierr)
  else
    irdorbs(iorb) = 2
    call rdgspr_cvb(recn,orbsao(1,iorb),iorb,nbas_mo,3,ierr)
  end if
  if (ierr /= 0) then
    call prtfid_cvb(' Error in orbital read from ',recn)
    write(6,'(a)') ' Orbital no :',iorb
    if (use_ao) then
      write(6,'(a)') ' AO basis ? : Yes'
    else
      write(6,'(a)') ' AO basis ? : No'
    end if
    call abend_cvb()
  end if
end do
call rdgspr_cvb(recn,cvb,1,nvb,2,ierr)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(kbasis1)

end subroutine mkstrtgs_cvb
