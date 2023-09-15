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

subroutine mkrdcas_cvb()

implicit real*8(a-h,o-z)
logical need_cas
! ... Make: up to date? ...
logical, external :: up2date_cvb
! ... Files/Hamiltonian available ...
logical, external :: valid_cvb, ifcasci_cvb, ifhamil_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "frag_cvb.fh"
#include "formats_cvb.fh"
#include "WrkSpc.fh"

if (ifinish == 0) then
  need_cas = (icrit == 1) .or. projcas
else
  need_cas = ifcasci_cvb() .and. (.not. variat)
end if
if (.not. need_cas) return
! Get CASSCF eigenvector
if (.not. ifcasci_cvb()) then
  if ((ip(1) >= 0) .and. valid_cvb(strtci)) call prtfid_cvb(' Warning: CI vector not found - no ',strtci)
  if (icrit == 1) then
    write(6,*) ' No optimization without CASSCF vector!'
    call abend_cvb()
  end if
else
  if (ip(3) >= 2) write(6,'(/,a)') ' Read CASSCF eigenvector:'
  call getci_cvb(work(lc(1)))
end if
call cinorm2_cvb(work(lc(1)),cnrm)
cnrm = one/cnrm
call ciscale2_cvb(work(lc(1)),cnrm,iscf,cscf)
if ((.not. up2date_cvb('CASCHECK')) .or. (ip(3) >= 2)) then
  call untouch_cvb('CASCHECK')
  ! Some checks
  if (abs(cnrm-one) > 1.d-3) then
    if (ip(3) >= 0) write(6,formE) ' WARNING: Norm of CI vector read differs from one :',cnrm
  else if (ip(3) >= 2) then
    write(6,formE) ' Norm of CI vector read ',cnrm
  end if
  if ((ip(3) >= 2) .and. (iscf /= 0)) then
    write(6,'(a,i6)') ' SCF determinant:',iscf
    write(6,formE) '     coefficient:',cscf
  end if
  if (ifhamil_cvb()) then
    call cicopy_cvb(work(lc(1)),work(lc(2)))
    call applyh_cvb(work(lc(2)))
    call cidot_cvb(work(lc(1)),work(lc(2)),eexp)
    if (ip(3) >= 1) write(6,formE) ' CASSCF energy :',eexp+corenrg
    if (ip(3) >= 1) write(6,'(a)') ' '
  end if
end if
if (.not. memplenty) call ciwr_cvb(work(lc(1)),61001.2d0)

return

end subroutine mkrdcas_cvb
