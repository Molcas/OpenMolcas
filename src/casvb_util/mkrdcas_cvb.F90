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

use casvb_global, only: civb1, civb2, corenrg, formE, icrit, ifhamil, ifinish, ipr, memplenty, projcas, strtci, variat
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iscf
real(kind=wp) :: cnrm, cscf, eexp
logical(kind=iwp) :: need_cas
logical(kind=iwp), external :: ifcasci_cvb, & ! ... Files available ...
                               up2date_cvb, & ! ... Make: up to date? ...
                               valid_cvb

if (ifinish == 0) then
  need_cas = (icrit == 1) .or. projcas
else
  need_cas = ifcasci_cvb() .and. (.not. variat)
end if
if (.not. need_cas) return
! Get CASSCF eigenvector
if (.not. ifcasci_cvb()) then
  if ((ipr(1) >= 0) .and. valid_cvb(strtci)) call prtfid_cvb(' Warning: CI vector not found - no ',strtci)
  if (icrit == 1) then
    write(u6,*) ' No optimization without CASSCF vector!'
    call abend_cvb()
  end if
else
  if (ipr(3) >= 2) write(u6,'(/,a)') ' Read CASSCF eigenvector:'
  call getci_cvb(civb1)
end if
call cinorm2_cvb(civb1,cnrm)
cnrm = One/cnrm
call ciscale2_cvb(civb1,cnrm,iscf,cscf)
if ((.not. up2date_cvb('CASCHECK')) .or. (ipr(3) >= 2)) then
  call untouch_cvb('CASCHECK')
  ! Some checks
  if (abs(cnrm-One) > 1.0e-3_wp) then
    if (ipr(3) >= 0) write(u6,formE) ' WARNING: Norm of CI vector read differs from one :',cnrm
  else if (ipr(3) >= 2) then
    write(u6,formE) ' Norm of CI vector read ',cnrm
  end if
  if ((ipr(3) >= 2) .and. (iscf /= 0)) then
    write(u6,'(a,i6)') ' SCF determinant:',iscf
    write(u6,formE) '     coefficient:',cscf
  end if
  if (ifhamil) then
    call cicopy_cvb(civb1,civb2)
    call applyh_cvb(civb2)
    call cidot_cvb(civb1,civb2,eexp)
    if (ipr(3) >= 1) write(u6,formE) ' CASSCF energy :',eexp+corenrg
    if (ipr(3) >= 1) write(u6,'(a)') ' '
  end if
end if
if (.not. memplenty) call ciwr_cvb(civb1,61001.2_wp)

return

end subroutine mkrdcas_cvb
