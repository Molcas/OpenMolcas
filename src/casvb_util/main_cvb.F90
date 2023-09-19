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

subroutine main_cvb()

implicit real*8(a-h,o-z)
! ... Make: up to date? ...
logical, external :: up2date_cvb
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
external loopcntr_iterate_cvb
logical loopcntr_iterate_cvb

if (service) return
ibase = mstackr_cvb(0)
if (variat) nmcscf = nmcscf+1
call stat1_cvb()

! ---------  make objects  ---------
call makefile_cvb()
if (nmcscf <= 1) call touch_cvb('WRITEGS')
! ----------------------------------

call change0_cvb()

call loopcntr_init_cvb(2,.true.)
call input_cvb()
call loopcntr_init_cvb(2,.false.)
do while (loopcntr_iterate_cvb())
  call input_cvb()

  if (variat .and. (.not. endvar) .and. (ip(6) < 2)) then
    ! Reduce output level for main variational iterations:
    do i=1,10
      ip(i) = -1
    end do
  end if

  if (endvar .and. (.not. up2date_cvb('PRTSUM'))) then
    ! End of variational calculation
    if (ip(1) >= 0) write(6,'(/,a)') ' CASVB -- summary of results :'
    if (ip(1) >= 0) write(6,'(a)') ' -----------------------------'
    call make_cvb('PRTSUM')
  end if

  call make_cvb('STAT')

  call touch_cvb('ORBFREE')
  call touch_cvb('CIFREE')

  if (ifinish <= 2) call change_cvb()

  call casinfoprint_cvb()
  call cnfprint_cvb()
  call prtopt_cvb()

  if (ifinish <= 2) call make_cvb('INIT')

  ! -------  make dependencies -------
  if (nort > 0) then
    call depend_cvb('ORBFREE','ORBS')
  else
    call undepend_cvb('ORBFREE','ORBS')
  end if
  call depend_cvb('CIFREE','CVB')
  ! ----------------------------------

  if (ifinish == 0) then
    call opt_cvb()
    call ncset_cvb(ioptc_new)
  else if ((ifinish == 1) .or. (ifinish == 2)) then
    call reprt_cvb()
    call ncset_cvb(0)
  end if
  call writegs_cvb()
end do

call stat2_cvb()

call mfreer_cvb(ibase)

return

end subroutine main_cvb
