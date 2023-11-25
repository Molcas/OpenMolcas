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

function loopcntr_iterate_cvb()

use casvb_global, only: icode, iopt2step, ioptim, ioptstep, ipos, istackrep, joptstep, loopstepmx, noptim, noptstep
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp) :: loopcntr_iterate_cvb
integer(kind=iwp) :: i, iend, ioptstep1, ioptstep2, italter, kk, kk2, ll, ll1, mxalter, nc_zeroed, nconvinone, nstep
logical(kind=iwp) :: done, done2, unmatched
logical(kind=iwp), external :: istkprobe_cvb

if (iopt2step(ioptim+1) == iopt2step(ioptim)) then
  ioptim = ioptim+1
else

  joptstep = 0
  do ll=1,loopstepmx
    if (joptstep == ioptstep) exit
    if ((icode(ll) == 1) .or. (icode(ll) == 3)) joptstep = joptstep+1
  end do
  ll1 = ll
  do
    ! First determine if end of multi-step optimization may have been reached:
    if (.not. istkprobe_cvb(istackrep)) exit
    call istkpop_cvb(istackrep,nc_zeroed)
    call istkpop_cvb(istackrep,nconvinone)
    call istkpop_cvb(istackrep,italter)
    call istkpop_cvb(istackrep,mxalter)
    call istkpop_cvb(istackrep,kk2)
    call istkpop_cvb(istackrep,ioptstep2)
    call istkpop_cvb(istackrep,ioptstep1)
    ! Number of steps is IOPTSTEP2-IOPTSTEP1+1
    if (nconvinone == ioptstep2-ioptstep1+1) then
      ioptstep = ioptstep2
      joptstep = ioptstep2
      ll1 = kk2+1
    else
      ! Restore loop information:
      call istkpush_cvb(istackrep,ioptstep1)
      call istkpush_cvb(istackrep,ioptstep2)
      call istkpush_cvb(istackrep,kk2)
      call istkpush_cvb(istackrep,mxalter)
      call istkpush_cvb(istackrep,italter)
      call istkpush_cvb(istackrep,nconvinone)
      call istkpush_cvb(istackrep,nc_zeroed)
      exit
    end if
  end do
  done2 = .false.
  do ll=ll1,loopstepmx
    if (joptstep == ioptstep) then
      ! Looking for next card after previous OPTIM/REPORT:
      if ((icode(ll) == 2) .or. (icode(ll) == 4)) then
        if (ll == 1) then
          unmatched = .true.
        else
          unmatched = (icode(ll)-icode(ll-1) /= 1)
        end if
        if (unmatched) then
          write(u6,'(a)') ' Unmatched END or closing bracket!'
          call abend_cvb()
        end if
      end if
      if ((icode(ll) == 1) .or. (icode(ll) == 3)) then
        ! OPTIM / REPORT:
        ioptstep = ioptstep+1
        done2 = .true.
        exit
      else if (icode(ll) == 5) then
        ! ALTERN:
        ! Scan rest of file for number of steps in this loop:
        iend = 0
        ioptstep2 = joptstep
        done = .false.
        do kk=ll+1,loopstepmx
          if ((icode(kk) == 1) .or. (icode(kk) == 3)) then
            ioptstep2 = ioptstep2+1
          else if (icode(kk) == 5) then
            iend = iend-1
          else if (icode(kk) == 6) then
            iend = iend+1
          end if
          if (iend == 1) then
            done = .true.
            exit
          end if
        end do
        if (.not. done) then
          write(u6,*) ' Run-away ENDALTERN or closing bracket!'
          call abend_cvb()
        end if

        italter = 1
        mxalter = ipos(ll)
        nconvinone = -1
        nc_zeroed = 0
        call istkpush_cvb(istackrep,ioptstep+1)
        call istkpush_cvb(istackrep,ioptstep2)
        call istkpush_cvb(istackrep,kk)
        call istkpush_cvb(istackrep,mxalter)
        call istkpush_cvb(istackrep,italter)
        call istkpush_cvb(istackrep,nconvinone)
        call istkpush_cvb(istackrep,nc_zeroed)
      else if (icode(ll) == 6) then
        ! ENDALTERN:
        call istkpop_cvb(istackrep,nc_zeroed)
        call istkpop_cvb(istackrep,nconvinone)
        call istkpop_cvb(istackrep,italter)
        call istkpop_cvb(istackrep,mxalter)
        call istkpop_cvb(istackrep,kk2)
        call istkpop_cvb(istackrep,ioptstep2)
        call istkpop_cvb(istackrep,ioptstep1)
        italter = italter+1
        nstep = ioptstep-ioptstep1+1
        if ((nstep > 0) .and. (italter <= mxalter)) then
          ! Next loop iteration:
          call istkpush_cvb(istackrep,ioptstep1)
          call istkpush_cvb(istackrep,ioptstep2)
          call istkpush_cvb(istackrep,kk2)
          call istkpush_cvb(istackrep,mxalter)
          call istkpush_cvb(istackrep,italter)
          call istkpush_cvb(istackrep,nconvinone)
          call istkpush_cvb(istackrep,nc_zeroed)
          ioptstep = ioptstep1
          done2 = .true.
          exit
        else if (nstep > 0) then
          write(u6,'(/,a,i4,a)') ' Exiting',nstep,'-step optimization.'
          write(u6,'(a,i4)') ' Maximum number of loop iterations reached :',mxalter
        end if
      end if
    end if
    if ((icode(ll) == 1) .or. (icode(ll) == 3)) joptstep = joptstep+1
  end do
  if (.not. done2) ioptstep = noptstep+1

  do i=1,noptim
    if (iopt2step(i) == ioptstep) exit
  end do
  ioptim = i

end if

loopcntr_iterate_cvb = (ioptim <= noptim)

return

end function loopcntr_iterate_cvb
