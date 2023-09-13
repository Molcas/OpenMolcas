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

logical function loopcntr_iterate_cvb()

implicit real*8(a-h,o-z)
#include "seth_cvb.fh"
#include "loopcntr_cvb.fh"
#include "initopt_cvb.fh"
logical unmatched
external istkprobe_cvb
logical istkprobe_cvb

if (iopt2step(ioptim+1) == iopt2step(ioptim)) then
  ioptim = ioptim+1
  goto 1100
end if

joptstep = 0
do ll=1,loopstepmx
  if (joptstep == ioptstep) goto 11
  if ((icode(ll) == 1) .or. (icode(ll) == 3)) joptstep = joptstep+1
end do
11 ll1 = ll
10 continue
! First determine if end of multi-step optimization may have been reached:
if (istkprobe_cvb(istackrep)) then
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
    goto 10
  end if
  ! Restore loop information:
  call istkpush_cvb(istackrep,ioptstep1)
  call istkpush_cvb(istackrep,ioptstep2)
  call istkpush_cvb(istackrep,kk2)
  call istkpush_cvb(istackrep,mxalter)
  call istkpush_cvb(istackrep,italter)
  call istkpush_cvb(istackrep,nconvinone)
  call istkpush_cvb(istackrep,nc_zeroed)
end if
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
        write(6,'(a)') ' Unmatched END or closing bracket!'
        call abend_cvb()
      end if
    end if
    if ((icode(ll) == 1) .or. (icode(ll) == 3)) then
      ! OPTIM / REPORT:
      ioptstep = ioptstep+1
      goto 1000
    else if (icode(ll) == 5) then
      ! ALTERN:
      ! Scan rest of file for number of steps in this loop:
      iend = 0
      ioptstep2 = joptstep
      do kk=ll+1,loopstepmx
        if ((icode(kk) == 1) .or. (icode(kk) == 3)) then
          ioptstep2 = ioptstep2+1
        else if (icode(kk) == 5) then
          iend = iend-1
        else if (icode(kk) == 6) then
          iend = iend+1
        end if
        if (iend == 1) goto 300
      end do
      write(6,*) ' Run-away ENDALTERN or closing bracket!'
      call abend_cvb()
300   continue

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
        goto 1000
      else
        if (nstep > 0) then
          write(6,'(/,a,i4,a)') ' Exiting',nstep,'-step optimization.'
          write(6,'(a,i4)') ' Maximum number of loop iterations reached :',mxalter
        end if
      end if
    end if
  end if
  if ((icode(ll) == 1) .or. (icode(ll) == 3)) joptstep = joptstep+1
end do
ioptstep = noptstep+1

1000 continue
do i=1,noptim
  if (iopt2step(i) == ioptstep) goto 1099
end do
1099 ioptim = i
1100 loopcntr_iterate_cvb = (ioptim <= noptim)

return

end function loopcntr_iterate_cvb
