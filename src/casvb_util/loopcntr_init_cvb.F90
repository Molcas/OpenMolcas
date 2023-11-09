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

subroutine loopcntr_init_cvb(inputmode1,initfalse)

use casvb_global, only: icode, icrit, initial, inputmode, iopt2step, ioptcode, ioptim, ioptstep, istackrep, joptstep, lfxvb, &
                        loopstep, loopstepmx, lzrvb, ndrot, nfxorb, nfxvb, nmcscf, noptim, noptstep, norb, norbrel, nort, &
                        nstackrep, nzrvb, ploc, strucopt
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: inputmode1
logical(kind=iwp), intent(in) :: initfalse
integer(kind=iwp) :: i, lll, noptkw, nrepkw
logical(kind=iwp) :: constrained_opt, guess_available, initial_opts, svbfirst
logical(kind=iwp), external :: ifcasci_cvb, & ! ... Files/Hamiltonian available ...
                               up2date_cvb    ! ... Make: up to date? ...

call istkinit_cvb(istackrep,nstackrep)
inputmode = inputmode1
ioptim = 0
ioptstep = 0
if (inputmode == 2) then
  loopstepmx = loopstep
  noptstep = joptstep

  ! Check for "special case" => initial optimizations
  initial_opts = .true.
  guess_available = .false.
  if (nmcscf >= 2) guess_available = .true.
  if (up2date_cvb('WRITEGS')) guess_available = .true.
  if (.not. up2date_cvb('STRTGS')) guess_available = .true.
  if (.not. up2date_cvb('INPGS')) guess_available = .true.

  if (guess_available) initial_opts = .false.
  if (noptstep > 0) initial_opts = .false.
  ! Are we doing a constrained optimization?:
  constrained_opt = .false.
  if (norbrel > 0) constrained_opt = .true.
  if (nort > 0) constrained_opt = .true.
  if (ndrot > 0) constrained_opt = .true.
  if (nfxorb > 0) constrained_opt = .true.
  if (ploc) constrained_opt = .true.
  if ((nfxvb > 0) .or. (lfxvb == 1)) constrained_opt = .true.
  if ((nzrvb > 0) .or. (lzrvb == 1)) constrained_opt = .true.
  ! If INIT/NOINIT keyword was encountered, override:
  if (initial == 0) initial_opts = .false.
  if (initial == 1) initial_opts = .true.
  ! Finally may be overridden by initfalse:
  if (initfalse) initial_opts = .false.
  if (initial_opts) then
    ! IOPTCODE bit 0 = REPORT
    !          bit 1 = OPTIM
    !          bit 2 = Svb
    !          bit 3 = freeze structure coefficients
    !          bit 4 = strong-orthogonality constraints

    ! Should we do Svb optimization first?:
    svbfirst = ifcasci_cvb()

    if (.not. constrained_opt) then
      ! Two first optimizations are SOPP & PP:
      noptim = 0
      if (svbfirst) then
        if (norb > 2) then
          noptim = noptim+1
          ioptcode(noptim) = ibset(ibset(ibset(0,1),2),4)
        end if
        if (strucopt) then
          noptim = noptim+1
          ioptcode(noptim) = ibset(ibset(ibset(0,1),2),3)
          if (noptim == 2) ioptcode(1) = ibset(ibclr(ioptcode(1),3),4)
        end if
      else
        if (norb > 2) then
          noptim = noptim+1
          ioptcode(noptim) = ibset(ibset(0,1),4)
        end if
        if (strucopt) then
          noptim = noptim+1
          ioptcode(noptim) = ibset(ibset(0,1),3)
          if (noptim == 2) ioptcode(1) = ibset(ibclr(ioptcode(1),3),4)
        end if
      end if
    end if
    ! Then a third Svb optimization if we are doing Evb:
    if ((icrit /= 1) .and. svbfirst) then
      noptim = noptim+1
      ioptcode(noptim) = ibset(ibset(0,1),2)
    end if
    ! Finally actual optimization:
    noptim = noptim+1
    ioptcode(noptim) = 2
    ioptcode(noptim) = ibset(0,1)
    ! Add "report":
    noptim = noptim+1
    ioptcode(noptim) = ibset(0,0)

    iopt2step(0) = 0
    iopt2step(1:noptim) = 1
    iopt2step(noptim+1) = noptstep+1
  else
    noptim = noptstep
    ioptcode(1:noptim) = 0
    do i=0,noptim
      iopt2step(i) = i
    end do
    ! Append OPTIM keyword if none present
    noptkw = 0
    do lll=1,loopstepmx
      if (icode(lll) == 1) noptkw = noptkw+1
    end do
    if (noptkw == 0) then
      noptim = noptim+1
      ioptcode(noptim) = ibset(0,1)
      iopt2step(noptim) = iopt2step(noptim-1)
    end if
    ! Append REPORT keyword if none present
    nrepkw = 0
    do lll=1,loopstepmx
      if (icode(lll) == 3) nrepkw = nrepkw+1
    end do
    if (nrepkw == 0) then
      noptim = noptim+1
      ioptcode(noptim) = ibset(0,0)
      iopt2step(noptim) = iopt2step(noptim-1)
    end if
    iopt2step(noptim+1) = noptstep+1
  end if
end if

return

end subroutine loopcntr_init_cvb
