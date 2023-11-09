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

subroutine trust_cvb(iopth,opth,maxize,fx,fxbest,expc,hh,dxnrm,ioptc,scalesmall1,close2conv,converged,skipupd)

use casvb_global, only: cpropt, delopth1, delopth2, dfxmin, formAD, hhaccfac, hhkeep, hhmax, hhopt, hhrejfac, hhtol, ipr, nopth1, &
                        nopth2, scalesmall, zzacclim, zzrejmax, zzrejmin
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: iopth, ioptc
logical(kind=iwp), intent(inout) :: opth
logical(kind=iwp), intent(in) :: maxize, close2conv, converged
real(kind=wp), intent(in) :: fx, fxbest, expc, dxnrm
real(kind=wp), intent(out) :: hh
logical(kind=iwp), intent(out) :: scalesmall1, skipupd
integer(kind=iwp) :: icprbst, icprbst2, idum, iop, ioptst, ipu, nopth
real(kind=wp) :: cprbst, dfx, dum, gap2, hh_min, hhlargest, oldstep, scl, zz
logical(kind=iwp) :: dfx_ok, zz_ok

call zz_cvb(dum,zz,fx,fxbest,expc,-1)
skipupd = .false.
if (.not. close2conv) then
  ipu = 1
else
  ipu = 2
end if
nopth = nopth1(ipu)+nopth2(ipu)
scalesmall1 = scalesmall(ipu)

! Trust region control, based on the quadratic model
iop = mod(iopth,nopth)
if (iop == 0) iop = iop+nopth
if (iopth > 0) cpropt(iop) = fx
! Set HHOPT to HH actually used (might have been scaled):
if (iopth > 0) hhopt(iop) = dxnrm
do
  if ((iopth > 0) .and. (mod(iopth,nopth) == 0) .and. opth) then
    ! Optimisation of trust region completed
    opth = .false.
    if (maxize) then
      call findmx_cvb(cpropt,nopth,cprbst,icprbst)
    else
      call findmn_cvb(cpropt,nopth,cprbst,icprbst)
    end if
    dfx = cprbst-fxbest
    dfx_ok = ((dfx > dfxmin(ipu)) .and. maxize) .or. ((dfx < -dfxmin(ipu)) .and. (.not. maxize))
    zz_ok = (zz > zzrejmin(ipu)) .and. (zz < zzrejmax(ipu))
    if (dfx_ok .and. zz_ok) then
      ! << Accepting update >>
      iopth = 0
      ! Restore HH as used before (so repeated update vector will be
      ! exactly the same):
      if (icprbst <= nopth1(ipu)) then
        hh = hhkeep*(One+(real(icprbst,kind=wp)-Half*real(nopth1(ipu)+1,kind=wp))*delopth1(ipu))
      else if (icprbst <= nopth) then
        ! IFG: nopth1 was used in these two calls, probably a bug
        if (maxize) then
          call findmx_cvb(cpropt,nopth,cprbst,icprbst2)
        else
          call findmn_cvb(cpropt,nopth,cprbst,icprbst2)
        end if
        hh = hhkeep*(One+(real(icprbst2,kind=wp)-Half*real(nopth1(ipu)+1,kind=wp))*delopth1(ipu))
        gap2 = hhkeep*delopth1(ipu)*delopth2(ipu)
        hh = hh+gap2*(real(icprbst-nopth1(ipu),kind=wp)-Half*real(nopth2(ipu)+1,kind=wp))
      end if
      hh = min(hh,hhmax(ipu))
      ! Scale trust region size according to ZZ:
      if (zz < zzacclim(1,ipu)) then
        scl = hhaccfac(1,ipu)
      else if (zz > zzacclim(4,ipu)) then
        scl = hhaccfac(5,ipu)
      else if (zz < zzacclim(2,ipu)) then
        scl = hhaccfac(2,ipu)
      else if (zz > zzacclim(3,ipu)) then
        scl = hhaccfac(4,ipu)
      else
        scl = hhaccfac(3,ipu)
      end if
      if (nopth > 1) then
        if ((hhopt(icprbst) > 1.0e-8_wp) .and. (hh/hhopt(icprbst) > Two)) then
          hhkeep = hh
        else
          ! Scale trust region size according to ZZ:
          hhkeep = hh*scl
        end if
      else
        oldstep = hhopt(icprbst)
        ! Scale trust region size according to ZZ:
        if (scl >= One) then
          hhkeep = max(hhkeep,oldstep*scl)
        else
          hhkeep = hhkeep*scl
        end if
      end if
      skipupd = (icprbst == nopth)
      exit
    else
      ! << Rejecting update >>
      if (converged) then
        iopth = 0
        hh = Zero
        return
      end if
      if (ipr(3) >= 1) write(u6,'(a)') ' Rejecting step.'
      call findmn_cvb(hhopt,nopth,hh_min,idum)
      hhkeep = min(hh_min,hhkeep)*hhrejfac(ipu)
      gap2 = hhkeep*delopth1(ipu)*delopth2(ipu)
      if (nopth2(ipu) == 0) gap2 = Zero
      hhlargest = hhkeep*(One+(real(nopth1(ipu),kind=wp)-Half*real(nopth1(ipu)+1,kind=wp))*delopth1(ipu))+ &
                  gap2*(real(nopth-nopth1(ipu),kind=wp)-Half*real(nopth2(ipu)+1,kind=wp))
      if (hhlargest < hhtol(ipu)) then
        if (ipr(3) >= 0) then
          write(u6,formAD) ' Trust region size smaller than tolerance !',hhlargest,hhtol(ipu)
          write(u6,'(a)') ' Calculation NOT converged!'
        end if
        ioptc = -2
        return
      end if
    end if
  else
    if ((iopth == 0) .and. (nopth > 1) .and. (ipr(3) >= 2)) write(u6,'(/,a)') ' Optimising trust region size :'
    opth = .true.
    iopth = iopth+1
    ioptst = mod(iopth,nopth)
    if (ioptst == 0) ioptst = ioptst+nopth

    if (ioptst <= nopth1(ipu)) then
      hh = hhkeep*(One+(real(ioptst,kind=wp)-Half*real(nopth1(ipu)+1,kind=wp))*delopth1(ipu))
    else if (ioptst <= nopth) then
      ! IFG: nopth1 was used in these two calls, probably a bug
      if (maxize) then
        call findmx_cvb(cpropt,nopth,cprbst,icprbst)
      else
        call findmn_cvb(cpropt,nopth,cprbst,icprbst)
      end if
      hh = hhkeep*(One+(real(icprbst,kind=wp)-Half*real(nopth1(ipu)+1,kind=wp))*delopth1(ipu))
      gap2 = hhkeep*delopth1(ipu)*delopth2(ipu)
      hh = hh+gap2*(real(ioptst-nopth1(ipu),kind=wp)-Half*real(nopth2(ipu)+1,kind=wp))
    end if
    hh = min(hh,hhmax(ipu))
    hhopt(ioptst) = hh
    exit
  end if
end do

return

end subroutine trust_cvb
