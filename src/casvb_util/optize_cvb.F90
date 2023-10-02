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

subroutine optize_cvb(fx,ioptc,iter,imethod,isadinp,mxiter,maxinp,corenrg,ipinp,ipdd1,ipdd2,strucopt)
!***********************************************************************
!*                                                                     *
!*  Routine for second-order optimization.                             *
!*                                                                     *
!*  Uses the following functions:                                      *
!*                                                                     *
!*   GETFREE:   Determine # free parameters                            *
!*    UPDATE:   Add dx to vector of free variables, x                  *
!*        FX:   Calculate f(x)                                         *
!*      GRAD:   Calculate gradient of f                                *
!*      HESS:   Calculate the action of the hessian on trial vector    *
!*   GETHESS:   Calculate full hessian                                 *
!*                                                                     *
!*  IOPTC is optimization control :                                    *
!*                                                                     *
!*  IOPTC=-3    Opt. terminated close to convergence (at request)      *
!*  IOPTC=-2    Optimization failed -- too small step size             *
!*  IOPTC=-1    Maximum number of iterations used                      *
!*  IOPTC= 0    Converged                                              *
!*  IOPTC= 1    Not complete                                           *
!*                                                                     *
!***********************************************************************

use casvb_global, only: expct, fxbest, hh, hhkeep, hhstart, ip, isaddle, ix, maxize
use casvb_interfaces, only: opta_sub, optb_sub
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: fx, corenrg
integer(kind=iwp) :: ioptc, iter, imethod, isadinp, mxiter, ipinp, ipdd1, ipdd2
logical(kind=iwp) :: maxinp, strucopt
#include "WrkSpc.fh"
integer(kind=iwp) :: ifollow, maxd, mxit, n_div, nfrdim, nfrdim_dav, nparm, nparm_dav
logical(kind=iwp) :: done, iter_is_1
integer(kind=iwp), external :: mstackr_cvb
procedure(opta_sub) :: dum_a_cvb, o10a_cvb, o123a_cvb, o12ea_cvb, o12sa_cvb, o7a_cvb
procedure(optb_sub) :: o10b_cvb, o123b_cvb, o12eb_cvb, o12sb_cvb, o5b_cvb, o7b_cvb, o8b_cvb

if (mxiter == 0) then
  ioptc = -1
  return
end if

! Initialize for new optimization - input parameters:
isaddle = isadinp
maxize = maxinp
ip = ipinp

! Parameters initialized:
expct = Zero
hh = hhstart
hhkeep = hh
ioptc = 1

if (maxize) then
  ifollow = 1
else
  ifollow = 2
end if

call fx_cvb(fx,.false.)
fxbest = fx
done = .false.
do iter=1,mxiter
  iter_is_1 = (iter == 1)
  call getfree_cvb(nparm,n_div,nfrdim,iter,fx)
  if (nfrdim <= 0) then
    if (ip >= 0) write(u6,'(a)') ' No free parameters!'
    if (ip >= 0) write(u6,'(a)') ' Calculation converged.'
    ioptc = 0
    return
  end if
  if ((fx < Zero) .and. maxize) then
    call chgsgn_cvb(fx)
    call getfree_cvb(nparm,n_div,nfrdim,iter,fx)
  end if

  if ((imethod == 1) .or. (imethod == 2) .or. (imethod == 3)) then
    ix(1) = mstackr_cvb(nparm)
    ix(2) = mstackr_cvb(nparm)
    ix(3) = mstackr_cvb(nparm*nparm)
    ix(4) = mstackr_cvb(nparm)
    ix(5) = mstackr_cvb(nparm)
    ix(6) = mstackr_cvb(nparm)
    ix(7) = mstackr_cvb(nparm)
    call optize2_cvb(fx,nparm,ioptc,work(ix(1)),work(ix(2)),iter_is_1,o123a_cvb,o123b_cvb)
    call mfreer_cvb(ix(1))
  else if (imethod == 5) then
    ix(1) = mstackr_cvb(nparm)
    ix(2) = mstackr_cvb(nparm)
    call optize2_cvb(fx,nparm,ioptc,work(ix(1)),work(ix(2)),iter_is_1,dum_a_cvb,o5b_cvb)
    call mfreer_cvb(ix(1))
  else if (imethod == 7) then
    ix(1) = mstackr_cvb(nparm+1)
    ix(2) = mstackr_cvb(nparm+1)
    maxd = min(nparm+1,200)
    mxit = 500
    call ddinit_cvb('AxEx',nparm+1,nfrdim+1,maxd,mxit,ifollow,isaddle,ipdd1,Zero,n_div)
    call asonC7init_cvb(ix(2),ipdd2)
    call optize2_cvb(fx,nparm,ioptc,work(ix(1)),work(ix(2)),iter_is_1,o7a_cvb,o7b_cvb)
    call mfreer_cvb(ix(1))
  else if (imethod == 8) then
    ix(1) = mstackr_cvb(nparm)
    ix(2) = mstackr_cvb(nparm)
    ix(3) = mstackr_cvb((nparm+1)*(nparm+1))
    ix(4) = mstackr_cvb(nparm+1)
    call optize2_cvb(fx,nparm,ioptc,work(ix(1)),work(ix(2)),iter_is_1,dum_a_cvb,o8b_cvb)
    call mfreer_cvb(ix(1))
  else if (imethod == 9) then
    call optize9_cvb(fx,nparm,ioptc)
  else if (imethod == 10) then
    ix(1) = mstackr_cvb(nparm)
    ix(2) = mstackr_cvb(nparm)
    maxd = min(nparm,200)
    mxit = 500
    call ddinit_cvb('AxExb',nparm,nfrdim,maxd,mxit,ifollow,isaddle,ipdd1,Zero,n_div)
    call asonc10init_cvb(ipdd2)
    call optize2_cvb(fx,nparm,ioptc,work(ix(1)),work(ix(2)),iter_is_1,o10a_cvb,o10b_cvb)
    call mfreer_cvb(ix(1))
  else if ((imethod == 12) .and. maxize) then
    if (strucopt) then
      nparm_dav = nparm
      nfrdim_dav = nfrdim
    else
      nparm_dav = nparm+1
      nfrdim_dav = nfrdim+1
    end if
    ix(1) = mstackr_cvb(nparm_dav)
    ix(2) = mstackr_cvb(nparm_dav)
    maxd = min(nparm_dav,200)
    mxit = 500
    call ddinit_cvb('Axb',nparm_dav,nfrdim_dav,maxd,mxit,ifollow,isaddle,ipdd1,Zero,0)
    call asonc12sinit_cvb(ipdd2)
    call optize2_cvb(fx,nparm_dav,ioptc,work(ix(1)),work(ix(2)),iter_is_1,o12sa_cvb,o12sb_cvb)
    call mfreer_cvb(ix(1))
  else if ((imethod == 12) .and. (.not. maxize)) then
    if (strucopt) then
      nparm_dav = nparm
      nfrdim_dav = nfrdim
    else
      nparm_dav = nparm+1
      nfrdim_dav = nfrdim+1
    end if
    ix(1) = mstackr_cvb(nparm_dav)
    ix(2) = mstackr_cvb(nparm_dav)
    maxd = min(nparm_dav,200)
    mxit = 500
    call ddinit_cvb('AxESx',nparm_dav,nfrdim_dav,maxd,mxit,ifollow,isaddle,ipdd1,corenrg,n_div)
    call asonc12einit_cvb(ipdd2)
    call optize2_cvb(fx,nparm_dav,ioptc,work(ix(1)),work(ix(2)),iter_is_1,o12ea_cvb,o12eb_cvb)
    call mfreer_cvb(ix(1))
  else
    write(u6,*) ' Unrecognized optimization algorithm!',imethod
    call abend_cvb()
  end if
  if (ioptc <= 0) then
    done = .true.
    exit
  end if
end do
if (.not. done) iter = iter-1
if (ioptc >= 0) call getfree_cvb(nparm,n_div,nfrdim,-1,fx)
if ((iter == mxiter) .and. (ioptc > 0)) ioptc = -1

return

end subroutine optize_cvb
