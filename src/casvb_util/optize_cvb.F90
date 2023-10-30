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

subroutine optize_cvb(fx,ioptc,iter,imethod,isadinp,mxiter,maxinp,ipinp,ipdd1,ipdd2,strucopt)
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

use casvb_global, only: corenrg, eigval, eigvec, expct, fxbest, hh, hhkeep, hhstart, ip, ipp10, ipp12e, ipp12s, ipp7, isaddleo, &
                        iter10, iter12e, iter12s, iter7, maxize, odx, odxp, ograd, ogradp, owrk
use casvb_interfaces, only: opta_sub, optb_sub
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: fx
integer(kind=iwp), intent(out) :: ioptc, iter
integer(kind=iwp), intent(in) :: imethod, isadinp, mxiter, ipinp, ipdd1, ipdd2
logical(kind=iwp), intent(in) :: maxinp, strucopt
integer(kind=iwp) :: ifollow, maxd, mxit, n_div, nfrdim, nfrdim_dav, nparm, nparm_dav
logical(kind=iwp) :: done, iter_is_1
procedure(opta_sub) :: dum_a_cvb, o10a_cvb, o123a_cvb, o12ea_cvb, o12sa_cvb, o7a_cvb
procedure(optb_sub) :: o10b_cvb, o123b_cvb, o12eb_cvb, o12sb_cvb, o5b_cvb, o7b_cvb, o8b_cvb

if (mxiter == 0) then
  ioptc = -1
  return
end if

! Initialize for new optimization - input parameters:
isaddleo = isadinp
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
    call mma_allocate(odx,nparm,label='odx')
    call mma_allocate(ograd,nparm,label='ograd')
    call mma_allocate(eigvec,nparm,nparm,label='eigvec')
    call mma_allocate(eigval,nparm,label='eigval')
    call mma_allocate(odxp,nparm,label='odxp')
    call mma_allocate(ogradp,nparm,label='ogradp')
    call mma_allocate(owrk,nparm,label='owrk')
    call optize2_cvb(fx,nparm,ioptc,iter_is_1,o123a_cvb,o123b_cvb)
    call mma_deallocate(odxp)
    call mma_deallocate(ogradp)
    call mma_deallocate(owrk)
  else if (imethod == 5) then
    call mma_allocate(odx,nparm,label='odx')
    call mma_allocate(ograd,nparm,label='ograd')
    call optize2_cvb(fx,nparm,ioptc,iter_is_1,dum_a_cvb,o5b_cvb)
  else if (imethod == 7) then
    call mma_allocate(odx,nparm+1,label='odx')
    call mma_allocate(ograd,nparm+1,label='ograd')
    maxd = min(nparm+1,200)
    mxit = 500
    call ddinit_cvb('AxEx',nparm+1,nfrdim+1,maxd,mxit,ifollow,isaddleo,ipdd1,Zero,n_div)
    iter7 = 0
    ipp7 = ipdd2
    call orthcvb_init_cvb()
    call optize2_cvb(fx,nparm,ioptc,iter_is_1,o7a_cvb,o7b_cvb)
    call ddclean_cvb()
  else if (imethod == 8) then
    call mma_allocate(odx,nparm,label='odx')
    call mma_allocate(ograd,nparm,label='ograd')
    call mma_allocate(eigvec,nparm+1,nparm+1,label='eigvec')
    call mma_allocate(eigval,nparm+1,label='eigval')
    call optize2_cvb(fx,nparm,ioptc,iter_is_1,dum_a_cvb,o8b_cvb)
  else if (imethod == 9) then
    call optize9_cvb(fx,nparm,ioptc)
  else if (imethod == 10) then
    call mma_allocate(odx,nparm,label='odx')
    call mma_allocate(ograd,nparm,label='ograd')
    maxd = min(nparm,200)
    mxit = 500
    call ddinit_cvb('AxExb',nparm,nfrdim,maxd,mxit,ifollow,isaddleo,ipdd1,Zero,n_div)
    iter10 = 0
    ipp10 = ipdd2
    call orthcvb_init_cvb()
    call optize2_cvb(fx,nparm,ioptc,iter_is_1,o10a_cvb,o10b_cvb)
    call ddclean_cvb()
  else if ((imethod == 12) .and. maxize) then
    if (strucopt) then
      nparm_dav = nparm
      nfrdim_dav = nfrdim
    else
      nparm_dav = nparm+1
      nfrdim_dav = nfrdim+1
    end if
    call mma_allocate(odx,nparm_dav,label='odx')
    call mma_allocate(ograd,nparm_dav,label='ograd')
    maxd = min(nparm_dav,200)
    mxit = 500
    call ddinit_cvb('Axb',nparm_dav,nfrdim_dav,maxd,mxit,ifollow,isaddleo,ipdd1,Zero,0)
    iter12s = 0
    ipp12s = ipdd2
    call optize2_cvb(fx,nparm_dav,ioptc,iter_is_1,o12sa_cvb,o12sb_cvb)
    call ddclean_cvb()
  else if ((imethod == 12) .and. (.not. maxize)) then
    if (strucopt) then
      nparm_dav = nparm
      nfrdim_dav = nfrdim
    else
      nparm_dav = nparm+1
      nfrdim_dav = nfrdim+1
    end if
    call mma_allocate(odx,nparm_dav,label='odx')
    call mma_allocate(ograd,nparm_dav,label='ograd')
    maxd = min(nparm_dav,200)
    mxit = 500
    call ddinit_cvb('AxESx',nparm_dav,nfrdim_dav,maxd,mxit,ifollow,isaddleo,ipdd1,corenrg,n_div)
    iter12e = 0
    ipp12e = ipdd2
    call optize2_cvb(fx,nparm_dav,ioptc,iter_is_1,o12ea_cvb,o12eb_cvb)
    call ddclean_cvb()
  else
    write(u6,*) ' Unrecognized optimization algorithm!',imethod
    call abend_cvb()
  end if
  if (allocated(odx)) call mma_deallocate(odx)
  if (allocated(ograd)) call mma_deallocate(ograd)
  if (allocated(eigval)) call mma_deallocate(eigval)
  if (allocated(eigvec)) call mma_deallocate(eigvec)
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
