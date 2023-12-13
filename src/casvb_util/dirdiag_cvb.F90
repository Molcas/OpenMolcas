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

subroutine dirdiag_cvb(ddasonc,ddsol,ddres,ddres2upd,c,axc,sxc,share,vec,res,rhs,ap,rhsp,solp,solp_res,symm,use_a,use_rhs,maxdav, &
                       n,nprmdim,nvguess,nvrestart,isaddle,ifollow,mxiter,resthr,orththr,nortiter,corenrg,ioptc,iter,fx,ip)
!***********************************************************************
!*                                                                     *
!*  Routine for "direct" diagonalization.                              *
!*                                                                     *
!*  Uses the following functions:                                      *
!*                                                                     *
!*  "External" :                                                       *
!*     DDASonC:    Apply A and S onto vector(s)                        *
!*     DDSOL:      Solve linear equation in Davidson subspace          *
!*     DDRES:      From solution, evaluate residual vector             *
!*     DDRES2UPD:  Convert residual into next Davidson vector          *
!*     DDRESTART:  Restart when Davidson dimension becomes max         *
!*                                                                     *
!*  Other :                                                            *
!*     ABEND:      Error exit                                          *
!*     MXATB:      Matrix multiply                                     *
!*     DDOT:       Level 1 BLAS                                        *
!*     DNRM2:      Level 1 BLAS                                        *
!*                                                                     *
!*  If no metric is used (S=1), C and SxC may share memory. If they    *
!*  do, SHARE should be set true.                                      *
!*                                                                     *
!*  IFOLLOW controls root selecting :                                  *
!*                                                                     *
!*     IFOLLOW=1   Maximization                                        *
!*     IFOLLOW=2   Minimization                                        *
!*     IFOLLOW=3   Overlap-based root following                        *
!*     IFOLLOW=4   Eigenvalue-based root following                     *
!*                                                                     *
!*  For IFOLLOW=1 the root chosen will be number ISADDLE+1 from the    *
!*  end, for IFOLLOW=2, number ISADDLE+1 from the beginning.           *
!*                                                                     *
!*  IOPTC is optimization control :                                    *
!*                                                                     *
!*     IOPTC=-3    Opt. terminated close to convergence (at request)   *
!*     IOPTC=-2    Optimization failed -- too small step size          *
!*     IOPTC=-1    Maximum number of iterations used                   *
!*     IOPTC= 0    Converged                                           *
!*     IOPTC= 1    Not complete                                        *
!*                                                                     *
!*  NVRESTART on entry: the number of Davidson vectors in a "restart"  *
!*  run. These are assumed to be obtained from a previous call to      *
!*  DIRDIAG, and the relevant C, AxC, SxC, AP and RHSP quantities      *
!*  must be unchanged.                                                 *
!*                                                                     *
!*  NVRESTART on exit: the final dimension of the Davidson space --    *
!*  can be used as input in a following restart run.                   *
!*                                                                     *
!*  NVGUESS represents the total number of guess vectors (thus         *
!*  NVGUESS>=NVRESTART).                                               *
!*                                                                     *
!***********************************************************************

use casvb_global, only: formAD, formAF
use casvb_interfaces, only: ddasonc_sub, ddsol_sub, ddres_sub, ddres2upd_sub
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
procedure(ddasonc_sub) :: ddasonc
procedure(ddsol_sub) :: ddsol
procedure(ddres_sub) :: ddres
procedure(ddres2upd_sub) :: ddres2upd
integer(kind=iwp), intent(in) :: maxdav, n, nprmdim, isaddle, ifollow, mxiter, nortiter, ip
real(kind=wp), intent(inout) :: c(n,maxdav), axc(n,maxdav), sxc(n,maxdav)
logical(kind=iwp), intent(in) :: share, symm, use_a, use_rhs
real(kind=wp), intent(out) :: vec(n), res(n), ap(maxdav,maxdav), rhsp(maxdav), solp(maxdav), solp_res(maxdav), fx
real(kind=wp), intent(in) :: rhs(n), resthr, orththr, corenrg
integer(kind=iwp), intent(inout) :: nvguess, nvrestart
integer(kind=iwp), intent(out) :: ioptc, iter
integer(kind=iwp) :: i, ifail, imacro, iorth, itdav, ivres, ndavvec, nroot
real(kind=wp) :: e1, eig, eig_res, facn, resnrm
logical(kind=iwp) :: done, is_converged
real(kind=wp), external :: ddot_, dnrm2_

e1 = Zero

if (ip >= 2) then
  write(u6,'(/,a)') ' Starting Davidson optimization.'
  write(u6,'(a)') ' -------------------------------'
end if
if (ip >= 1) write(u6,'(a,i5)') ' Maximum dimension of Davidson subspace:',maxdav

nroot = max(1,isaddle+1)
ifail = -1
if ((ifollow <= 2) .and. (maxdav < max(nroot,min(2*nroot,nprmdim)))) then
  ifail = max(nroot,min(2*nroot,nprmdim))
else if (maxdav < min(2,nprmdim)) then
  ifail = min(2,nprmdim)
end if
if (ifail /= -1) then
  write(u6,'(a)') ' Davidson dimension too small!'
  write(u6,'(a,i3,a)') ' Need storage for at least',ifail,' vectors.'
  call abend_cvb()
end if

do ivres=1,nvrestart
  if (use_a) then
    do i=1,ivres
      ap(i,ivres) = ddot_(n,c(:,i),1,axc(:,ivres),1)
      if (.not. symm) then
        ap(ivres,i) = ddot_(n,c(:,ivres),1,axc(:,i),1)
      else
        ap(ivres,i) = ap(i,ivres)
      end if
    end do
  end if
  if (use_rhs) rhsp(ivres) = ddot_(n,c(:,ivres),1,rhs,1)
end do

iter = 0
itdav = 0
! MXITER max possible no of macro iterations, normally opt will skip
outer: do imacro=1,mxiter
  do itdav=max(nvrestart,1),min(maxdav,mxiter-iter)
    if (itdav > nvrestart) then
      iter = iter+1
      ! Ensure accurate orthogonalization:
      facn = dnrm2_(n,c(:,itdav),1)
      done = .false.
      do iorth=1,nortiter
        c(:,itdav) = c(:,itdav)/facn
        call schmidtd2_cvb(c,sxc,itdav-1,c(:,itdav),1,n)
        call ddproj_cvb(c(:,itdav),n)
        facn = dnrm2_(n,c(:,itdav),1)
        if (abs(One-facn) < orththr) then
          done = .true.
          exit
        end if
      end do
      if (.not. done) then
        write(u6,*) ' Not able to achieve orthonormality in max number of attempts:',nortiter
        call abend_cvb()
      end if
      c(:,itdav) = c(:,itdav)/facn

      call ddasonc(c(:,itdav),axc(:,itdav),sxc(:,itdav),1,n)

      facn = sqrt(ddot_(n,c(:,itdav),1,sxc(:,itdav),1))
      c(:,itdav) = c(:,itdav)/facn
      if (.not. share) sxc(:,itdav) = sxc(:,itdav)/facn
      if (use_a) then
        axc(:,itdav) = axc(:,itdav)/facn

        do i=1,itdav
          ap(i,itdav) = ddot_(n,c(:,i),1,axc(:,itdav),1)
          if (.not. symm) then
            ap(itdav,i) = ddot_(n,c(:,itdav),1,axc(:,i),1)
          else
            ap(itdav,i) = ap(i,itdav)
          end if
        end do
      end if
    end if
    if (use_rhs) rhsp(itdav) = ddot_(n,c(:,itdav),1,rhs,1)

    call ddsol(ap,rhsp,itdav,maxdav,nprmdim,solp,solp_res,eig,eig_res)

    if (ip >= 2) write(u6,formAF) ' Optimal eigenvalue :',eig+corenrg
    if (e1 == Zero) e1 = eig

    if ((ip >= 3) .or. ((n <= 100) .and. (ip == 2))) then
      write(u6,'(a)') ' Current vector :'
      call vecprint_cvb(c(:,itdav),n)
    end if
    if (.not. (itdav < nvguess)) then
      call ddres(axc,sxc,rhs,res,solp_res,maxdav,n,itdav,eig_res,is_converged)

      call ddproj_cvb(res,n)
      if (is_converged) then
        resnrm = dnrm2_(n,res,1)
        is_converged = (resnrm < resthr)
      end if
      if (is_converged) then
        if (ip >= 0) write(u6,formAD) ' Converged ... residual norm:',resnrm
      else
        if (ip >= 2) then
          write(u6,'(a)') ' '
          write(u6,formAD) ' Residual norm:',resnrm
          write(u6,'(a)') ' '
        end if
      end if
      if ((ip >= 3) .or. ((n <= 100) .and. (ip == 2))) then
        write(u6,'(a)') ' Residual vector :'
        call vecprint_cvb(res,n)
      end if
      if (is_converged) then
        ioptc = 0
        exit outer
      end if
      if (itdav <= maxdav-1) call ddres2upd(res,c(:,itdav+1),n)
    end if
  end do
  if (iter >= mxiter) then
    if (ip >= 0) write(u6,'(2a,i5,a)') ' Davidson optimization not converged in ',mxiter,' iterations'
    ioptc = -1
    exit outer
  end if
  call ddrestart_cvb(c,axc,vec,ap,solp,maxdav,n,nvguess,nvrestart)
end do outer
ndavvec = min(itdav,maxdav)
call mxatb_cvb(c,solp,n,ndavvec,1,vec)

if (ip >= 2) write(u6,formAD) ' Total eigenvalue change in Davidson :',eig-e1
if ((ip >= 3) .and. (n <= 500)) then
  write(u6,'(a)') ' Final projected solution vector :'
  call vecprint_cvb(solp,n)
  write(u6,'(a)') ' Final solution vector :'
  call vecprint_cvb(vec,n)
end if
fx = eig+corenrg
! Prepare for restart:
nvguess = 0
nvrestart = itdav

return

end subroutine dirdiag_cvb
