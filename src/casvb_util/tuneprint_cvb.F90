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

subroutine tuneprint_cvb()

use casvb_global, only: alftol, cnrmtol, delopth1, delopth2, dfx, dfxmin, dfxtol, dx, eigwrngtol, endwhenclose, exp12tol, follow, &
                        grd, grdwrngtol, hhaccfac, hhmax, hhrejfac, hhstart, hhtol, imethod, ipr, mxdav, nopth1, nopth2, nortiter, &
                        orththr, resthr, safety, scalesmall, sgn, signtol, singul, zzacclim, zzmax, zzmin, zzrejmax, zzrejmin
use Definitions, only: u6

implicit none

if (ipr(3) < 3) return
if (imethod /= 4) then
  write(u6,'(/,a)') ' -------- Details of parameters used by 2nd-order optimizer: -------------'
  write(u6,'(/,a,/)') ' General parameters:'
  call fout_cvb(safety,'SAFETY','Alpha safety in denominator, (H - alpha * I):')
  call fout_cvb(cnrmtol,'CNRMTOL','Tolerance for size of update:')
  call fout_cvb(signtol,'SIGNTOL','Tolerance for sign of Hessian eigenvalues:')
  call fout_cvb(alftol,'ALFTOL','Convergence criterion on alpha:')
  call fout_cvb(dfxtol,'DFXTOL','DFX tolerance for act/exp ratio:')
  call fout_cvb(exp12tol,'EXP12TOL','Criterion on expected change of f(x):')
  call fout_cvb(grdwrngtol,'GRDWRNGTOL','Gradient tol. for scaling small updates:')
  call fout_cvb(eigwrngtol,'EIGWRNGTOL','Eigenvalue tol. for scaling small updates:')
  call lout_cvb(endwhenclose,'ENDIFCLOSE','Exit if optimization close to convergence?')
  write(u6,'(/,a)') ' Convergence criteria:'
  write(u6,'(/,a)') ' Elements of arrays:'
  write(u6,'(a)') ' (1) ... Optimization is in global region.'
  write(u6,'(a)') ' (2) ... Optimization is in local region.'
  write(u6,'(a)') ' (3) ... Optimization is close to wrong stationary point.'
  call fouti_cvb(singul,3,'SINGUL(3)','Thresholds for sing. Hessian (max abs eig):')
  write(u6,'(/,a)') ' Elements of arrays:'
  write(u6,'(a)') ' (*,1) ... Global region, non-singular Hessian.'
  write(u6,'(a)') ' (*,2) ... Global region, singular Hessian.'
  write(u6,'(a)') ' (*,3) ... Local region, non-singular Hessian.'
  write(u6,'(a)') ' (*,4) ... Local region, singular Hessian.'
  write(u6,'(a)') ' (*,5) ... Wrong stationary point, non-singular Hessian.'
  write(u6,'(a)') ' (*,6) ... Wrong stationary point, singular Hessian.'
  call fouti_cvb(sgn,6,'SIGN(6)','Threshold for sign of Hessian eigenvalues:')
  call fouti_cvb(zzmin,6,'ZZMIN(6)','Mininum allowed act/exp ratio:')
  call fouti_cvb(zzmax,6,'ZZMAX(6)','Maximum allowed act/exp ratio:')
  call fouti_cvb(dfx,6,'DFX(6)','Maximum allowed change in f(x):')
  write(u6,'(/,a)') ' Elements of arrays:'
  write(u6,'(a)') ' (1,*) ... Use maximum absolute value in vector.'
  write(u6,'(a)') ' (2,*) ... Use norm of vector.'
  write(u6,'(a)') ' (3,*) ... Use RMS of elements in vector.'
  call foutij_cvb(dx,3,6,'DX(3,6)','Maximum allowed change in variables:')
  call foutij_cvb(grd,3,6,'GRD(3,6)','Maximum allowed gradient:')
  write(u6,'(/,a,/)') ' Trust region control:'
  call fout_cvb(hhstart,'HHSTART','Initial trust region size:')
  call iout_cvb(nopth1(1),'NOPTH1(1)','Number of steps (primary trust size opt):')
  call iout_cvb(nopth2(1),'NOPTH2(1)','Number of steps (secondary trust size opt):')
  call iout_cvb(nopth1(2),'NOPTH1(2)','Number of steps (primary trust size opt):')
  call iout_cvb(nopth2(2),'NOPTH2(2)','Number of steps (secondary trust size opt):')
  call fouti_cvb(delopth1,2,'DELOPTH1(2)','Primary change of trust region size:')
  call fouti_cvb(delopth2,2,'DELOPTH2(2)','Secondary change of trust region size:')
  call fouti_cvb(hhmax,2,'HHMAX(2)','Maximum allowed trust region size:')
  call fouti_cvb(zzrejmin,2,'ZZREJMIN(2)','Minimum allowed act/exp ratio:')
  call fouti_cvb(zzrejmax,2,'ZZREJMAX(2)','Maximum allowed act/exp ratio:')
  call fouti_cvb(dfxmin,2,'DFXMIN(2)','Minimum allowed change in f(x):')
  call fouti_cvb(hhrejfac,2,'HHREJFAC(2)','Trust region size scale factor for rejections:')
  call foutij_cvb(zzacclim,4,2,'ZZACCLIM(4,2)','Act/exp regions for scaling accepted steps:')
  call foutij_cvb(hhaccfac,5,2,'HHACCFAC(5,2)','Trust scale factors for accepted steps:')
  call fouti_cvb(hhtol,2,'HHTOL(2)','Minimum allowed trust region size:')
  call lout_cvb(scalesmall(1),'SCALESMALL(1)','Scale predicted steps smaller than trust?')
  call lout_cvb(scalesmall(2),'SCALESMALL(2)','Scale predicted steps smaller than trust?')
  write(u6,'(/,a)') ' -------------------------------------------------------------------------'
else if (imethod == 4) then
  write(u6,'(/,a,/)') ' -------- Details of parameters used by Davidson optimizer: --------------'
  call iout_cvb(mxdav,'MXDAV','Maxium dimension of Davidson subspace:')
  call fout_cvb(resthr,'RESTHR','Convergence criterion on residual norm:')
  call lout_cvb(follow,'FOLLOW','Root following (for excited states):')
  call fout_cvb(orththr,'ORTHTHR','Tolerance for orthogonality between vectors:')
  call iout_cvb(nortiter,'NORTITER','Maximum number of orthogonalization attempts:')
  write(u6,'(/,a)') ' -------------------------------------------------------------------------'
end if

return

end subroutine tuneprint_cvb
