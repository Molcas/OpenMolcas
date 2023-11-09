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

subroutine ddinit_cvb(method,nparm1,nfrdim1,maxd1,mxit1,ifollow1,isaddle1,ip1,corenrg1,n_div1)

use casvb_global, only: ap, axc, c, corenrg, ifollow, ipdd, isaddledd, maxd, mxit, mxrhs, n_div, nfrdim, nortiterdd, nparm, nroot, &
                        orththrdd, res, resthrdd, rhs, rhsp, solp, solp_res, sxc
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: method
integer(kind=iwp), intent(in) :: nparm1, nfrdim1, maxd1, mxit1, ifollow1, isaddle1, ip1, n_div1
real(kind=wp), intent(in) :: corenrg1

! Input parameters:
nparm = nparm1
nfrdim = nfrdim1
maxd = maxd1
mxit = mxit1
ifollow = ifollow1
isaddledd = isaddle1
nroot = max(1,isaddle1+1)
ipdd = ip1
corenrg = corenrg1
n_div = n_div1
! Defaults:
resthrdd = 1.0e-5_wp
! Local DIRDIAG parameters:
orththrdd = 1.0e-10_wp
nortiterdd = 50

call mma_allocate(c,nparm,maxd,label='c')
call mma_allocate(res,nparm,label='res')
call mma_allocate(solp,maxd,label='solp')
call mma_allocate(solp_res,maxd,label='solp_res')
if (method == 'AxEx') then
  mxrhs = 0
  ! Arrays:  C       AxC     RES
  !          AP      SOLP    SOLP_RES
  call mma_allocate(axc,nparm,maxd,label='axc')
  call mma_allocate(ap,maxd,maxd,label='ap')
else if (method == 'AxESx') then
  mxrhs = 0
  ! Arrays:  C       AxC     SxC     RES
  !          AP      SOLP    SOLP_RES
  call mma_allocate(axc,nparm,maxd,label='axc')
  call mma_allocate(sxc,nparm,maxd,label='sxc')
  call mma_allocate(ap,maxd,maxd,label='ap')
else if (method == 'Axb') then
  mxrhs = 1
  ! Arrays:  C       SxC     RES     RHS
  !          RHSP    SOLP    SOLP_RES
  call mma_allocate(sxc,nparm,maxd,label='sxc')
  call mma_allocate(rhs,nparm,mxrhs,label='rhs')
  call mma_allocate(rhsp,maxd,label='rhsp')
else if (method == 'AxExb') then
  mxrhs = 1
  ! Arrays:  C       AxC     RES     RHS
  !          AP      RHSP    SOLP    SOLP_RES
  call mma_allocate(axc,nparm,maxd,label='axc')
  call mma_allocate(rhs,nparm,mxrhs,label='rhs')
  call mma_allocate(ap,maxd,maxd,label='ap')
  call mma_allocate(rhsp,maxd,label='rhsp')
end if

return

end subroutine ddinit_cvb
