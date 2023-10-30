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

subroutine setjobiph_cvb(nel_c,norb_c,i2s_c,isym_c,neltot_c)

use casvb_global, only: iorclos_c, iorcore_c, iorocc_c, istms2_c, istnel_c, istsy_c, mcore_c, mxstt_ci, nstats_c, nstsym_c, weight_c
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nel_c, norb_c, i2s_c, isym_c, neltot_c
#include "rasdim.fh"
#include "jobiph_j.fh"
integer(kind=iwp) :: i, j
real(kind=wp) :: wgt

! Orbitals  --  OCC, CLOSED and CORE cards
iorcore_c(:) = nfro_j(:)
iorclos_c(:) = nish_j(:)
iorocc_c(:) = nrs2_j(:)
! States  --  WF, STATE and WEIGHT cards
nstsym_c = 1
weight_c(:,:) = Zero
do i=1,lroots_j
  wgt = Zero
  do j=1,nroots_j
    if (iroot_j(j) == i) wgt = weight_j(j)
  end do
  if (wgt /= Zero) then
    if (i > mxstt_ci) then
      write(u6,*) ' Root number too large in casrecov_cvb :',i,mxstt_ci
      call abend_cvb()
    end if
  end if
  weight_c(i,1) = wgt
end do

istnel_c(1) = nactel_j
istsy_c(1) = lsym_j
istms2_c(1) = ispin_j-1
nstats_c(1) = lroots_j
! Set derived info
nel_c = nactel_j
i2s_c = ispin_j-1
isym_c = lsym_j
norb_c = sum(nrs2_j(:))
mcore_c = sum(nfro_j(:)+nish_j(:))
neltot_c = nel_c+2*mcore_c
! MO common block:
call setmocom_cvb()

return

end subroutine setjobiph_cvb
