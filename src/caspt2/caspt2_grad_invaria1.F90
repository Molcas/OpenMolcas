!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine caspt2_grad_invaria1(NDPT2,DPT2)
! Put zero to wrong (incomplete) density matrix elements
! If the MRPT2 energy is non-invariant with respect to rotations among the inactive and secondary orbital spaces,
! we cannot determine the off-diagonal elements (in the diagonal block) as in TRDNS2D.
! Therefore, remove the off-diagonal elements. They are computed in caspt2_grad_invaria2 using orbital derivatives.
! The diagonal elements are correct.

use caspt2_module, only: nAsh, nIsh, nOrb, nSsh, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDPT2
real(kind=wp), intent(inout) :: DPT2(NDPT2)
integer(kind=iwp) :: ia, ib, idab, idij, idtu, ii, ij, IOFDAB(8), IOFDIJ(8), is, isym, na, ni, no, ns

IDIJ = 0
do IS=1,NSYM
  NI = NISH(IS)
  NA = NASH(IS)
  NO = NORB(IS)
  IDTU = IDIJ+NO*NI+NI
  IDAB = IDTU+NO*NA+NA
  IOFDIJ(IS) = IDIJ
  IOFDAB(IS) = IDAB
  IDIJ = IDIJ+NO*NO
end do

do isym=1,nsym
  NI = NISH(ISYM)
  NS = NSSH(ISYM)
  NO = NORB(ISYM)
  do ii=1,ni
    do ij=1,ni
      if (ii == ij) cycle
      IDIJ = IOFDIJ(ISYM)+II+NO*(IJ-1)
      DPT2(IDIJ) = Zero
    end do
  end do
  do ia=1,ns
    do ib=1,ns
      if (ia == ib) cycle
      IDAB = IOFDAB(ISYM)+IA+NO*(IB-1)
      DPT2(IDAB) = Zero
    end do
  end do
end do

end subroutine caspt2_grad_invaria1
