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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!> S = H C through CI_Lib_Sigma with the integrals of MCLR: Int1 over all
!> orbitals in square symmetry blocks (block (iS,Mul(iS,iOp)), nOrb(iS) x
!> nOrb(Mul(iS,iOp)), column-major, in the order of iS), of which only the
!> active part enters; Int2 over the active orbitals, already in the layout
!> of CI_Lib_Sigma (the full triangle, or for Square the square pair
!> index). Int1 is converted here to TU of CI_Lib_Sigma: for Square the
!> full square matrix, otherwise the triangle (element (t,u), t >= u).
subroutine MCLR_CI_Sigma(iOpSpin,iCSym,iSSym,Square,Have_2_el,nInt1,Int1,nInt2,Int2,nC,C,nS,S)

use CI_Operations, only: CI_Lib_Sigma
use Symmetry_Info, only: Mul
use Orbital_Spaces, only: nAsh, nIsh, nOrb, nSym
use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iOpSpin, iCSym, iSSym, nInt1, nInt2, nC, nS
logical(kind=iwp), intent(in) :: Square, Have_2_el
real(kind=wp), intent(in) :: Int1(nInt1), Int2(nInt2), C(nC)
real(kind=wp), intent(out) :: S(nS)
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: iOff(8), iOp, iS, it, iu, jS, nA, nAOff(8), nTU, nTUVX, t, u
real(kind=wp), allocatable :: TU(:)

iOff(:) = 0
nAOff(:) = 0
iOp = Mul(iCSym,iSSym)
nA = 0
do iS=1,nSym
  if (iS > 1) iOff(iS) = iOff(iS-1)+nOrb(iS-1)*nOrb(Mul(iS-1,iOp))
  nAOff(iS) = nA
  nA = nA+nAsh(iS)
end do
if (iOff(nSym)+nOrb(nSym)*nOrb(Mul(nSym,iOp)) > nInt1) call Abend()
if (Square) then
  nTU = nA**2
  nTUVX = nTri_Elem(nA**2)
else
  nTU = nTri_Elem(nA)
  nTUVX = nTri_Elem(nTU)
end if
if (.not. Have_2_el) nTUVX = 1
if (nInt2 < nTUVX) call Abend()

call mma_allocate(TU,max(1,nTU),Label='TU')
TU(:) = Zero
do iS=1,nSym
  jS = Mul(iS,iOp)
  do iu=1,nAsh(jS)
    u = nAOff(jS)+iu
    do it=1,nAsh(iS)
      t = nAOff(iS)+it
      if (Square) then
        TU((u-1)*nA+t) = Int1(iOff(iS)+(nIsh(jS)+iu-1)*nOrb(iS)+nIsh(iS)+it)
      else if (t >= u) then
        TU(iTri(t,u)) = Int1(iOff(iS)+(nIsh(jS)+iu-1)*nOrb(iS)+nIsh(iS)+it)
      end if
    end do
  end do
end do

call CI_Lib_Sigma(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,max(1,nTU),TU,nTUVX,Int2(1:nTUVX),nC,C,nS,S)

call mma_deallocate(TU)

end subroutine MCLR_CI_Sigma
