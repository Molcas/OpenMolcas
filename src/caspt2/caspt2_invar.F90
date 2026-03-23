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
!
! Put zero to wrong (incomplete) density matrix elements
! If the MRPT2 energy is non-invariant with respect to rotations among the inactive and secondary orbital spaces,
! we cannot determine the off-diagonal elements (in the iagonal block) as in TRDNS2D.
! Therefore, remove the off-diagonal elements. They are computed in caspt2_grad_invaria2 using orbital derivatives.
! The diagonal elements are correct.
!
  use Constants, only: Zero
  use definitions, only: iwp,wp
  use caspt2_module, only: nSym, nIsh, nAsh, nSsh, nOrb

  implicit none

  integer(kind=iwp), intent(in) :: NDPT2
  real(kind=wp), intent(inout) :: DPT2(NDPT2)
  integer(kind=iwp) :: IOFDIJ(8), IOFDAB(8)

  integer(kind=iwp) :: idij, is, ni, na, ns, no, idtu, idab, isym, ii, ij, ia, ib

  IDIJ = 0
  DO IS = 1,NSYM
    NI = NISH(IS)
    NA = NASH(IS)
    NO = NORB(IS)
    IDTU = IDIJ+NO*NI+NI
    IDAB = IDTU+NO*NA+NA
    IOFDIJ(IS) = IDIJ
    IOFDAB(IS) = IDAB
    IDIJ = IDIJ+NO*NO
  END DO

  do isym = 1, nsym
    NI = NISH(ISYM)
    NS = NSSH(ISYM)
    NO = NORB(ISYM)
    do ii = 1, ni
      do ij = 1, ni
        if (ii == ij) cycle
        IDIJ = IOFDIJ(ISYM)+II+NO*(IJ-1)
        DPT2(IDIJ) = Zero
      end do
    end do
    do ia = 1, ns
      do ib = 1, ns
        if (ia == ib) cycle
        IDAB = IOFDAB(ISYM)+IA+NO*(IB-1)
        DPT2(IDAB) = Zero
      end do
    end do
  end do

end subroutine caspt2_grad_invaria1
!
!-----------------------------------------------------------------------
!
subroutine caspt2_grad_invaria2(NDPT2,nOLag,DPT2,OLag)
!
! Compute pseudo-density in the inactive and secondary orbital blocks for MRPT2 methods that are not invariant with respect to
! orbital rotations among inactive and secondary orbitals using the canonical condition of MOs
! See the IPEA-shift implementation for the active block non-invariance, and the sigma^P implementation for the non-invariance
! with respect to rotations in the internally contracted basis
!
  use Constants, only: Half
  use definitions, only: iwp,wp
  use caspt2_module, only: nSym, nDel, nIsh, EPSI, nSsh, nAsh, EPSE, nBas, nFro

  implicit none

  integer(kind=iwp), intent(in) :: NDPT2, nOLag
  real(kind=wp), intent(inout) :: DPT2(NDPT2)
  real(kind=wp), intent(in)    :: OLag(nOLag)

  integer(kind=iwp) :: iMO, iSym, nOrbI, nFroI, nIshI, nAshI, nSshI, iOrb, jOrb
  real(kind=wp) :: Tmp

  iMO  = 1
  DO iSym = 1, nSym
    nOrbI = nBas(iSym)-nDel(iSym)
    nFroI = nFro(iSym)
    nIshI = nIsh(iSym)
    If (nIshI > 0) Then
      Do iOrb = nFroI+1, nFroI+nIshI
        Do jOrb = iOrb+1, nFroI+nIshI
          Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))-OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))/(EPSI(iOrb-nFroI)-EPSI(jOrb-nFroI))
          DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Tmp
          DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Tmp
        End Do
      End Do
    End If
    nSshI = nSsh(iSym)
    nAshI = nAsh(iSym)
    If (nSshI > 0) Then
      Do iOrb = nOrbI-nSshI+1, nOrbI
        Do jOrb = iOrb+1, nOrbI
          Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))-OLag(iMO+jOrb-1+nOrbI*(iOrb-1))) &
              /(EPSE(iOrb-nFroI-nIshI-nAshI)-EPSE(jOrb-nFroI-nIshI-nAshI))
          DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Tmp
          DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Tmp
        End Do
      End Do
    End If
    iMO  = iMO  + nOrbI*nOrbI
  End Do

end subroutine caspt2_grad_invaria2
