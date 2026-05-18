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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine OLagNS2(iSym,NBSQT,lT2AO,DPT2C,T2AO)

use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: wp, iwp
use caspt2_module, only: NSYM, NACTEL, NFRO, NISH, NASH, NSSH, NDEL, NBAS

implicit none
integer(kind=iwp), intent(in) :: iSym, NBSQT, lT2AO
real(kind=wp), intent(inout) :: DPT2C(NBSQT), T2AO(lT2AO)
real(kind=wp), allocatable :: Int1(:), Scr1(:), Amp1(:)
integer(kind=iwp) :: nMaxOrb, jSym, lInt, iSymI, iSymJ, iSymIJ, iSymA, iSymB, iSymAB, iSymIJAB, iCase

!! orbital Lagrangian from the T-amplitude
!! See the loop structure in rhs_mp2
!! and the helper subroutine in rhs_mp2_help1/2

nMaxOrb = 0
do jSym=1,nSym
  nMaxOrb = max(nMaxOrb,nBas(jSym))
end do
lInt = nMaxOrb*nMaxOrb

call mma_allocate(Int1,lInt,Label='Int1') !! for (ia|jb)
call mma_allocate(Scr1,lInt,Label='Scr1') !! work space
call mma_allocate(Amp1,lInt,Label='Amp1') !! for amplitude

!! (ia|jb)
do iSymI=1,nSym !! Symmetry of occupied (docc+act) orbitals
  !! Check, in particular nFro
  if (nFro(iSymI)+nIsh(iSymI)+nAsh(iSymI) == 0) cycle
  do iSymJ=1,iSymI
    if (nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ) == 0) cycle
    iSymIJ = Mul(iSymI,iSymJ)
    do iSymA=1,nSym !! Symmetry of non-filled (act+virt) orbs
      if (nAsh(iSymA)+nSsh(iSymA)+nDel(iSymA) == 0) cycle
      do iSymB=1,iSymA
        if (nAsh(iSymB)+nSsh(iSymB)+nDel(iSymB) == 0) cycle
        iSymAB = Mul(iSymA,iSymB)
        iSymIJAB = Mul(iSymIJ,iSymAB)
        if (iSym /= iSymIJAB) cycle
        do iCase=1,13
          call OLagNs_Hel2(iCase,NBSQT,lT2AO,iSym,iSymA,iSymB,iSymI,iSymJ,nMaxOrb,Int1,Amp1,Scr1,DPT2C,T2AO)
        end do
      end do
    end do
  end do
end do

DPT2C(1:NBSQT) = DPT2C(1:NBSQT)/real(max(1,NACTEL),kind=wp)

call mma_deallocate(Int1)
call mma_deallocate(Scr1)
call mma_deallocate(Amp1)

end subroutine OLagNS2
