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
      SUBROUTINE OLagNS2(iSym,NBSQT,lT2AO,DPT2C,T2AO)

      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NACTEL, NFRO, NISH, NASH, NSSH,    &
     &                         NDEL, NBAS

      implicit none

      integer(kind=iwp), intent(in) :: iSym, NBSQT, lT2AO
      real(kind=wp), intent(inout) :: DPT2C(NBSQT), T2AO(lT2AO)

      real(kind=wp), allocatable :: Int1(:), Scr1(:), Amp1(:)
      integer(kind=iwp) :: nMaxOrb, jSym, lInt, iSymI, iSymJ, iSymIJ,   &
     &  iSymA, iSymB, iSymAB, iSymIJAB, iCase

      !! orbital Lagrangian from the T-amplitude
      !! See the loop structure in rhs_mp2.f
      !! and the helper subroutine in rhs_mp2_help1/2

      nMaxOrb=0
      Do jSym = 1, nSym
        nMaxOrb = Max(nMaxOrb,nBas(jSym))
      End Do
      lInt = nMaxOrb*nMaxOrb

      call mma_allocate(Int1,lInt,Label='Int1') !! for (ia|jb)
      call mma_allocate(Scr1,lInt,Label='Scr1') !! work space
      call mma_allocate(Amp1,lInt,Label='Amp1') !! for amplitude

      !! (ia|jb)
      Do iSymI = 1, nSym !! Symmetry of occupied (docc+act) orbitals
        !! Check, in particular nFro
        If (nFro(iSymI)+nIsh(iSymI)+nAsh(iSymI) == 0) Cycle
        Do iSymJ = 1, iSymI
          If (nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ) == 0) Cycle
          iSymIJ = 1 + iEor(iSymI-1,iSymJ-1)
          Do iSymA = 1, nSym !! Symmetry of non-filled (act+virt) orbs
            If (nAsh(iSymA)+nSsh(iSymA)+nDel(iSymA) == 0) Cycle
            Do iSymB = 1, iSymA
              If (nAsh(iSymB)+nSsh(iSymB)+nDel(iSymB) == 0) Cycle
              iSymAB = 1 + iEor(iSymA-1,iSymB-1)
              iSymIJAB = 1 + iEor(iSymIJ-1,iSymAB-1)
              If (iSym /= iSymIJAB) Cycle
              Do iCase = 1, 13
                Call OLagNs_Hel2(iCase,NBSQT,lT2AO,iSym,iSymA,iSymB,    &
     &                           iSymI,iSymJ,nMaxOrb,Int1,Amp1,Scr1,    &
     &                           DPT2C,T2AO)
              End Do
            End Do
          End Do
        End Do
      End Do

      DPT2C(1:NBSQT) = DPT2C(1:NBSQT)/real(max(1,NACTEL),kind=wp)

      call mma_deallocate(Int1)
      call mma_deallocate(Scr1)
      call mma_deallocate(Amp1)

      END SUBROUTINE OLagNS2
