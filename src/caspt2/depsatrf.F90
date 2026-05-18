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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      Subroutine DEPSATrf(NBSQT,nAshT,DEPSA,FPT2,WRK1,WRK2)

      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NASH,          &
     &                         NBAS, NBAST
      use Constants, only: Zero, Half

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: NBSQT, nAshT
      real(kind=wp), intent(in) :: DEPSA(nAshT,nAshT)
      real(kind=wp), intent(_OUT_) :: FPT2(NBSQT), WRK1(NBSQT),         &
     &                                WRK2(NBSQT)

      real(kind=wp),allocatable :: DAO(:),DMO(:)
      integer(kind=iwp) :: iSym, iSymA, iSymI, iSymB, iSymJ,            &
     &  nCorI, nBasI, iAsh, jAsh, iAshI, iOrb, jAshI, jOrb

      FPT2(1:nBasT**2) = Zero

      iSym = 1
      iSymA= 1
      iSymI= 1
      iSymB= 1
      iSymJ= 1

!     If (nFroT /= 0.and.IfChol) Then
      If (IfChol) Then
        !! DEPSA(MO) -> DEPSA(AO) -> G(D) in AO -> G(D) in MO
        !! The Cholesky vectors do not contain frozen orbitals...
        call mma_allocate(DAO,NBSQT,Label='DAO')
        call mma_allocate(DMO,NBSQT,Label='DMO')
        !! First, MO-> AO transformation of DEPSA
        Do iSym = 1, nSym
          DMO(:) = Zero
          nCorI = nFro(iSym)+nIsh(iSym)
          nBasI = nBas(iSym)
          Do iAsh = 1, nAsh(iSym)
            Do jAsh = 1, nAsh(iSym)
              DMO(nCorI+iAsh+nBasI*(nCorI+jAsh-1)) = DEPSA(iAsh,jAsh)
            End Do
          End Do
          Call OLagTrf(1,iSym,NBSQT,CMOPT2,DMO,DAO,WRK1)
        End Do
        !! Compute G(D)
        WRK1(1:NBSQT) = Zero
        DMO(:) = Zero
        !! it's very inefficient
        Call OLagFro4(NBSQT,1,1,1,1,1,                                  &
     &                DAO,WRK1,DMO,WRK1,WRK2)
        !! G(D) in AO -> G(D) in MO
        Do iSym = 1, nSym
          Call OLagTrf(2,iSym,NBSQT,CMOPT2,FPT2,DMO,WRK1)
        End Do
        call mma_deallocate(DAO)
        call mma_deallocate(DMO)
      Else
        nCorI = nFro(iSym)+nIsh(iSym)
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI

            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            FPT2(1:nBasT**2) = FPT2(1:nBasT**2)                         &
     &        + DEPSA(iAshI,jAshI)*WRK1(1:nBasT**2)

            Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            FPT2(1:nBasT**2) = FPT2(1:nBasT**2)                         &
     &        - Half*DEPSA(iAshI,jAshI)*WRK1(1:nBasT**2)
          End Do
        End Do
      End If

      Return

      End Subroutine DEPSATrf
