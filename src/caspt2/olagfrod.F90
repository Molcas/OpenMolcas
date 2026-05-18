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
      Subroutine OLagFroD(NBSQT,NASHT,DIA,DI,RDMSA,Trf)

      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NBAS
      use Constants, only: Zero, One, Two

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: NBSQT, nAshT
      real(kind=wp), intent(_OUT_) :: DIA(NBSQT), DI(NBSQT)
      real(kind=wp), intent(in) :: RDMSA(nAshT**2), Trf(NBSQT)

      real(kind=wp),allocatable :: WRK1(:),WRK2(:)

      integer(kind=iwp) :: iAOtr, iAOsq, iSym, nFroI, nIshI, nAshI,     &
     &  nBasI, nCorI

      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')

      iAOtr = 0
      iAOsq = 1
      Do iSym = 1, nSym
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nBasI = nBas(iSym)
        nCorI = nFroI + nIshI                                           &

        !! full density matrix
      ! Call SQUARE(WRK1(1+iAOtr),DIA(iAOsq),1,nBasI,nBasI)
      ! !! off-diagonal elements have to be halved
      ! Do Mu = 1, nBasI
      !   Do Nu = 1, nBasI
      !     If (Mu == Nu) Cycle
      !     DIA(iAOsq+Mu-1+nBasI*(Nu-1))
     &!       = Half*DIA(iAOsq+Mu-1+nBasI*(Nu-1))
      !   End Do
      ! End Do

        !! inactive density matrix
        Call DGEMM_('N','T',nBasI,nBasI,nCorI,                          &
     &              Two,CMOPT2,nBasI,CMOPT2,nBasI,                      &
     &              Zero,DI(iAOsq),nBasI)

        !! inactive+active density matrix
        !! Somehow, the above density matrix obtained by calling
        !! Get_D1AO is incorrect... at least, cannot be used.
        ! 1) inactive part
        DIA(1:nBasI**2) = DI(1:nBasI**2)
        ! 2)  RDMSA is defined in CASSCF orbitals, so transform RDMSA to
        !     CASPT2 orbital basis
        Call DGemm_('T','N',nAshI,nAshI,nAshI,                          &
     &              One,Trf(1+nCorI+nBasI*nCorI),nBasI,RDMSA,nAshI,     &
     &              Zero,WRK2,nAshI)
        Call DGemm_('N','N',nAshI,nAshI,nAshI,                          &
     &              One,WRK2,nAshI,Trf(1+nCorI+nBasI*nCorI),nBasI,      &
     &              Zero,WRK1,nAshI)
        ! 3) Finally, add the active part
        Call DGemm_('N','N',nBasI,nAshI,nAshI,                          &
     &              One,CMOPT2(1+nBasI*nCorI),nBasI,WRK1,nAshI,         &
     &              Zero,WRK2,nBasI)
        Call DGemm_('N','T',nBasI,nBasI,nAshI,                          &
     &              One,WRK2,nBasI,CMOPT2(1+nBasI*nCorI),nBasI,         &
     &              One,DIA,nBasI)

        iAOtr = iAOtr + nBasI*(nBasI+1)/2
        iAOsq = iAOsq + nBasI*nBasI
      End Do

      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)

      Return

      End Subroutine OLagFroD
