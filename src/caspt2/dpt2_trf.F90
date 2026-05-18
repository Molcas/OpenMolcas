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
      SUBROUTINE DPT2_Trf(NBSQT,nAshT,DPT2,DPT2AO,CMO,DEPSA,DSUM)

      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use Constants, only: Zero, One, Half
      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NORB,            &
     &                         NDEL, NBAS

      implicit none

      integer(kind=iwp), intent(in) :: NBSQT, nAshT
      real(kind=wp), intent(inout) :: DPT2(NBSQT), DPT2AO(NBSQT),       &
     &                                DSUM(NBSQT)
      real(kind=wp), intent(in) :: CMO(NBSQT), DEPSA(nAshT,nAshT)

      real(kind=wp), allocatable :: WRK(:)
      real(kind=wp) :: Val
      integer(kind=iwp) :: iCMO, iAO, iMO, iSym, nBasI, nOrbI, iOrb0,   &
     &  iOrb, jOrb0, jOrb

      !! DPT2 transformation
      !! Just transform DPT2 (in MO, block-squared) to DPT2AO (in AO,
      !! block-squared). Also, for DPT2C which couples with the inactive
      !! density matrix.
      call mma_allocate(WRK,NBSQT,Label='WRK')

      !! MO -> AO back transformation
      iCMO =1
      iAO = 1
      iMO = 1
      DO iSym = 1, nSym
        iCMO = iCMO  + nBas(iSym)*nFro(iSym)
        If (nORB(ISYM) > 0) Then
          nBasI = nBas(iSym)
          nOrbI = nOrb(iSym)
          !! Add active orbital density
          Do iOrb0 = 1, nAsh(iSym)
            iOrb = nIsh(iSym)+iOrb0
            ! iOrb2= nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              jOrb = nIsh(iSym)+jOrb0
              ! jOrb2= nFro(iSym)+nIsh(iSym)+jOrb0
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))                           &
     &          = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) + DEPSA(iOrb0,jOrb0)
              DSUM(iMO+iOrb-1+nOrbI*(jOrb-1))                           &
     &          = DSUM(iMO+iOrb-1+nOrbI*(jOrb-1)) + DEPSA(iOrb0,jOrb0)
            End Do
          End Do
          !! Symmetrize DPT2 (for shift)
          Do iOrb = 1, nOrb(iSym)
            Do jOrb = 1, iOrb
              Val =(DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))                     &
     &             +DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*Half
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
            End Do
          End Do
          !! First, DPT2 -> DPT2AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,                        &
     &                 One,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,             &
     &                 Zero,WRK,nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,                        &
     &                 One,WRK,nBasI,CMO(iCMO),nBasI,                   &
     &                 Zero,DPT2AO(iAO),nBasI)
        END IF
        iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
        iAO  = iAO  + nBasI*nBasI
        iMO  = iMO  + nBasI*nBasI
      End Do

      call mma_deallocate(WRK)

      END SUBROUTINE DPT2_Trf
