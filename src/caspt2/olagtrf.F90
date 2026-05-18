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
! MO->AO or AO->MO transformation of 1-RDM
      Subroutine OLagTrf(mode,iSym,NBSQT,CMO,DPT2,DPT2AO,WRK)

      use caspt2_module, only: NFRO, NORB, NDEL, NBAS
      use Constants, only: Zero, One, Half
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: mode, iSym, NBSQT
      real(kind=wp), intent(in) :: CMO(NBSQT)
      real(kind=wp), intent(inout) :: DPT2(NBSQT), DPT2AO(NBSQT)
      real(kind=wp), intent(_OUT_) :: WRK(NBSQT)

      real(kind=wp) :: Val
      integer(kind=iwp) :: iCMO, iAO, iMO, jSym, nBasI, nOrbI, iBas,    &
     &  jBas

      !! Mode = 1: MO -> AO transformation
      !! Mode = 2: AO -> MO transformation
      iCMO =1
      iAO = 1
      iMO = 1
      Do jSym = 1, iSym-1
        iCMO = iCMO + nBas(jSym)**2 !! ??
        iAO  = iAO  + nBas(jSym)**2
        iMO  = iMO  + (nOrb(jSym)+nFro(jSym))**2
      End Do

      If (nOrb(iSym)+nFro(iSym) > 0) Then
        nBasI = nBas(iSym)
        nOrbI = nBas(iSym)-nDel(iSym)
        If (Mode == 1) Then
          !! MO -> AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,                        &
     &                One,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,              &
     &                Zero,WRK,nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,                        &
     &                One,WRK,nBasI,CMO(iCMO),nBasI,                    &
     &                Zero,DPT2AO(iAO),nBasI)
          !! Symmetrize, just in case
          Do iBas = 1, nBasI
            Do jBas = 1, iBas-1
              Val =(DPT2AO(iAO+iBas-1+nBasI*(jBas-1))                   &
     &            + DPT2AO(iAO+jBas-1+nBasI*(iBas-1)))*Half
              DPT2AO(iAO+iBas-1+nBasI*(jBas-1)) = Val
              DPT2AO(iAO+jBas-1+nBasI*(iBas-1)) = Val
            End Do
          End Do
        Else If (Mode == 2) Then
          !! AO -> MO
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,                        &
     &                One,CMO(iCMO),nBasI,DPT2AO(iAO),nBasI,            &
     &                Zero,WRK,nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,                        &
     &                One,WRK,nOrbI,CMO(iCMO),nBasI,                    &
     &                Zero,DPT2(iMO),nOrbI)
        End If
      END IF

      Return

      End Subroutine OLagTrf
