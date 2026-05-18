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
      Subroutine DPT2_TrfStore(Scal,NBSQT,DPT2q,DPT2n,Trf,WRK)

      use caspt2_module, only: NSYM, NORB, NDEL, NBAS
      use Constants, only: Zero, One
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: NBSQT
      real(kind=wp), intent(in) :: Scal, DPT2q(NBSQT), Trf(NBSQT)
      real(kind=wp), intent(inout) :: DPT2n(NBSQT)
      real(kind=wp), intent(_OUT_) :: WRK(NBSQT)

      integer(kind=iwp) :: iMO, iSym, nOrbI

      iMO = 1
      Do iSym = 1, nSym
        If (nOrb(iSym) > 0) Then
          nOrbI = nBas(iSym)-nDel(iSym)
          !! Quasi-canonical -> natural transformation of DPT2
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,                        &
     &                One,Trf(iMO),nOrbI,DPT2q(iMO),nOrbI,              &
     &                Zero,WRK,nOrbI)
          Call DGemm_('N','T',nOrbI,nOrbI,nOrbI,                        &
     &                Scal,WRK,nOrbI,Trf(iMO),nOrbI,                    &
     &                One,DPT2n(iMO),nOrbI)
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do

      Return

      End Subroutine DPT2_TrfStore
