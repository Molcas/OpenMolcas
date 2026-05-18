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
      SUBROUTINE AddDEPSA(nDPT2,nAshT,DPT2,DEPSA)

      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NORB,            &
     &                         NDEL, NBAS
      use Constants, only: Half
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: nDPT2, nAshT
      real(kind=wp), intent(inout) :: DPT2(nDPT2)
      real(kind=wp), intent(in) :: DEPSA(nAshT,nAshT)

      integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, iOrb0,     &
     &  iOrb2, jOrb0, jOrb2, iOrb, jOrb
      real(kind=wp) :: Val

      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI2 > 0) Then
          !! Add active orbital density
          !! Probably incorrect if symmetry
          Do iOrb0 = 1, nAsh(iSym)
            ! iOrb1 = nIsh(iSym)+iOrb0
            iOrb2 = nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              ! jOrb1 = nIsh(iSym)+jOrb0
              jOrb2 = nFro(iSym)+nIsh(iSym)+jOrb0
              DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))                       &
     &          = DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))                   &
     &          + DEPSA(iOrb0,jOrb0)
            End Do
          End Do
          !! Symmetrize DPT2 (for shift)
          Do iOrb = 1, nBas(iSym)-nDel(iSym)
            Do jOrb = 1, iOrb-1
              Val =(DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1))                   &
     &             +DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)))*Half
              DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1)) = Val
              DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do

      End Subroutine AddDEPSA
