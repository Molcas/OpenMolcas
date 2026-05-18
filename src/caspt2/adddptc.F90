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
      SUBROUTINE AddDPTC(nDPTC,nDSUM,DPTC,DSUM)

      use caspt2_module, only: NSYM, NFRO, NORB, NBAS
      use Constants, only: Half
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: nDPTC, nDSUM
      real(kind=wp), intent(in) :: DPTC(nDPTC)
      real(kind=wp), intent(inout) :: DSUM(nDSUM)

      integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, iOrb0,     &
     &  iOrb1, jOrb0, jOrb1, iOrb, jOrb
      real(kind=wp) :: Val

      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)!-nDel(iSym)
        If (nOrbI2 > 0) Then
          !! Add active orbital density
          !! Probably incorrect if symmetry
          Do iOrb0 = 1, nOrb(iSym)
            iOrb1 = nFro(iSym) + iOrb0
            Do jOrb0 = 1, nOrb(iSym)
              jOrb1 = nFro(iSym) + jOrb0
              DSUM(iMO1+iOrb0-1+nOrbI1*(jOrb0-1))                       &
     &          = DSUM(iMO1+iOrb0-1+nOrbI1*(jOrb0-1))                   &
     &          + DPTC(iMO2+iOrb1-1+nOrbI2*(jOrb1-1))
            End Do
          End Do
          !! Symmetrize DSUM
          Do iOrb = 1, nOrb(iSym)
            Do jOrb = 1, iOrb-1
              Val =(DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1))                   &
     &             +DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)))*Half
              DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1)) = Val
              DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)) = Val
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do

      End Subroutine AddDPTC
