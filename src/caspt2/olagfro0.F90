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
      Subroutine OLagFro0(NOSQT,NBSQT,DPT2_ori,DPT2)

      use caspt2_module, only: NSYM, NFRO, NORB, NDEL, NBAS
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: NOSQT, NBSQT
      real(kind=wp), intent(in) :: DPT2_ori(NOSQT)
      real(kind=wp), intent(inout) :: DPT2(NBSQT)

      integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, nFroI,     &
     &  iOrb, iOrb1, iOrb2, jOrb, jOrb1, jOrb2

      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI1 > 0) Then
          nFroI = nFro(iSym)
          !! Do for all orbitals
          Do iOrb = 1, nOrbI1
            iOrb1 = iOrb
            iOrb2 = iOrb+nFroI
            Do jOrb = 1, nOrbI1
              jOrb1 = jOrb
              jOrb2 = jOrb+nFroI
              DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))                       &
     &          = DPT2_ori(iMO1+iOrb1-1+nOrbI1*(jOrb1-1))
              DPT2(iMO2+jOrb2-1+nOrbI2*(iOrb2-1))                       &
     &          = DPT2_ori(iMO1+jOrb1-1+nOrbI1*(iOrb1-1))
            End Do
          End Do
        End If
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do

      End Subroutine OLagFro0
