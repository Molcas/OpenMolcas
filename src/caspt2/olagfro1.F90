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
      Subroutine OLagFro1(NBSQT,nOLag,DPT2,OLag)

      use caspt2_global, only: FIFA_all
      use caspt2_module, only: NSYM, NFRO, NISH, NBAS, NDEL
      use Constants, only: Zero, Half
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: NBSQT, nOLag
      real(kind=wp), intent(inout) :: DPT2(NBSQT), OLag(nOLag)

      integer(kind=iwp) :: iMO, iSym, nOrbI, nFroI, nIshI, nBasI, iOrb, &
     &  jOrb
      real(kind=wp) :: Tmp

      iMO  = 1
      DO iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym)
        nFroI = nFro(iSym)
        If (nOrbI > 0 .and. nFroI > 0) Then
          nIshI = nIsh(iSym)
          nBasI = nBas(iSym)
          !! Make sure that the frozen orbital derivative is zero
          !! (it does not appear in the PT2 energy)
          OLag(1:nOrbI*nFroI) = Zero
          Do iOrb = 1, nFroI
            Do jOrb = nFroI+1, nFroI+nIshI
              Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))              &
     &                    -OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))             &
     &            /(FIFA_all(iOrb+nBasI*(iOrb-1))                       &
     &             -FIFA_all(jOrb+nBasI*(jOrb-1)))
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))                           &
     &          = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) + Tmp
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1))                           &
     &          = DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) + Tmp
            End Do
          End Do
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do
!     write(u6,*) 'DPT2 after frozen orbital'
!     call sqprt(dpt2,nbast)

      End Subroutine OLagFro1
