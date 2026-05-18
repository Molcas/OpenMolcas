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
      Subroutine OLagFro2(NBSQT,DPT2,FPT2,ERI,Scr)

      use caspt2_module, only: NSYM, NFRO, NISH, NDEL, NBAS
      use Constants, only: Half
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: NBSQT
      real(kind=wp), intent(in) :: DPT2(NBSQT)
      real(kind=wp), intent(inout) :: FPT2(NBSQT)
      real(kind=wp), intent(_OUT_) :: ERI(NBSQT), Scr(NBSQT)

      integer(kind=iwp) :: iMO, iSymI, iSymJ, iSymA, iSymB, iSym, nOrbI,&
     &  nFroI, nIshI, iOrb, jOrb
      real(kind=wp) :: Scal, Val

!     write(u6,*) 'FPT2 before frozen orbital'
!     call sqprt(fpt2,nbast)
      iMO = 1
      isymi = 1
      isymj = 1
      isyma = 1
      isymb = 1
      DO iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        !! Fpq = ((pq|rs)-1/2(pr|qs))*Drs
        Do iOrb = 1, nFroI
          Do jOrb = nFroI+1, nFroI+nIshI
            Scal = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,ERI,Scr)
            FPT2(iMO:iMO+nOrbI*nOrbI-1) = FPT2(iMO:iMO+nOrbI*nOrbI-1)   &
     &        + Scal*ERI(1:nOrbI*nOrbI)
            Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,ERI,Scr)
            FPT2(iMO:iMO+nOrbI*nOrbI-1) = FPT2(iMO:iMO+nOrbI*nOrbI-1)   &
     &        - Half*Scal*ERI(1:nOrbI*nOrbI)
          End Do
        End Do

        !! Symmetrize FPT2
        Do iOrb = 1, nOrbI
          Do jOrb = 1, iOrb-1
            Val = (FPT2(iMO+iOrb-1+nOrbI*(jOrb-1))                      &
     &            +FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*Half
            FPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
            FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
          End Do
        End Do
        iMO = iMO + nOrbI*nOrbI
      End Do
!     write(u6,*) 'FPT2 after frozen orbital'
!     call sqprt(fpt2,nbast)

      End Subroutine OLagFro2
