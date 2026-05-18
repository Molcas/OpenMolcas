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
      Subroutine OLagFro3(NBSQT,FIFA,FIMO,WRK1,WRK2)

      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NDEL, NBAS, NBTRI

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: NBSQT
      real(kind=wp), intent(inout) :: FIFA(NBSQT), FIMO(NBSQT)
      real(kind=wp), intent(_OUT_) :: WRK1(NBSQT), WRK2(NBSQT)

      Character(Len=8) :: Label
      real(kind=wp), allocatable :: WFLT(:)
      integer(kind=iwp) :: IRC, IOPT, ICOMP, ISYLBL, iAO, iAOtr, iCMO,  &
     &  iMO, iSym, nBasI, nOrbI

      !! Read H_{\mu \nu}
      call mma_allocate(WFLT,NBTRI,Label='WFLT')
      IRC=-1
      IOPT=6
      ICOMP=1
      ISYLBL=1
      Label='OneHam  '
      CALL RDONE(IRC,IOPT,Label,ICOMP,WFLT,ISYLBL)

      !! AO -> MO transformation
      iAO   = 1
      iAOtr = 1
      iCMO  = 1
      iMO   = 1
      DO iSym = 1, nSym
        nBasI = nBas(iSym)
        nOrbI = nBas(iSym)-nDel(iSym)

        !! FIFA
        !! WRK1 = G(D)
        WRK1(1:nBasI*nBasI) = FIFA(iAO:iAO+nBasI*nBasI-1)
        !! WRK1 = H+G(D)
        Call Square(WFLT(iAOtr),WRK2,1,nBasI,nBasI)
        WRK1(1:nBasI*nBasI) = WRK1(1:nBasI*nBasI) + WRK2(1:nBasI*nBasI)
        !! AO -> MO transformation of H+G(D)
        Call OLagTrf(2,iSym,NBSQT,CMOPT2(iCMO),FIFA(iMO),WRK1,WRK2)

        !! FIMO
        !! WRK1 = G(D)
        WRK1(1:nBasI*nBasI) = FIMO(iAO:iAO+nBasI*nBasI-1)
        !! WRK1 = H+G(D)
        Call Square(WFLT(iAOtr),WRK2,1,nBasI,nBasI)
        WRK1(1:nBasI*nBasI) = WRK1(1:nBasI*nBasI) + WRK2(1:nBasI*nBasI)
        !! AO -> MO transformation of H+G(D)
        Call OLagTrf(2,iSym,NBSQT,CMOPT2(iCMO),FIMO(iMO),WRK1,WRK2)

        iAO   = iAO   + nBasI*nBasI
        iAOtr = iAOtr + nBasI*(nBasI+1)/2
        iCMO  = iCMO  + nBasI*nOrbI !?
        iMO   = iMO   + nOrbI*nOrbI
      End Do
!     write(u6,*) 'FIFA'
!     call sqprt(fifa,nbast)
!     write(u6,*) 'FIMO'
!     call sqprt(fimo,nbast)

      call mma_deallocate(WFLT)

      End Subroutine OLagFro3
