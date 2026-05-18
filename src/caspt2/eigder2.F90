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
      SUBROUTINE EigDer2(NBSQT,nAshT,RDMEIG,Trf,FIFA,RDMSA,DEPSA,WRK1,  &
     &                   WRK2)

      use caspt2_global, only: OLag
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NDEL, NBAS, NBAST
      use Constants, only: Zero, One, Two

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: NBSQT, nAshT
      real(kind=wp), intent(_OUT_) :: RDMEIG(NBSQT), WRK1(NBSQT),       &
     &                                WRK2(NBSQT)
      real(kind=wp), intent(in) :: Trf(NBSQT), FIFA(NBSQT),             &
     &                             RDMSA(nAshT**2), DEPSA(nAshT**2)

      real(kind=wp),allocatable :: FPT2_loc(:),RDMqc(:)
      integer(kind=iwp) :: iSQ, iSym, nOrbI, nFroI, nIshI, nAshI, nCor, &
     &  iSQA, iT, iTabs, iU, iUabs, iTU, iTUA

      call mma_allocate(FPT2_loc,NBSQT,Label='FPT2_loc')

      !! Compute G(D), where D=DEPSA
      Call DEPSATrf(NBSQT,nAshT,DEPSA,FPT2_loc,WRK1,WRK2)
      FPT2_loc(:) = Two*FPT2_loc(:)

      iSQ = 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        ! nSshI = nSsh(iSym)
        ! nDelI = nDel(iSym)
        nCor  = nFroI + nIshI
        !! Inactive orbital contributions: (p,q) = (all,inact)
        OLAG(iSQ:iSQ+nOrbI*nCor-1) = OLAG(iSQ:iSQ+nOrbI*nCor-1)         &
     &    + Two*FPT2_loc(iSQ:iSQ+nOrbI*nCor-1)
        !! Active orbital contributions: (p,q) = (all,act)
        call mma_allocate(RDMqc,nAshI**2,Label='RDMqc')
        RDMqc(1:nAshT**2) = RDMSA(1:nAshT**2)
        Call DGemm_('T','N',nAshT,nAshT,nAshT,                          &
     &              One,Trf(iSQ+nBasT*nCor+nCor),nBasT,RDMqc,nAshT,     &
     &              Zero,WRK1,nAshT)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,                          &
     &              One,WRK1,nAshT,Trf(iSQ+nBasT*nCor+nCor),nBasT,      &
     &              Zero,RDMqc,nAshT)
        !  Then just multiply with G(DPT2)
        CALL DGEMM_('N','N',nOrbI,nAshI,nAshI,                          &
     &              One,FPT2_loc(iSQ+nOrbI*nCor),nOrbI,RDMqc,nAshI,     &
     &              One,OLAG(iSQ+nOrbI*nCor),nOrbI)
        call mma_deallocate(RDMqc)
        !! From the third term of U_{ij}
        !  FIFA is already in quasi-canonical basis
        CALL DGEMM_('N','T',nOrbI,nAshI,nAshI,                          &
     &              Two,FIFA(1+nOrbI*nCor),nOrbI,DEPSA,nAshI,           &
     &              One,OLAG(iSQ+nOrbI*nCor),nOrbI)
        iSQ = iSQ + nOrbI*nOrbI
      End Do
!
!     ----- CASSCF density derivative contribution in active space
!
      iSQ = 1
      iSQA= 1
      Do iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym) !! nOrb(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nCor  = nFroI + nIshI
        Do iT = 1, nAshI
          iTabs = nCor + iT
          Do iU = 1, nAshI
            iUabs = nCor + iU
            iTU = iTabs-1 + nOrbI*(iUabs-1)
            iTUA= iT   -1 + nAshI*(iU   -1)
            RDMEIG(iSQA+iTUA) = FPT2_loc(iSq+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do

      call mma_deallocate(FPT2_loc)

      END SUBROUTINE EigDer2
