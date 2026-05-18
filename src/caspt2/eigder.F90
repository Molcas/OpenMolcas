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
      SUBROUTINE EigDer(NBSQT,nAshT,DPT2,DPT2C,FPT2AO,FPT2CAO,RDMEIG,   &
     &                  CMO,Trf,FPT2,FPT2C,FIFA,FIMO,RDMSA)

      use caspt2_global, only: OLag
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NFRO, NFROT, NISH, NASH,           &
     &                         NORB, NDEL, NBAS, NBAST
      use Constants, only: Zero, One, Two

      implicit none

      integer(kind=iwp), intent(in) :: NBSQT, nAshT
      real(kind=wp), intent(in) :: DPT2(NBSQT), DPT2C(NBSQT),           &
     &  FPT2AO(NBSQT), FPT2CAO(NBSQT), CMO(NBSQT), Trf(NBSQT),          &
     &  FPT2(NBSQT), FPT2C(NBSQT), FIFA(NBSQT), FIMO(NBSQT),            &
     &  RDMSA(nAshT**2)
      real(kind=wp), intent(inout) :: RDMEIG(nAshT**2)

      real(kind=wp),allocatable :: WRK1(:),FPT2_loc(:),FPT2C_loc(:),    &
     &                             RDMqc(:)

      integer(kind=iwp) :: iCMO, iAO, iSym, nBasI, nOrbI, iSQ, nFroI,   &
     &  nIshI, nAshI, nCor, iSQA, iT, iTabs, iU, iUabs, iTU, iTUA

      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(FPT2_loc,NBSQT,Label='FPT2_loc')
      call mma_allocate(FPT2C_loc,NBSQT,Label='FPT2C_loc')

      !! AO -> MO transformation
      iCMO =1
      iAO = 1
      if (nfrot /= 0) then
        FPT2_loc(1:NBSQT) = FPT2(1:NBSQT)
        FPT2C_loc(1:NBSQT) = FPT2C(1:NBSQT)
      else
        DO iSym = 1, nSym
          iCMO = iCMO  + nBas(iSym)*nFro(iSym)
!         iOFF = iWTMP + nBas(iSym)*nBas(iSym)
          If (nOrb(iSym) > 0) Then
            nBasI = nBas(iSym)
            nOrbI = nOrb(iSym)
            !! First, FPT2(AO) -> FPT2(MO)
            CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,                      &
     &                   One,CMO(iCMO),nBasI,FPT2AO(iAO),nBasI,         &
     &                   Zero,WRK1,nOrbI)
            CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,                      &
     &                   One,WRK1,nOrbI,CMO(iCMO),nBasI,                &
     &                   Zero,FPT2_loc(iAO),nOrbI)
            !! Second, FPT2C(AO) -> FPT2C(MO)
            CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,                      &
     &                   One,CMO(iCMO),nBasI,FPT2CAO(iAO),nBasI,        &
     &                   Zero,WRK1,nOrbI)
            CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,                      &
     &                   One,WRK1,nOrbI,CMO(iCMO),nBasI,                &
     &                   Zero,FPT2C_loc(iAO),nOrbI)
          END IF
          iCMO = iCMO + nBas(iSym)*(nOrb(iSym)+nDel(iSym))
          iAO  = iAO  + nBasI*nBasI
!         iMO  = iMO  + nBasI*nBasI
        End Do
      end if

      FPT2_loc(:)  = Two*FPT2_loc(:)
      FPT2C_loc(:) = Two*FPT2C_loc(:)
      !! construct Fock in MO

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
        !  Construct the active density of the orbital energy
        !  Assume the state-averaged density (SS- and XMS-CASPT2)
!       nSeq = 0
!       Call DCopy_(nAshI*nAshI,[Zero],0,WRK1,1)
!       Do iState = 1, nState
!         Wgt  = DWgt(iState,iState)
!         Wgt  = One/nState
!         Call DaXpY_(nDRef,Wgt,DMIX(:,iState),1,WRK1,1)
!       End Do
        !  RDM of CASSCF
        !  RDMSA is defined by a set of natural orbitals.
        !  Here, we have to transform to a set of quasi-canonical
        !  orbitals (RDMqc), so forward transformation is appropriate.
!       Call SQUARE(WRK1,RDMqc,1,nAshI,nAshI)
        RDMqc(1:nAshI**2) = RDMSA(1:nAshI**2)
        !! nbast?
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
!       If (nFroI == 0) Then
!         Call SQUARE(FIFA(iSQ),WRK1,1,nOrbI,nOrbI)
!       Else
!         Call OLagFroSq(iSym,NBSQT,FIFA(iSQ),WRK1)
!       End If
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,                          &
!    *              Two,WRK1,nOrbI,DPT2(iSQ),nOrbI,
     &              Two,FIFA(iSQ),nOrbI,DPT2(iSQ),nOrbI,                &
     &              One,OLAG,nOrbI)

        !! explicit derivative of the effective Hamiltonian
        !! dfpq/da = d/da(C_{mu p} C_{nu q} f_{mu nu})
        !!         = f_{mu nu}^a + (C_{mu m} U_{mp} C_{nu q}
        !!                       + C_{mu p} C_{nu m} U_{mq}) * f_{mu nu}
        !!         = f_{mu nu}^a + U_{mp} f_{mq} + U_{mq} f_{pm}
        !! U_{pq}  = f_{pm} df_{qm} + f_{mp} df_{mq}
        CALL DGEMM_('N','T',nOrbI,nOrbI,nOrbI,                          &
!    *              One,WRK1,nOrbI,DPT2C(iSQ),nOrbI,
     &              One,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,               &
     &              One,OLAG,nOrbI)
        CALL DGEMM_('T','N',nOrbI,nOrbI,nOrbI,                          &
!    *              One,WRK1,nOrbI,DPT2C(iSQ),nOrbI,
     &              One,FIMO(iSQ),nOrbI,DPT2C(iSQ),nOrbI,               &
     &              One,OLAG,nOrbI)
!       End If
        !! Implicit derivative of inactive orbitals (DPT2C)
        OLAG(iSQ:iSQ+nOrbI*nCor-1) = OLAG(iSQ:iSQ+nOrbI*nCor-1)         &
     &    + Two*FPT2C_loc(iSQ:iSQ+nOrbI*nCor-1)
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
            RDMEIG(iSQA+iTUA)                                           &
     &        = RDMEIG(iSQA+iTUA) + FPT2_loc(iSq+iTU)
!           write(u6,'(2i3,f20.10)') it,iu,FPT2(iSq+iTU)
          End Do
        End Do
        iSQ = iSQ + nOrbI*nOrbI
        iSQA= iSQA+ nAshI*nAshI
      End Do

      call mma_deallocate(WRK1)
      call mma_deallocate(FPT2_loc)
      call mma_deallocate(FPT2C_loc)

      END SUBROUTINE EigDer
