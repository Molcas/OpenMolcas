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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************
      SUBROUTINE Tw_corr_drv(EOrb,nEO,CMO,nCMO,Ecorr)
      use InfSCF, only: nnOc, nSym, nOcc, nDel, nOrb, nFro, nBas
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nEO, nCMO
      Real*8 EOrb(nEO), CMO(nCMO), Ecorr

      Integer i, iOff, ipEOkk, ipEVir, iRC, iSym, jOff, jOkk, jOrb, jVir, kOff, nExt, nOkk
      Real*8, Allocatable :: Eov(:)

      Call mma_Allocate(Eov,nEO,Label='Eov')

      ipEOkk=1
      ipEVir=ipEOkk+nnOc
      iOff=0
      jOff=0
      kOff=0
      Do iSym=1,nSym
         nOkk=nFro(iSym)+nOcc(iSym,1)
         nExt=nBas(iSym)-nDel(iSym)-nOkk
         jOrb=1+jOff
         jOkk=ipEOkk+iOff
         Do i=0,nOkk-1
            Eov(jOkk+i)=EOrb(jOrb+i)
         End Do
         jOrb=jOrb+nOkk
         jVir=ipEVir+kOff
         Do i=0,nExt-1
            Eov(jVir+i)=EOrb(jOrb+i)
         End Do
         iOff=iOff+nOkk
         jOff=jOff+nOrb(iSym)
         kOff=kOff+nExt
      End Do

      Call Tw_corr(irc,Ecorr,CMO,Eov(:ipEVir-1),Eov(ipEVir:))

      Call mma_deallocate(Eov)
      Return
      End SUBROUTINE Tw_corr_drv
!****************************************************************************
!                                                                           *
!****************************************************************************
      SUBROUTINE Tw_corr(irc,DeTW,CMOI,EOcc,EVir)
      use InfSCF, only: nBT, nSym, nFro, nOcc, nDel, nBas
      use ChoMP2, only: ChoAlg, DoDens
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, Half
      Implicit None
      Integer iRC
      Real*8 DeTW, CMOI(*), EOcc(*), EVir(*)
!
      Integer nExt(8)
      Integer i, nElk
      Real*8 TW, TW0
      Real*8 Grad(1)
      Real*8, Allocatable :: DMAT(:,:), F_DFT(:)

      DoDens = .false.
      ChoAlg = 2
!
      CALL mma_allocate(DMAT,nBT,2,Label='DMAT')

      nElk=0
      Do i=1,nSym
         nExt(i)=nBas(i)-nDel(i)-nOcc(i,1)-nFro(i)
         nElk=nElk+2*(nFro(i)+nOcc(i,1))
      End Do

      CALL DM_FNO_RHF(irc,nSym,nBas,nFro,nOcc(1,1),nExt,nDel,CMOI,EOcc,EVir,DMAT(:,2),DMAT(:,1))
      If (irc .ne. 0) Then
         Write(6,*) 'DM_FNO_RHF returned ',irc
         Call SysAbendMsg('DM_FNO_RHF','Non-zero return code from DM_FNO_RHF',' ')
      EndIf

      CALL mma_allocate(F_DFT,nBT,Label='F_DFT')
!
      Call Fold_tMat(nSym,nBas,DMAT(:,1),DMAT(:,1))
      call dscal_(nBT,Half,DMAT(:,1),1)
      Call Fold_tMat(nSym,nBas,DMAT(:,2),DMAT(:,2))
      call dscal_(nBT,Half,DMAT(:,2),1)
      Grad=Zero

      Call wrap_DrvNQ('HUNTER',F_DFT,1,TW,DMAT(:,1),nBT,1,.false.,Grad,1,'SCF ')

      Call wrap_DrvNQ('HUNTER',F_DFT,1,TW0,DMAT(:,2),nBT,1,.false.,Grad,1,'SCF ')
      DeTW=(TW-TW0)/dble(nElk)
!
      Call mma_deallocate(F_DFT)
      Call mma_deallocate(DMAT)
!
      Return
      End SUBROUTINE Tw_corr
!****************************************************************************
!                                                                           *
!****************************************************************************

      SUBROUTINE DM_FNO_RHF(irc,nSym,nBas,nFro,nIsh,nSsh,nDel,CMOI,EOcc,EVir,DM0,DM)
!****************************************************************************
!                                                                           *
!     Purpose:  setup of FNO density matrix calculation (RHF-based)         *
!                                                                           *
!     Author:   F. Aquilante  (Geneva, Sep 2010)                            *
!                                                                           *
!****************************************************************************
      use ChoMP2, only: MP2_small
      use Constants, only: Zero, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
#include "Molcas.fh"
      Integer iRC, nSym
      Integer nBas(nSym),nFro(nSym),nIsh(nSym),nSsh(nSym),nDel(nSym)
      Real*8  CMOI(*), EOcc(*), EVir(*), DM0(*), DM(*)
!
      Integer i, ifr, ioff, iSkip, iSym, iTo, jD, jOcc, jp,jTo, jVir, kDM, kfr, kij, kOff, kTo, lij, lOff, nBasT, nBmx, nCMO, nOA, &
              nOkk, nOrb, nSQ, nTri, nVV, j, jOff
      Real*8 SqOcc, tmp, Dummy
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
      Real*8, Allocatable:: CMO(:,:), EOrb(:,:), DMAT(:)
!
      irc=0
      MP2_small=.false.
!
!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS, etc. AND CHECK LIMITS   *
!----------------------------------------------------------------------*
!
      nBasT=0
      ntri=0
      nSQ=0
      nBmx=0
      nOrb=0
      nVV=0
      Do i=1,nSym
        nBasT=nBasT+nBas(i)
        nOrb=nOrb+nFro(i)+nIsh(i)+nSsh(i)+nDel(i)
        ntri=ntri+nBas(i)*(nBas(i)+1)/2
        nSQ=nSQ+nBas(i)**2
        nVV=nVV+nSsh(i)**2
        nBmx=Max(nBmx,nBas(i))
      End Do
      IF(nBasT.GT.mxBas) then
       Write(6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
       Call Abend()
      Endif
!
      NCMO=nSQ
      Call mma_allocate(CMO,nCMO,2,Label='CMO')
      CALL DCOPY_(NCMO,CMOI,1,CMO(:,1),1)
!
      nOA=0
      Do iSym=1,nSym  ! setup info
         lnFro(iSym)=nFro(iSym)
         lnOcc(iSym)=nIsh(iSym)
         nOA=nOA+lnOcc(iSym)
         lnVir(iSym)=nSsh(iSym)
         lnOrb(iSym)=lnOcc(iSym)+lnVir(iSym)
         lnDel(iSym)=nDel(iSym)
      End Do
!
      Call mma_Allocate(EOrb,nOrb,4,Label='EOrb')
      jOff=0
      kOff=0
      lOff=0
      Do iSym=1,nSym
         jp=1+lOff+nFro(iSym)
         jOcc=jOff+1
         call dcopy_(nIsh(iSym),EOcc(jOcc),1,EOrb(jp,1),1)
         jVir=kOff+1
         jp=jp+nIsh(iSym)
         call dcopy_(nSsh(iSym),EVir(jVir),1,EOrb(jp,1),1)
         jOff=jOff+nIsh(iSym)
         kOff=kOff+nSsh(iSym)
         lOff=lOff+nBas(iSym)
      End Do
      ioff=0
      joff=0
      koff=0
      Do iSym=1,nSym
         ifr=1+ioff+nFro(iSym)
         ito=1+joff
         call dcopy_(nIsh(iSym),EOrb(ifr,1),1,EOrb(ito,3),1)
         ifr=1+ioff+nFro(iSym)+nIsh(iSym)
         ito=1+koff
         call dcopy_(nSsh(iSym),EOrb(ifr,1),1,EOrb(ito,4),1)
         ioff=ioff+nBas(iSym)
         joff=joff+nIsh(iSym)
         koff=koff+nSsh(iSym)
      End Do
      Call mma_Allocate(DMAT,nVV+nOA,Label='DMAT')
      DMAT(:)=Zero

      Call FnoSCF_putInf(nSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)

      CMO(:,2)=Zero
      iOff=0
      Do iSym=1,nSym
         kfr=1+iOff+nBas(iSym)*nFro(iSym)
         kto=1+iOff+nBas(iSym)*lnFro(iSym)
         call dcopy_(nBas(iSym)*lnOcc(iSym),CMO(kfr,1),1,CMO(kto,2),1)
         kfr=1+iOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
         kto=kto+nBas(iSym)*lnOcc(iSym)
         call dcopy_(nBas(iSym)*lnVir(iSym),CMO(kfr,1),1,CMO(kto,2),1)
         iOff=iOff+nBas(iSym)**2
      End Do
      Call Check_Amp_SCF(nSym,lnOcc,lnVir,iSkip)
      If (iSkip.gt.0) Then
         Call ChoMP2_Drv(irc,Dummy,CMO(:,2),EOrb(:,3),EOrb(:,4), &
                         DMAT(1:nVV),DMAT(nVV+1:))
         If(irc.ne.0) then
           write(6,*) 'MP2 pseudodensity calculation failed !'
           Call Abend
         Endif
      Else
         write(6,*)
         write(6,*)'There are ZERO amplitudes T(ai,bj) with the given '
         write(6,*)'combinations of inactive and virtual orbitals !! '
         write(6,*)'Check your input and rerun the calculation! Bye!!'
         Call Abend
      Endif
!
!
!     Compute the correlated density in AO basis
!     -------------------------------------------------------------
      jOcc=1+nVV
!           write(6,*) ' Occ    : ',(DMAT(jOcc+j),j=0,nOA-1)
!           write(6,*) ' Sum    : ',ddot_(nOA,One,0,DMAT(jOcc),1)
      call dscal_(nOA,Two,DMAT(jOcc),1)
      Call daxpy_(nOA,Two,[One],0,DMAT(jOcc),1)
!
      iOff=0
      jOff=0
      kDM=1
      Do iSym=1,nSym
!
         kto=1+jOff
         nOkk=nFro(iSym)+nIsh(iSym)
         Call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk,      &
                            Two,CMO(kto,1),nBas(iSym),         &
                                  CMO(kto,1),nBas(iSym),         &
                            Zero,DM0(kDM),nBas(iSym))
!
         sqocc=sqrt(Two)
         call dscal_(nBas(iSym)*nFro(iSym),sqocc,CMO(kto,1),1)
         Do j=0,nIsh(iSym)-1
             sqocc=sqrt(DMAT(jOcc+j))
             ito=kto+nBas(iSym)*j
             call dscal_(nBas(iSym),sqocc,CMO(ito,1),1)
         End Do
         Call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk,      &
                            One,CMO(kto,1),nBas(iSym),         &
                                  CMO(kto,1),nBas(iSym),         &
                            Zero,DM(kDM),nBas(iSym))
!
         if (nSsh(iSym).gt.0) then
           jD=1+iOff
!     Eigenvectors will be in increasing order of eigenvalues
           Call Eigen_Molcas(nSsh(iSym),DMAT(jD),EOrb(:,2),Eorb(:,1))
!     Reorder to get relevant eigenpairs first
           Do j=1,nSsh(iSym)/2
              Do i=1,nSsh(iSym)
                 lij=jD-1+nSsh(iSym)*(j-1)+i
                 kij=jD-1+nSsh(iSym)**2-(nSsh(iSym)*j-i)
                 tmp=DMAT(lij)
                 DMAT(lij)=DMAT(kij)
                 DMAT(kij)=tmp
              End Do
              tmp=EOrb(j,2)
              EOrb(j,2)=EOrb(nSsh(iSym)-j,2)
              EOrb(nSsh(iSym)-j,2)=tmp
           End Do
!
!     Compute new MO coeff. : X=C*U
           kfr=1+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
           kto=1+jOff+nBas(iSym)*(nFro(iSym)+nIsh(iSym))
           Call DGEMM_('N','N',nBas(iSym),nSsh(iSym),nSsh(iSym),  &
                              One,CMO(kfr,2),nBas(iSym),        &
                                    DMAT(jD),nSsh(iSym),          &
                              Zero,CMO(kto,1),nBas(iSym))

!          write(6,*) ' Occ_vir: ',(EOrb(j,2),j=1,nSsh(iSym))
!          write(6,*) ' Sum_vir: ',ddot_(nSsh(iSym),One,0,EOrb(:,2),1)
           Do j=0,nSsh(iSym)-1
              sqocc=sqrt(Two*EOrb(1+j,2))
              jto=kto+nBas(iSym)*j
              call dscal_(nBas(iSym),sqocc,CMO(jto,1),1)
           End Do
           Call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nSsh(iSym),       &
                              One,CMO(kto,1),nBas(iSym),                &
                                    CMO(kto,1),nBas(iSym),                &
                              One,DM(kDM),nBas(iSym))

           iOff=iOff+nSsh(iSym)**2
         endif
         jOff=jOff+nBas(iSym)**2
         kDM =kDM +nBas(iSym)*(nBas(iSym)+1)/2
         jOcc=jOcc+nIsh(iSym)
      End Do
!
      Call mma_deAllocate(EOrb)
      Call mma_deAllocate(DMAT)
      Call mma_deallocate(CMO)
!
      Return
      End SUBROUTINE DM_FNO_RHF
!***********************************************************************
!                                                                      *
!***********************************************************************
      Subroutine Check_Amp_SCF(nSym,nOcc,nVir,iSkip)

      Implicit None
      Integer nSym, nOcc(nSym), nVir(nSym), iSkip

      Integer nT1amTot, nT1am(8)
      Integer MulD2h, i, j, iSym, iSyma, iSymi

      MulD2h(i,j)=iEor(i-1,j-1) + 1

      iSkip=0
      nT1amTot=0
      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            nT1am(iSym) = nT1am(iSym) + nVir(iSyma)*nOcc(iSymi)
         End Do
         nT1amTot = nT1amTot + nT1am(iSym)
      End Do

      If (nT1amTot .gt. 0) iSkip=1
      Return
      End Subroutine Check_Amp_SCF
!***********************************************************************
!                                                                      *
!***********************************************************************
      SubRoutine FnoSCF_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
!
!     Purpose: put info in MP2 common blocks.
!
      Use ChoMP2, only: DoFNO, l_Dii
      Implicit None
      Integer mSym
      Integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
!
#include "corbinf.fh"
!
      Integer iSym
!
      nSym = mSym
!
      Do iSym = 1,nSym
         nOrb(iSym) = lnOrb(iSym)
         nOcc(iSym) = lnOcc(iSym)
         nFro(iSym) = lnFro(iSym)
         nDel(iSym) = lnDel(iSym)
         nExt(iSym) = lnVir(iSym)
      End Do
!
      DoFNO=.true.
      l_Dii=nOcc(1)
      Do iSym=2,nSym
         l_Dii=l_Dii+nOcc(iSym)
      End Do
!
      Return
      End SubRoutine FnoSCF_putInf
