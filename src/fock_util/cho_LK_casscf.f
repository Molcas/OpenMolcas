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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
      SUBROUTINE CHO_LK_CASSCF(DLT,FLT,MSQ,W_PWXY,FactXI,nFIorb,nAorb,  &
     &                         nChM,Ash,DoActive,nScreen,dmpk,dFmat,    &
     &                         CMO,ExFac)

!*********************************************************************
!  Author : F. Aquilante
!
!           This routine makes use of the "Local K" scheme for
!                computing the exchange terms in both the Inactive
!                and Active Fock matrix
!
! *************** INACTIVE AO-BASIS FOCK MATRIX **********************
!
!   FI(ab) = 2 * sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
!
! ***************   ACTIVE AO-BASIS FOCK MATRIX **********************
!
!   FA(ab) = sum_J  Lab,J * U(J)  -  0.5 * sum_Jw  Lwa,J * Lwb,J
!
! ***************   (WA|XY) integrals     ****************************
!
!   (WA|XY) = sum_J  L(wa,J) * L(xy,J)
!
!*********************************************************************
!
!      V(J) = sum_gd  Lgd,J * DI(gd)     DI=DLT(1)
!      U(J) = sum_gd  Lgd,J * DA(gd)     DA=DLT(2)
!
!      a,b,g,d:  AO-index
!      k:        MO-index   belonging to (Frozen+Inactive)
!      u,w,x,y:  MO-indeces belonging to (Active)
!
!*********************************************************************

#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      use ChoArr, only: nBasSh, nDimRS
      use ChoSwp, only: nnBstRSh, InfVec, IndRed
      use Data_Structures, only: SBA_Type, Allocate_SBA,                &
     &                           Deallocate_SBA
      use Data_Structures, only: DSBA_Type, Allocate_DSBA,              &
     &                           Deallocate_DSBA
      use Data_Structures, only: twxy_Type
      use Data_Structures, only: Allocate_twxy, Deallocate_twxy
      use Data_Structures, only: NDSBA_Type, Allocate_NDSBA,            &
     &                           Deallocate_NDSBA
      use Data_Structures, only: Allocate_L_Full, Deallocate_L_Full
      use Data_Structures, only: L_Full_Type
      use Data_Structures, only: Allocate_Lab, Deallocate_Lab,          &
     &                           Lab_Type
      Implicit Real*8 (a-h,o-z)

      Real*8 W_PWXY(*)
      Integer   kOff(8,2),nnA(8,8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2)
      Real*8    tmotr(2),tscrn(2)

      Type (DSBA_Type)   MSQ, CMO
      Type (DSBA_Type)   FLT(2), DLT(2), Ash(2)
      Type (SBA_Type) Laq(1), Lxy
      Type (twxy_Type) Scr
      Type (NDSBA_Type) DIAH
      Type (L_Full_Type) L_Full
      Type (Lab_Type) Lab

      Integer   nFIorb(8),nAorb(8),nChM(8)
#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   DoTraInt,DoActive,DoScreen
      Real*8    FactC(2),FactX(2),tau(2),xtau(2),thrv(2)
      Real*8    dmpK, dFmat,ExFac
      Character*50 CFmt
      Character(LEN=13), Parameter:: SECNAM = 'CHO_LK_CASSCF'
#include "chotime.fh"
#include "choscreen.fh"

      Logical, Parameter:: DoRead = .false.

#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Real*8, parameter:: xone = -one

      Type (DSBA_Type) :: KLT(2)
      Logical add
      Character(LEN=6) mode
      Real*8   LKThr
      Real*8, External :: Cho_LK_ScreeningThreshold
      Integer, External :: Cho_LK_MaxVecPerBatch

      Real*8, Allocatable:: Lrs(:,:), Drs(:,:), Frs(:,:)
      Real*8, Allocatable:: VJ(:)

      Integer, Allocatable:: nnBfShp(:,:), kOffSh(:,:),                 &
     &                       iShp_rs(:), Indx(:,:)
      Real*8, Allocatable:: SvShp(:,:), Diag(:), AbsC(:), SumAClk(:,:), &
     &                      Ylk(:,:), MLk(:,:), Faa(:), Fia(:)
#if defined (_MOLCAS_MPP_)
      Real*8, Allocatable:: DiagJ(:)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface

        Subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,    &
     &                         iSwap,IREDC,nDen,kDen,MOs,ChoT,          &
     &                         DoRead)
        use Data_Structures, only: DSBA_Type, SBA_Type
        Integer irc, lRedVec
        Real*8 RedVec(lRedVec)
        Integer IVEC1,NUMV,ISYM,iSwap,IREDC
        Integer   nDen,kDen

        Type (DSBA_Type) MOs(nDen)
        Type (SBA_Type) Chot(nDen)

        Logical   DoRead
        End Subroutine Cho_X_getVtra

        subroutine dgemv_(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
          Character(LEN=1) TRANS
          Integer M, N
          Real*8 ALPHA, BETA
          Integer LDA, INCX, INCY
          Real*8  A(lda,*), X(*), Y(*)
        End subroutine dgemv_

      End Interface
!                                                                      *
!***********************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
!*****
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
!***********************************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in CASSCF-debug
#endif

      DoTraInt = .false.
      IREDC = -1  ! unknown reduced set in core

      FactC(:) = [ one, one ]
      FactX(:) = [ FactXI*ExFac, -0.5D0*ExFac ]

      nDen=2  ! inactive and active density, respectively

      If (.not.DoActive) nDen=1

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      ! 1 --> CPU   2 --> Wall
      tread(:) = zero  !time read/transform vectors
      tcoul(:) = zero  !time for computing Coulomb
      texch(:) = zero  !time for computing Exchange
      tintg(:) = zero  !time for computing (tw|xy) integrals
      tmotr(:) = zero  !time for the half-transf of vectors
      tscrn(:) = zero  !time for screening overhead

! ==================================================================
      Call set_nnA(nSym,nAorb,nnA)

! --- Various offsets
! --------------------
      MaxB=nBas(1)
      DO ISYM=2,NSYM
        MaxB=Max(MaxB,nBas(iSym))
      END DO

!*************************************************

      nnO=0
      nIt=0
      DO jDen=1,nDen

         kOff(1,jDen) = nnO

         DO ISYM=2,NSYM

            nnO = nnO + nFIorb(iSym-1)*(2-jDen)                         &
     &                + nChM(iSym-1)*(jDen-1)

            kOff(iSym,jDen) = nnO

         END DO

         nnO = nnO + nFIorb(nSym)*(2-jDen)                              &
     &             + nChM(nSym)*(jDen-1)

         nIt = nIt + nnO*(2-jDen) ! tot # of (inactive+frozen) orbitals

      END DO

      nAt = nnO - nIt  ! tot # of active orb (decomp. Active density)

! --- Define max number of vectors to be treated in core at once

      MaxVecPerBatch=Cho_LK_MaxVecPerBatch()

! --- Define the basic screening threshold

      LKThr=Cho_LK_ScreeningThreshold(dFmat)

! --- Adjust the damping according to the Abs(Max BLB)
!tbp, may 2013: adjustment moved to Cho_LK_ScreeningThreshold
      fcorr = dmpK
!tbp  If (dFmat.gt.zero) Then
!tbp     If (dFmat.lt.1.0d3*LKThr) Then
!tbp        fcorr=dmpk*1.0d-2
!tbp     EndIf
!tbp     If (dFmat.le.LKThr) fcorr=fcorr*1.0d-2
!tbp  EndIf

      tau(1) = (LKThr/Max(1,nIt))*fcorr !Inactive fock matrix screening
      tau(2) = (LKThr/Max(1,nAt))*fcorr !Active fock matrix screening

      MaxRedT=MaxRed
      Call GAIGOP_SCAL(MaxRedT,'+')

      If (Estimate) Then
         Do i=1,nDen
            tau(i)=tau(i)/MaxRedT
         End Do
      EndIf

      xtau(1) = sqrt(tau(1))
      xtau(2) = sqrt(tau(2))

! --- Vector MO transformation screening thresholds
      NumVT=NumChT
      Call GAIGOP_SCAL(NumVT,'+')
      thrv(1) = ( sqrt(LKThr/(Max(1,nIt)*NumVT)) )*fcorr
      thrv(2) = ( sqrt(LKThr/(Max(1,nAt)*NumVT)) )*fcorr

      CALL mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')

#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par()) Then
         NNBSTMX=0
         Do i=1,nSym
            NNBSTMX = Max(NNBSTMX,NNBSTR(i,1))
         End Do
         Call mma_allocate(diagJ,NNBSTMX,Label='diagJ')
         diagJ(:)=Zero
      EndIf
#endif
      Do jDen = 1, nDen
         Call Allocate_DSBA(KLT(jDen),nBas,nBas,nSym,aCase='TRI')
         KLT(jDen)%A0(:)=Zero
      End Do

! *************** Read the diagonal integrals (stored as 1st red set)
      If (Update) CALL CHO_IODIAG(DIAG,2) ! 2 means "read"

! --- allocate memory for sqrt(D(a,b)) stored in full (squared) dim
      Call Allocate_NDSBA(DIAH,nBas,nBas,nSym)
      DIAH%A0(:)=Zero

! --- allocate memory for the abs(C(l)[k])
      Call mma_allocate(AbsC,MaxB,Label='AbsC')

! --- allocate memory for the Y(l)[k] vectors
      Call mma_allocate(Ylk,MaxB,nnO,Label='Ylk')

! --- allocate memory for the ML[k] lists of largest elements
! --- in significant shells
      Call mma_allocate(MLk,nShell,nnO,Label='MLk')

! --- allocate memory for the list of  S:= sum_l abs(C(l)[k])
! --- for each shell
      Call mma_allocate(SumAClk,nShell,nnO,Label='SumAClk')

! --- allocate memory for the Index arrays
      Call mma_allocate(Indx,[0,nShell],[1,nnO],Label='Indx')

! --- allocate memory for kOffSh
      Call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')

! --- allocate memory for nnBfShp
      Call mma_allocate(nnBfShp,nnShl_tot,nSym,Label='nnBfShp')

! --- allocate memory for iShp_rs
      Call mma_allocate(iShp_rs,nnShl_tot,Label='iShp_rs')

! --- allocate memory for the shell-pair Frobenius norm of the vectors
      Call mma_allocate(SvShp,nnShl,2,Label='SvShp')


! *** Compute Shell Offsets ( MOs and transformed vectors)

      MxBasSh = 0

      Do iSyma=1,nSym

         LKsh=0

         Do iaSh=1,nShell    ! kOffSh(iSh,iSym)

            kOffSh(iaSh,iSyma) = LKsh

            LKsh = LKsh + nBasSh(iSyma,iaSh)

            MxBasSh = Max(MxBasSh,nBasSh(iSyma,iaSh))

         End Do

      End Do

! --- allocate memory for the Diagonal of the Fock matrix
      Call mma_allocate(Fia,MxBasSh,Label='Fia')
      Call mma_allocate(Faa,nShell,Label='Faa')
      Faa(:)=Zero
      Fia(:)=Zero

! *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
      Do jDen=1,nDen
         Do kSym=1,nSym

            If (jDen.eq.2) Then
               Do jK=1,nChM(kSym)
                  jK_a = jK + kOff(kSym,jDen)

               Do iaSh=1,nShell


                  SKsh=zero
                  iS = kOffSh(iaSh,kSym) + 1
                  iE = kOffSh(iaSh,kSym) + nBasSh(kSym,iaSh)
                  Do ik=iS, iE
                     SKsh = SKsh + Ash(2)%SB(kSym)%A2(ik,jK)**2
                  End Do

                  SumAClk(iaSh,jk_a) = SKsh

               End Do
               End Do

            Else
               Do jK=1,nFIorb(kSym)
                  jK_a = jK + kOff(kSym,jDen)


               Do iaSh=1,nShell


                  SKsh=zero
                  iS = kOffSh(iaSh,kSym) + 1
                  iE = kOffSh(iaSh,kSym) + nBasSh(kSym,iaSh)
                  Do ik=iS, iE
                     SKsh = SKsh + MSQ%SB(kSym)%A2(ik,jK)**2
                  End Do

                  SumAClk(iaSh,jk_a) = SKsh

               End Do
               End Do

            End If

         End Do
      End Do

! *** Compute Shell-pair Offsets in the K-matrix

         Do iSyma=1,nSym

            LKshp=0

            Do iaSh=1,nShell

             Do ibSh=1,iaSh

               iShp = iaSh*(iaSh-1)/2 + ibSh

               nnBfShp(iShp,iSyma) = LKShp

               LKShp = LKShp + nBasSh(iSyma,iaSh)*nBasSh(iSyma,ibSh)    &
     &                       - (1-Min((iaSh-ibSh),1))*nBasSh(iSyma,iaSh)&
     &                       * (nBasSh(iSyma,iaSh) - 1)/2

             End Do

            End Do

         End Do

! *** Mapping shell pairs from the full to the reduced set

      Call Mk_iShp_rs(iShp_rs,nShell)

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

        NumCV=NumCho(jSym)
        Call GAIGOP_SCAL(NumCV,'max')
        If (NumCV .lt. 1) Cycle

        JNUM=1
        Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,      &
     &                       Memory=LFULL)

        iCase = 1  ! (wa|xy)
        Call Allocate_twxy(Scr,nAorb,nBas,JSYM,nSym,iCase)

! ****************     MEMORY MANAGEMENT SECTION    *****************
! ------------------------------------------------------------------
! --- compute memory needed to store at least 1 vector of JSYM
! --- and do all the subsequent calculations
! ------------------------------------------------------------------
         mTvec = 0
         mTvec1=0
         mTvec2=0
         MxB=0
         do l=1,nSym
            k=Muld2h(l,JSYM)
            If ((nFIorb(k)+nChM(k)).gt.0) MxB=Max(MxB,nBas(l))
            mTvec1= mTvec1+ nAorb(k)*nBas(l)
            If (k.le.l) mTvec2= mTvec2+ nnA(k,l)
         end do
         mTvec = mTvec1 + mTvec2

         LFMAX = Max(mTvec,LFULL) ! re-use memory for the active vec
         mTvec = Max(MxB,1) ! mem for storing half-transformed vec

! ------------------------------------------------------------------
! ------------------------------------------------------------------

         iLoc = 3 ! use scratch location in reduced index arrays

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
#if defined (_MOLCAS_MPP_)
         myJRED1=JRED1 ! first red set present on this node
         ntv0=0
#endif
         myJRED2=JRED2 ! last  red set present on this node

! --- entire red sets range for parallel run
         Call GAIGOP_SCAL(JRED1,'min')
         Call GAIGOP_SCAL(JRED2,'max')

         kscreen=1
         DoScreen=.true.

         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
               call Abend
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED)
!           !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero rc = ',irc
              call Abend
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then
             Call mma_allocate(Drs,nRS,nDen,Label='Drs')
             Call mma_allocate(Frs,nRS,nDen,Label='Frs')
             Drs(:,:)=Zero
             Frs(:,:)=Zero
            EndIf

            Call mma_maxDBLE(LWORK)

            nVec=min(LWORK/(nRS+mTvec+LFMAX),min(nVrs,MaxVecPerBatch))

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec+LFMAX
               WRITE(6,*) 'jsym= ',jsym
               WRITE(6,*) ' nRS = ',nRS
               WRITE(6,*) ' mTvec = ',mTvec
               WRITE(6,*) ' LFMAX = ',LFMAX
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

            If(JSYM.eq.1)Then
! --- Transform the densities to reduced set storage
               mode = 'toreds'
               add  = .false.
               Call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,                &
     &                           DLT,Drs,mode,add)
            EndIf

! --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch

               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               else
                  JNUM = nVec
               endif

               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CWTIME(TCR1,TWR1)

               CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,                &
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               If(JSYM.eq.1)Then

                 CALL CWTIME(TCC1,TWC1)

                 Call mma_allocate(VJ,JNUM,Label='VJ')

                 Fact = dble(min(jVec-iVrs,1))

                 Do jDen=1,nDen
! ************ COULOMB CONTRIBUTIONS  *********************
!
!  jDen=1   ---> inactive fock matrix
!  jDen=2   ---> active fock matrix
!
! --- Contraction with the density matrix
! ---------------------------------------
! --- V{#J} =  sum_rs  L(rs,{#J}) * DI(rs)
!==========================================================
!
                   CALL DGEMV_('T',nRS,JNUM,                            &
     &                  ONE,Lrs,nRS,                                    &
     &                  Drs(:,jDen),1,ZERO,VJ,1)

! --- F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
!===============================================================

                   CALL DGEMV_('N',nRS,JNUM,                            &
     &                 FactC(jDen),Lrs,nRS,                             &
     &                 VJ,1,Fact,Frs(:,jDen),1)

                 End Do
                 Call mma_deallocate(VJ)

                 CALL CWTIME(TCC2,TWC2)
                 tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                 tcoul(2) = tcoul(2) + (TWC2 - TWC1)


               EndIf  ! Coulomb contribution


! *************** EXCHANGE CONTRIBUTIONS  ***********************

               CALL CWTIME(TCS1,TWS1)
! ---------------------------------------------------------------------
! --- Estimate the diagonals :   D(a,b) = sum_J (Lab,J)^2
!
               If (Estimate) Then

                  Call Fzero(Diag(1+iiBstR(jSym,1)),NNBSTR(jSym,1))

                  Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        Diag(jrs) = Diag(jrs) + Lrs(krs,jvc)**2

                     End Do

                  End Do

               EndIf

               CALL CWTIME(TCS2,TWS2)
               tscrn(1) = tscrn(1) + (TCS2 - TCS1)
               tscrn(2) = tscrn(2) + (TWS2 - TWS1)
!                                                                      *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
               Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,    &
     &                              nSym)
               mDen=1
               Call Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,mDen)

               CALL CWTIME(TCX1,TWX1)

! *** Reorder vectors to Full-dimensions
! ***
! *** Vectors are returned in the storage LaJ,b with the restriction:
! ***
! ***    Sym(a).ge.Sym(b)
! ***
! *** and blocked in shell pairs

               CALL CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,L_Full,     &
     &                            SvShp,nnShl,iShp_rs,nnShl_tot)

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)


               IF (DoScreen) THEN

                  CALL CWTIME(TCS1,TWS1)

! --- Compute DH(a,b)=sqrt(D(a,b)) from the updated diagonals.
! ---                              Only the symmetry blocks with
! ---                              compound symmetry JSYM are computed
! --------------------------------------------------------------------
                   ired1 = 1 ! location of the 1st red set
                   Call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,          &
     &                               DIAH,DIAG)

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               ENDIF


               Do jDen=1,nDen

                 Do kSym=1,nSym

                    lSym=MulD2h(JSYM,kSym)

                    nkOrb = nFIorb(kSym)*(2-jDen)                       &
     &                    + nChM(kSym)*(jDen-1)

                    If (nBas(lsym).eq.0) Cycle

                    Do jK=1,nkOrb

                       jK_a = jK + kOff(kSym,jDen)

                       Lab%A0(1:nBas(lSym)*JNUM)=Zero

                     IF (DoScreen) THEN

                        CALL CWTIME(TCS1,TWS1)
!------------------------------------------------------------------
! --- Setup the screening
!------------------------------------------------------------------

                        If (jDen.eq.2) Then
                           Do ik=1,nBas(kSym)
                              AbsC(ik)=abs(Ash(2)%SB(kSym)%A2(ik,jK))
                           End Do
                        Else
                           Do ik=1,nBas(kSym)
                              AbsC(ik) = abs(MSQ%SB(kSym)%A2(ik,jK))
                           End Do
                        Endif

                        If (lSym.ge.kSym) Then
! --------------------------------------------------------------
! --- Y(l)[k] = sum_n  DH(l,n) * |C(n)[k]|
!===============================================================
                           Mode(1:1)='N'
                           n1 = nBas(lSym)
                           n2 = nBas(kSym)

                        Else
! --------------------------------------------------------------
! --- Y(l)[k] = sum_n  DH(n,l) * |C(n)[k]|
!===============================================================
                           Mode(1:1)='T'
                           n1 = nBas(kSym)
                           n2 = nBas(lSym)

                        EndIf

                        If (n1>0)                                       &
     &                  CALL DGEMV_(Mode(1:1),n1,n2,                    &
     &                             ONE,DIAH%SB(lSym,kSym)%A2,n1,        &
     &                                 AbsC,1,                          &
     &                            ZERO,Ylk(1,jK_a),1)


! --- List the shells present in Y(l)[k] by the largest element
                        Do ish=1,nShell
                           YshMax=zero
                           Do ibs=1,nBasSh(lSym,ish)
                              YshMax = Max(YshMax,                      &
     &                          Ylk(koffSh(ish,lSym)+ibs,jK_a))
                           End Do
                           MLk(ish,jK_a) = YshMax
                        End Do

! --- Sort the lists ML[k]
                        Do ish=1,nShell
                           Indx(ish,jk_a) = ish
                        End Do

! **** Sort the list
                        numSh=0  ! # of significant shells
                        jml=1
                        Do while (jml.le.nShell)

                           YMax=MLk(jml,jK_a)
                           jmlmax=jml

                           Do iml=jml+1,nShell  ! get the max
                              If (MLk(iml,jK_a).gt.YMax) then
                                 YMax = MLk(iml,jK_a)
                                 jmlmax = iml
                              Endif
                           End Do

                           If(jmlmax.ne.jml) then  ! swap positions
                             xTmp = MLk(jml,jK_a)
                             iTmp = Indx(jml,jk_a)
                             MLk(jml,jK_a) = YMax
                             Indx(jml,jk_a) = Indx(jmlmax,jk_a)
                             MLk(jmlmax,jK_a) = xTmp
                             Indx(jmlmax,jk_a) = iTmp
                           Endif

! --- Exact bounds (quadratic scaling of the MO transformation)
!
!                           If(MLk(jml,jK_a)*MLk(1,jK_a)
!     &                                         .ge. tau(jDen))then
!
! --- Here we use a non-exact bound for the exchange matrix to achieve
! --- linear scaling. The positive definiteness of the exchange matrix
! --- combined with the structure of the density matrix makes this
! --- bound acceptable and likely to be almost exact for what concerns
! --- the exchange energy
                           If ( MLk(jml,jK_a) .ge. xtau(jDen) ) then
                             numSh = numSh + 1
                           else
                             jml=nShell  ! exit the loop
                           endif

                           jml=jml+1

                        End Do

                        Indx(0,jk_a) = numsh

                        CALL CWTIME(TCS2,TWS2)
                        tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                        tscrn(2) = tscrn(2) + (TWS2 - TWS1)
!------------------------------------------------------------------
                     ENDIF    ! Screening setup



! --- Transform vectors for shells in the list ML[k]
!
! --- Screening based on the Frobenius norm: sqrt(sum_ij  A(i,j)^2)
!
! ---   || La,J[k] ||  .le.  || Lab,J || * || Cb[k] ||

                     CALL CWTIME(TCT1,TWT1)

                     Do iSh=1,Indx(0,jk_a)

                        iaSh = Indx(iSh,jk_a)

                        Lab%Keep(iaSh,1)=.True.

                        ibcount=0

                        Do ibSh=1,nShell

                           iOffShb = kOffSh(ibSh,kSym)

                           iShp = iTri(iaSh,ibSh)

                           If ( iShp_rs(iShp)<=0) Cycle

                           If ( nnBstRSh(JSym,iShp_rs(iShp),iLoc)*      &
     &                       nBasSh(lSym,iaSh)*                         &
     &                       nBasSh(kSym,ibSh) .gt. 0                   &
     &                       .and. sqrt(abs(SumAClk(ibSh,jK_a)*         &
     &                                  SvShp(iShp_rs(iShp),1) ))       &
     &                       .ge. thrv(jDen) )Then

                             ibcount = ibcount + 1

                             IF (lSym.ge.kSym) Then

                                l1=1
                                If (iaSh<ibSh) l1=2

! ---  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
! ---------------------------------------
                                Mode(1:1)='N'
                                n1=nBasSh(lSym,iaSh)*JNUM
                                n2=nBasSh(kSym,ibSh)

                                If (JDen.eq.2) Then
                                  CALL DGEMV_(Mode(1:1),n1,n2,          &
     &                     One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1,&
     &                              Ash(2)%SB(kSym)%A2(1+ioffShb:,jK),1,&
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                Else
                                  CALL DGEMV_(Mode(1:1),n1,n2,          &
     &                     One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1,&
     &                              MSQ%SB(kSym)%A2(1+ioffShb:,jK),1,   &
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                End If

                             Else   ! lSym < kSym

                                l1=1
                                If (ibSh<iaSh) l1=2


! ---  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
! ---------------------------------------
                                Mode(1:1)='T'
                                n1=nBasSh(kSym,ibSh)
                                n2=JNUM*nBasSh(lSym,iaSh)

                                If (JDen.eq.2) Then
                                  CALL DGEMV_(Mode(1:1),n1,n2,          &
     &              One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1,       &
     &                              Ash(2)%SB(kSym)%A2(1+ioffShb:,jK),1,&
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                Else
                                  CALL DGEMV_(Mode(1:1),n1,n2,          &
     &              One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1,       &
     &                              MSQ%SB(kSym)%A2(1+ioffShb:,jK),1,   &
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                End If
                             EndIf


                           End If

                         End Do

! --- The following re-assignement is used later on to check if the
! --- iaSh vector LaJ[k] can be neglected because identically zero

                         If (ibcount==0) Lab%Keep(iaSh,1) = .False.

                      End Do

                     CALL CWTIME(TCT2,TWT2)
                     tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                     tmotr(2) = tmotr(2) + (TWT2 - TWT1)

! --- Prepare the J-screening

                      CALL CWTIME(TCS1,TWS1)

                      Do iSh=1,Indx(0,jk_a)

                         iaSh = Indx(iSh,jk_a)

                         If (.NOT.Lab%Keep(iaSh,1)) Cycle

                         IF (lSym.ge.kSym) Then

! ---  Faa,[k] = sum_J  (LaJ[k])**2
! ----------------------------------
                            Inc = nBasSh(lSym,iaSh)
                            n1 = 1

                         Else   ! lSym < kSym

! ---  Faa,[k] = sum_J  (LJa[k])**2
! ----------------------------------
                            Inc = 1
                            n1 = JNUM

                         End If

                         Tmp=Zero
                         Do ia=1,nBasSh(lSym,iaSh)
                            Fia(ia)=DDot_(JNUM,                         &
     &                        Lab%SB(iaSh,lSym,1)%A(1+n1*(ia-1):),Inc,  &
     &                        Lab%SB(iaSh,lSym,1)%A(1+n1*(ia-1):),Inc)
                            Tmp=Max(Abs(Fia(ia)),Tmp)
                         End Do

                         Faa(iaSh)=Tmp

                      End Do

                      CALL CWTIME(TCS2,TWS2)
                      tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                      tscrn(2) = tscrn(2) + (TWS2 - TWS1)


!------------------------------------------------------------
! --- Compute exchange matrix for the interacting shell pairs
!------------------------------------------------------------
                     CALL CWTIME(TCX1,TWX1)


                     Do lSh=1,Indx(0,jk_a)

                        iaSh = Indx(lSh,jk_a)

                        iaSkip=Merge(1,0,Lab%Keep(iaSh,1))

                        mSh = 1

                        Do while (mSh.le.Indx(0,jk_a))

                           ibSh = Indx(mSh,jk_a)

                           ibSkip=Merge(1,0,Lab%Keep(ibSh,1))

                           iShp = iTri(iaSh,ibSh)

                           iOffShb = kOffSh(ibSh,lSym)

                           iOffAB = nnBfShp(iShp,lSym)

                           xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                           If (MLk(mSh,jK_a)*MLk(lSh,jK_a)              &
     &                         .lt. tau(jDen))Then

                                mSh = Indx(0,jk_a) !skip rest

                           ElseIf(iaSh.eq.ibSh                          &
     &                            .and.xFab.ge.tau(jDen)/MaxRedT        &
     &                            .and.iaSkip.eq.1)Then

                               IF (lSym.ge.kSym) Then

! ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(a,J)[k] * L(b,J)[k]
! -------------------------------------------------------------------
                               nBs = nBasSh(lSym,iaSh)

                               CALL DGEMM_Tri('N','T',                  &
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,&
     &                FActX(jDen),Lab%SB(iaSh,lSym,1)%A,nBs,            &
     &                            Lab%SB(iaSh,lSym,1)%A,nBs,            &
     &                        ONE,KLT(jDen)%SB(lSym)%A1(iOffAB+1:),nBs)

                               ELSE   ! lSym < kSym

! ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(J,a)[k] * L(J,b)[k]
! -------------------------------------------------------------------

                               nBs = nBasSh(lSym,iaSh)

                               CALL DGEMM_Tri('T','N',                  &
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,&
     &                FActX(jDen),Lab%SB(iaSh,lSym,1)%A,JNUM,           &
     &                            Lab%SB(iaSh,lSym,1)%A,JNUM,           &
     &                        ONE,KLT(jDen)%SB(lSym)%A1(iOffAB+1:),nBs)

                               End If


                           ElseIf (iaSh.gt.ibSh                         &
     &                             .and.xFab.ge.tau(jDen)/MaxRedT       &
     &                             .and. iaSkip*ibSkip.eq.1) Then

                                 IF (lSym.ge.kSym) Then
! ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(a,J)[k] * L(b,J)[k]
! -------------------------------------------------------------------
                                  nBsa = nBasSh(lSym,iaSh)
                                  nBsb = nBasSh(lSym,ibSh)

                                  CALL DGEMM_('N','T',                  &
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,&
     &                FActX(jDen),Lab%SB(iaSh,lSym,1)%A,nBsa,           &
     &                            Lab%SB(ibSh,lSym,1)%A,nBsb,           &
     &                        ONE,KLT(jDen)%SB(lSym)%A1(iOffAB+1:),nBsa)

                               ELSE   ! lSym < kSym

! ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(J,a)[k] * L(J,b)[k]
! -------------------------------------------------------------------

                                  nBs = nBasSh(lSym,iaSh)

                                  CALL DGEMM_('T','N',                  &
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,&
     &                FActX(jDen),Lab%SB(iaSh,lSym,1)%A,JNUM,           &
     &                            Lab%SB(ibSh,lSym,1)%A,JNUM,           &
     &                        ONE,KLT(jDen)%SB(lSym)%A1(iOffAB+1:),nBs)

                               End If


                           EndIf


                           mSh = mSh + 1  ! update shell counter


                        End Do

                     End Do

                     CALL CWTIME(TCX2,TWX2)
                     texch(1) = texch(1) + (TCX2 - TCX1)
                     texch(2) = texch(2) + (TWX2 - TWX1)


                    End Do  ! loop over k MOs

                 End Do   ! loop over MOs symmetry

               End Do   ! loop over densities

               Call Deallocate_Lab(Lab)
               Call Deallocate_L_Full(L_Full)
!                                                                      *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
               DoScreen=.false. ! avoid redo screening inside batch loop


! ************  END EXCHANGE CONTRIBUTION  ****************

! --- Diagonals updating. It only makes sense if Nscreen > 0

               If (Update .and. Nscreen .gt. 0) Then

                  CALL CWTIME(TCS1,TWS1)
! ---------------------------------------------------------------------
! --- update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2
!
! --- subtraction is done in the 1st reduced set
#if defined (_MOLCAS_MPP_)
                  If (nProcs .gt. 1 .and. Is_Real_Par()) then

                   Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) - iiBstR(JSYM,1)

                     Do jvc=1,JNUM

                        DiagJ(jrs) = DiagJ(jrs) + Lrs(krs,jvc)**2
                     End Do

                   End Do

                  Else

                   Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        Diag(jrs) = Diag(jrs) - Lrs(krs,jvc)**2
                     End Do

                   End Do

                  EndIf

#else
                  Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        Diag(jrs) = Diag(jrs) - Lrs(krs,jvc)**2
                     End Do

                  End Do
#endif

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               EndIf

               ! Lvw,J , strictly LT storage
               iSwap = 5
               Call Allocate_SBA(Lxy,nAorb,nAorb,nVec,JSYM,nSym,iSwap)
               iSwap = 0  ! Lvb,J are returned
               Call Allocate_SBA(Laq(1),nAorb,nBas,nVec,JSYM,nSym,iSwap)
! --------------------------------------------------------------------
! --- First half Active transformation  Lvb,J = sum_a  C(v,a) * Lab,J
! --------------------------------------------------------------------

               CALL CWTIME(TCINT1,TWINT1)

! --- Set up the skipping flags
! --- The memory used before for the full-dimension AO-vectors
! ---     is now re-used to store half and full transformed
! ---     vectors in the active space
! -------------------------------------------------------------

               kMOs = 1  !
               nMOs = 1  ! Active MOs (1st set)

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,              &
     &                          JSYM,iSwap,IREDC,nMOs,kMOs,Ash,         &
     &                          Laq(1),DoRead)

               if (irc.ne.0) then
                  RETURN
               endif


! --------------------------------------------------------------------
! --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
! --------------------------------------------------------------------
               If (JSYM.eq.1) Then

                  Do iSymb=1,nSym

                     NAv = nAorb(iSymb)

                     If(NAv.gt.0)Then

                      Do JVC=1,JNUM


                       CALL DGEMM_Tri('N','T',NAv,NAv,NBAS(iSymb),      &
     &                            One,Laq(1)%SB(iSymb)%A3(:,:,JVC),NAv, &
     &                                Ash(1)%SB(iSymb)%A2,NAv,          &
     &                           Zero,Lxy%SB(iSymb)%A2(:,JVC),NAv)

                      End Do

                     EndIf

                  End Do

               Else

! --------------------------------------------------------------------
! --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
! --------------------------------------------------------------------
                  Do iSymb=1,nSym

                     iSymv = MulD2h(JSYM,iSymb)
                     NAv = nAorb(iSymv)
                     NAw = nAorb(iSymb) ! iSymb=iSymw

                     If(NAv*NAw.ne.0 .and. iSymv.lt.iSymb)Then

                      Do JVC=1,JNUM


                       CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),         &
     &                            One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv, &
     &                                Ash(1)%SB(iSymb)%A2,NAw,          &
     &                           Zero,Lxy%SB(iSymv)%A2(:,JVC),NAv)

                      End Do

                     EndIf

                  End Do

               EndIf
!
!
! *************** EVALUATION OF THE (WA|XY) INTEGRALS ***********

               DoTraInt = JRED.eq.myJRED2.and.iBatch.eq.nBatch

               CALL CHO_eval_waxy(irc,Scr,Laq(1),Lxy,W_PWXY,nAorb,      &
     &                            JSYM,JNUM,DoTraInt,CMO)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  RETURN
               endif

               Call Deallocate_SBA(Lxy)
               Call Deallocate_SBA(Laq(1))

! --------------------------------------------------------------------
! --------------------------------------------------------------------

            END DO  ! end batch loop


            If(JSYM.eq.1)Then
! --- backtransform fock matrix to full storage
               mode = 'tofull'
               add  = .true.
               Call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,                &
     &                           FLT,Frs,mode,add)
            EndIf

! --- free memory
            Call mma_deallocate(Lrs)

            If(JSYM.eq.1)Then
              Call mma_deallocate(Frs)
              Call mma_deallocate(Drs)
            EndIf

999         Continue

! --- Screening control section
            DoScreen = kscreen.eq.Nscreen

            if (.not.DoScreen) then
                kscreen = kscreen + 1
            else
                kscreen = 1
            endif

#if defined (_MOLCAS_MPP_)
            If (nProcs.gt.1 .and. Update .and. DoScreen                 &
     &          .and. Is_Real_Par()) Then
               Call GaDsum(DiagJ,nnBSTR(JSYM,1))
               Call Daxpy_(nnBSTR(JSYM,1),xone,DiagJ,1,                 &
     &                    Diag(1+iiBstR(JSYM,1)),1)
               Call Fzero(DiagJ,nnBSTR(JSYM,1))
            EndIf
!--- Need to activate the screening to setup the contributing shell
!--- indeces the first time the loop is entered .OR. whenever other nodes
!--- have performed screening in the meanwhile
            If (nProcs.gt.1 .and. .not.DoScreen .and. nVrs.eq.0         &
     &          .and. Is_Real_Par()) Then
               ntv0=ntv0+1
               DoScreen = (JRED.lt.myJRED1 .or. ntv0.ge.Nscreen)
               if (DoScreen) ntv0=0
            EndIf
#endif

         END DO   ! loop over red sets

         Call Deallocate_twxy(Scr)

      END DO  ! loop over JSYM

! --- Accumulate Coulomb and Exchange contributions
      Do jDen=1,nDen

        Do iSym=1,nSym

         Do iaSh=1,nShell

            ioffa = kOffSh(iaSh,iSym)

! --- ibSh < iaSh
! ---------------
            Do ibSh=1,iaSh-1

               iShp = iaSh*(iaSh-1)/2 + ibSh

               iOffAB = nnBfShp(iShp,iSym)

               ioffb = kOffSh(ibSh,iSym)

               Do ib=1,nBasSh(iSym,ibSh)

                Do ia=1,nBasSh(iSym,iaSh)

                  iab = nBasSh(iSym,iaSh)*(ib-1) + ia

                  iag = ioffa + ia
                  ibg = ioffb + ib

                  FLT(jDen)%SB(iSym)%A1(iTri(iag,ibg)) =                &
     &                FLT(jDen)%SB(iSym)%A1(iTri(iag,ibg))              &
     &              + KLT(jDen)%SB(iSym)%A1(iOffAB+iab)

                End Do

               End Do

            End Do

! --- ibSh = iaSh
! ---------------
            iShp = iaSh*(iaSh+1)/2

            iOffAB = nnBfShp(iShp,iSym)

            Do ib=1,nBasSh(iSym,iaSh)

             Do ia=ib,nBasSh(iSym,iaSh)

               iab = ia*(ia-1)/2 + ib

               iag = ioffa + ia
               ibg = ioffa + ib

               FLT(jDen)%SB(iSym)%A1(iTri(iag,ibg)) =                   &
     &             FLT(jDen)%SB(iSym)%A1(iTri(iag,ibg))                 &
     &           + KLT(jDen)%SB(iSym)%A1(iOffAB+iab)

             End Do

            End Do


         End Do

        End Do

      End Do


      Call mma_deallocate(Faa)
      Call mma_deallocate(Fia)
      Call mma_deallocate(SvShp)
      Call mma_deallocate(iShp_rs)
      Call mma_deallocate(nnBfShp)
      Call mma_deallocate(kOffSh)
      Call mma_deallocate(Indx)
      Call mma_deallocate(SumAClk)
      Call mma_deallocate(MLk)
      Call mma_deallocate(Ylk)
      Call mma_deallocate(AbsC)
      Call Deallocate_NDSBA(DIAH)
      Do jDen = 1, nDen
         Call Deallocate_DSBA(KLT(jDen))
      End Do
#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par())                 &
     &    Call mma_deallocate(DiagJ)
#endif
      Call mma_deallocate(Diag)

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

!
!---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky RASSCF timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '  &
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '  &
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'SCREENING OVERHEAD               '  &
     &                           //'         ',tscrn(1),tscrn(2)
         Write(6,'(2x,A26,2f10.2)')'MO HALF-TRANSFORM VECTORS        '  &
     &                           //'         ',tmotr(1),tmotr(2)
         Write(6,'(2x,A26,2f10.2)')'EXCHANGE                         '  &
     &                           //'         ',texch(1),texch(2)
         Write(6,'(2x,A26,2f10.2)')'(PU|VX) INTEGRALS                '  &
     &                           //'         ',tintg(1),tintg(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '  &
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in CASSCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      WRITE(6,'(6X,A)')'***** INACTIVE FOCK MATRIX ***** '
      DO ISYM=1,NSYM
        IF( NBAS(ISYM).GT.0 ) THEN
          WRITE(6,'(6X,A)')
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call TRIPRT('','',FLT(1)%SB(ISYM)%A1,NBAS(ISYM))
        ENDIF
      END DO
      IF(DoActive)THEN
        WRITE(6,'(6X,A)')
        WRITE(6,'(6X,A)')'***** ACTIVE FOCK MATRIX ***** '
        DO ISYM=1,NSYM
         IF( NBAS(ISYM).GT.0 ) THEN
           WRITE(6,'(6X,A)')
           WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
           call TRIPRT('','',FLT(2)%SB(ISYM)%A2,NBAS(ISYM))
         ENDIF
        END DO
      END IF

      endif

#endif

      Return
      END
