************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE CHO_LK_CASSCF(ipDI,ipDA1,ipFI,ipFA,ipKLT,ipMSQ,ipInt,
     &             FactXI,nFIorb,nAorb,nChM,Ash,DoActive,
     &             nScreen,dmpk,dFmat,ExFac)

**********************************************************************
*  Author : F. Aquilante
*
*           This routine makes use of the "Local K" scheme for
*                computing the exchange terms in both the Inactive
*                and Active Fock matrix
*
C *************** INACTIVE AO-BASIS FOCK MATRIX **********************
C
C   FI(ab) = 2 * sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
C
C ***************   ACTIVE AO-BASIS FOCK MATRIX **********************
C
C   FA(ab) = sum_J  Lab,J * U(J)  -  0.5 * sum_Jw  Lwa,J * Lwb,J
C
C ***************   (WA|XY) integrals     ****************************
C
C   (WA|XY) = sum_J  L(wa,J) * L(xy,J)
C
**********************************************************************
C
C      V(J) = sum_gd  Lgd,J * DI(gd)
C      U(J) = sum_gd  Lgd,J * DA(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index   belonging to (Frozen+Inactive)
C      u,w,x,y:  MO-indeces belonging to (Active)
C
**********************************************************************

#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      use ChoArr, only: nBasSh, nDimRS
      use ChoSwp, only: nnBstRSh, InfVec, IndRed
      use Data_Structures, only: DSBA_Type, SBA_Type
      use Data_Structures, only: Allocate_SBA, Deallocate_SBA
      use Data_Structures, only: twxy_Type
      use Data_Structures, only: Allocate_twxy, Deallocate_twxy
      use Data_Structures, only: NDSBA_Type, Allocate_NDSBA,
     &                           Deallocate_NDSBA
      use Data_Structures, only: Allocate_L_Full, Deallocate_L_Full
      use Data_Structures, only: L_Full_Type
      use Data_Structures, only: Allocate_Lab, Deallocate_Lab,
     &                           Lab_Type
      Implicit Real*8 (a-h,o-z)

      Integer   ipDLT(2),ipFLT(2),ipKLT(2)
      Integer   kOff(8,2),nnA(8,8)
      Integer   ISTLT(8),ISTSQ(8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2)
      Real*8    tmotr(2),tscrn(2)
      Integer   ipDD(2)

      Type (DSBA_Type)   Ash(2)
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
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Real*8, parameter:: xone = -one

      Logical add
      Character(LEN=6) mode
      Real*8   LKThr
      Real*8, External :: Cho_LK_ScreeningThreshold
      Integer, External :: Cho_LK_MaxVecPerBatch

      Real*8, Allocatable:: Lrs(:,:), Drs(:,:), Frs(:,:)
      Real*8, Allocatable:: VJ(:)

      Integer, Allocatable:: nnBfShp(:,:), kOffSh(:,:),
     &                       iShp_rs(:), Indx(:,:)
      Real*8, Allocatable:: SvShp(:,:), Diag(:), AbsC(:), SumAClk(:,:),
     &                      Ylk(:,:), MLk(:,:), Faa(:), Fia(:)
#if defined (_MOLCAS_MPP_)
      Real*8, Allocatable:: DiagJ(:)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Interface

        Subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,
     &                         iSwap,IREDC,nDen,kDen,MOs,ChoT,
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
*                                                                      *
************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in CASSCF-debug
#endif

      ipMO = 0 ! Dummy initiate
      DoTraInt = .false.
      IREDC = -1  ! unknown reduced set in core

      ipDLT(1) = ipDI    ! some definitions
      ipDLT(2) = ipDA1
      ipFLT(1) = ipFI
      ipFLT(2) = ipFA

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

C ==================================================================
        Call set_nnA(nSym,nAorb,nnA)

c --- Various offsets
c --------------------
        MaxB=nBas(1)
        ISTLT(1)=0
        ISTSQ(1)=0
      DO ISYM=2,NSYM
        MaxB=Max(MaxB,nBas(iSym))
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        NBQ=NBAS(ISYM-1)**2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive D and F matrices
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NBQ ! Diagonal integrals in full
      END DO

**************************************************

      nnO=0
      nIt=0
      DO jDen=1,nDen

         kOff(1,jDen) = nnO

         DO ISYM=2,NSYM

            nnO = nnO + nFIorb(iSym-1)*(2-jDen)
     &                + nChM(iSym-1)*(jDen-1)

            kOff(iSym,jDen) = nnO

         END DO

         nnO = nnO + nFIorb(nSym)*(2-jDen)
     &             + nChM(nSym)*(jDen-1)

         nIt = nIt + nnO*(2-jDen) ! tot # of (inactive+frozen) orbitals

      END DO

      nAt = nnO - nIt  ! tot # of active orb (decomp. Active density)

C --- Define max number of vectors to be treated in core at once

      MaxVecPerBatch=Cho_LK_MaxVecPerBatch()

C --- Define the basic screening threshold

      LKThr=Cho_LK_ScreeningThreshold(dFmat)

C --- Adjust the damping according to the Abs(Max BLB)
Ctbp, may 2013: adjustment moved to Cho_LK_ScreeningThreshold
      fcorr = dmpK
Ctbp  If (dFmat.gt.zero) Then
Ctbp     If (dFmat.lt.1.0d3*LKThr) Then
Ctbp        fcorr=dmpk*1.0d-2
Ctbp     EndIf
Ctbp     If (dFmat.le.LKThr) fcorr=fcorr*1.0d-2
Ctbp  EndIf

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

C --- Vector MO transformation screening thresholds
      NumVT=NumChT
      Call GAIGOP_SCAL(NumVT,'+')
      thrv(1) = ( sqrt(LKThr/(Max(1,nIt)*NumVT)) )*fcorr
      thrv(2) = ( sqrt(LKThr/(Max(1,nAt)*NumVT)) )*fcorr

      CALL mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')
      ipDD(1)= ip_of_Work(DIAG(1))

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

C *************** Read the diagonal integrals (stored as 1st red set)
      If (Update) CALL CHO_IODIAG(DIAG,2) ! 2 means "read"

c --- allocate memory for sqrt(D(a,b)) stored in full (squared) dim
      Call Allocate_NDSBA(DIAH,nBas,nBas,nSym)
      DIAH%A0(:)=Zero

c --- allocate memory for the abs(C(l)[k])
      Call mma_allocate(AbsC,MaxB,Label='AbsC')
      ipAbs = ip_of_Work(AbsC(1))

c --- allocate memory for the Y(l)[k] vectors
      Call mma_allocate(Ylk,MaxB,nnO,Label='Ylk')

c --- allocate memory for the ML[k] lists of largest elements
c --- in significant shells
      Call mma_allocate(MLk,nShell,nnO,Label='MLk')

c --- allocate memory for the list of  S:= sum_l abs(C(l)[k])
c --- for each shell
      Call mma_allocate(SumAClk,nShell,nnO,Label='SumAClk')

c --- allocate memory for the Index arrays
      Call mma_allocate(Indx,[0,nShell],[1,nnO],Label='Indx')

c --- allocate memory for kOffSh
      Call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')

c --- allocate memory for nnBfShp
      Call mma_allocate(nnBfShp,nnShl_tot,nSym,Label='nnBfShp')

c --- allocate memory for iShp_rs
      Call mma_allocate(iShp_rs,nnShl_tot,Label='iShp_rs')

c --- allocate memory for the shell-pair Frobenius norm of the vectors
      Call mma_allocate(SvShp,nnShl,2,Label='SvShp')


C *** Compute Shell Offsets ( MOs and transformed vectors)

      MxBasSh = 0

      Do iSyma=1,nSym

         LKsh=0

         Do iaSh=1,nShell    ! kOffSh(iSh,iSym)

            kOffSh(iaSh,iSyma) = LKsh

            LKsh = LKsh + nBasSh(iSyma,iaSh)

            MxBasSh = Max(MxBasSh,nBasSh(iSyma,iaSh))

         End Do

      End Do

c --- allocate memory for the Diagonal of the Fock matrix
      Call mma_allocate(Fia,MxBasSh,Label='Fia')
      Call mma_allocate(Faa,nShell,Label='Faa')
      Faa(:)=Zero
      Fia(:)=Zero

C *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
      Do jDen=1,nDen
         Do kSym=1,nSym

            If (jDen.eq.2) Then
               Do jK=1,nChM(kSym)
                  jK_a = jK + kOff(kSym,jDen)

               Do iaSh=1,nShell

                  ipMsh =  kOffSh(iaSh,kSym)

                  SKsh=zero
                  Do ik=1,nBasSh(kSym,iaSh)
                     SKsh = SKsh + Ash(2)%SB(kSym)%A2(ipMsh+ik,jK)**2
                  End Do

                  SumAClk(iaSh,jk_a) = SKsh

               End Do
               End Do

            Else
               Do jK=1,nFIorb(kSym)
                  jK_a = jK + kOff(kSym,jDen)

               ipMO = ipMSQ + ISTSQ(kSym) + nBas(kSym)*(jK-1)

               Do iaSh=1,nShell

                  ipMsh = ipMO + kOffSh(iaSh,kSym)

                  SKsh=zero
                  Do ik=0,nBasSh(kSym,iaSh)-1
                     SKsh = SKsh + Work(ipMsh+ik)**2
                  End Do

                  SumAClk(iaSh,jk_a) = SKsh

               End Do
               End Do

            End If

         End Do
      End Do

C *** Compute Shell-pair Offsets in the K-matrix

         Do iSyma=1,nSym

            LKshp=0

            Do iaSh=1,nShell

             Do ibSh=1,iaSh

               iShp = iaSh*(iaSh-1)/2 + ibSh

               nnBfShp(iShp,iSyma) = LKShp

               LKShp = LKShp + nBasSh(iSyma,iaSh)*nBasSh(iSyma,ibSh)
     &                       - (1-Min((iaSh-ibSh),1))*nBasSh(iSyma,iaSh)
     &                       * (nBasSh(iSyma,iaSh) - 1)/2

             End Do

            End Do

         End Do

C *** Mapping shell pairs from the full to the reduced set

      Call Mk_iShp_rs(iShp_rs,nShell)

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

        NumCV=NumCho(jSym)
        Call GAIGOP_SCAL(NumCV,'max')
        If (NumCV .lt. 1) Cycle

        JNUM=1
        Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,
     &                       Memory=LFULL)

        iCase = 1  ! (wa|xy)
        Call Allocate_twxy(Scr,nAorb,nBas,JSYM,nSym,iCase)

C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
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

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         iLoc = 3 ! use scratch location in reduced index arrays

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
#if defined (_MOLCAS_MPP_)
         myJRED1=JRED1 ! first red set present on this node
         ntv0=0
#endif
         myJRED2=JRED2 ! last  red set present on this node

c --- entire red sets range for parallel run
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
c           !set index arrays at iLoc
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
C --- Transform the densities to reduced set storage
               mode = 'toreds'
               add  = .false.
               Call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,
     &                           [ipDLT],Drs,mode,add)
            EndIf

C --- BATCH over the vectors ----------------------------

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

               CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,
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
C ************ COULOMB CONTRIBUTIONS  *********************
C
C  jDen=1   ---> inactive fock matrix
C  jDen=2   ---> active fock matrix
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} =  sum_rs  L(rs,{#J}) * DI(rs)
C==========================================================
C
                   CALL DGEMV_('T',nRS,JNUM,
     &                  ONE,Lrs,nRS,
     &                  Drs(:,jDen),1,ZERO,VJ,1)

C --- F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
C===============================================================

                   CALL DGEMV_('N',nRS,JNUM,
     &                 FactC(jDen),Lrs,nRS,
     &                 VJ,1,Fact,Frs(:,jDen),1)

                 End Do
                 Call mma_deallocate(VJ)

                 CALL CWTIME(TCC2,TWC2)
                 tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                 tcoul(2) = tcoul(2) + (TWC2 - TWC1)


               EndIf  ! Coulomb contribution


C *************** EXCHANGE CONTRIBUTIONS  ***********************

               CALL CWTIME(TCS1,TWS1)
C ---------------------------------------------------------------------
C --- Estimate the diagonals :   D(a,b) = sum_J (Lab,J)^2
C
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
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
               Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,
     &                              nSym)
               mDen=1
               Call Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,mDen)

               CALL CWTIME(TCX1,TWX1)

C *** Reorder vectors to Full-dimensions
C ***
C *** Vectors are returned in the storage LaJ,b with the restriction:
C ***
C ***    Sym(a).ge.Sym(b)
C ***
C *** and blocked in shell pairs

               CALL CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,L_Full,
     &                            SvShp,nnShl,iShp_rs,nnShl_tot)

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)


               IF (DoScreen) THEN

                  CALL CWTIME(TCS1,TWS1)

c --- Compute DH(a,b)=sqrt(D(a,b)) from the updated diagonals.
c ---                              Only the symmetry blocks with
c ---                              compound symmetry JSYM are computed
c --------------------------------------------------------------------
                   ired1 = 1 ! location of the 1st red set
                   Call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,
     &                               DIAH,Work(ipDD(1)))

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               ENDIF


               Do jDen=1,nDen

                 Do kSym=1,nSym

                    lSym=MulD2h(JSYM,kSym)

                    nkOrb = nFIorb(kSym)*(2-jDen)
     &                    + nChM(kSym)*(jDen-1)

                    If (nBas(lsym).eq.0) Cycle

                    Do jK=1,nkOrb

                       jK_a = jK + kOff(kSym,jDen)

                       Lab%A0(1:nBas(lSym)*JNUM)=Zero

                     If (jDen.eq.1) Then

                     ipMO = ipMSQ + ISTSQ(kSym) + nBas(kSym)*(jK-1)

                     End If

                     IF (DoScreen) THEN

                        CALL CWTIME(TCS1,TWS1)
C------------------------------------------------------------------
C --- Setup the screening
C------------------------------------------------------------------

                        If (jDen.eq.2) Then
                        Do ik=1,nBas(kSym)
                           AbsC(ik)=abs(Ash(2)%SB(kSym)%A2(ik,jK))
                        End Do
                        Else
                        Do ik=0,nBas(kSym)-1
                           AbsC(1+ik) = abs(Work(ipMO+ik))
                        End Do
                        Endif

                        If (lSym.ge.kSym) Then
c --------------------------------------------------------------
C --- Y(l)[k] = sum_n  DH(l,n) * |C(n)[k]|
C===============================================================
                           Mode(1:1)='N'
                           n1 = nBas(lSym)
                           n2 = nBas(kSym)

                        Else
c --------------------------------------------------------------
C --- Y(l)[k] = sum_n  DH(n,l) * |C(n)[k]|
C===============================================================
                           Mode(1:1)='T'
                           n1 = nBas(kSym)
                           n2 = nBas(lSym)

                        EndIf

                        If (n1>0)
     &                  CALL DGEMV_(Mode(1:1),n1,n2,
     &                             ONE,DIAH%SB(lSym,kSym)%A2,n1,
     &                                 AbsC,1,
     &                            ZERO,Ylk(1,jK_a),1)


C --- List the shells present in Y(l)[k] by the largest element
                        Do ish=1,nShell
                           YshMax=zero
                           Do ibs=1,nBasSh(lSym,ish)
                              YshMax = Max(YshMax,
     &                          Ylk(koffSh(ish,lSym)+ibs,jK_a))
                           End Do
                           MLk(ish,jK_a) = YshMax
                        End Do

C --- Sort the lists ML[k]
                        Do ish=1,nShell
                           Indx(ish,jk_a) = ish
                        End Do

C **** Sort the list
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

c --- Exact bounds (quadratic scaling of the MO transformation)
c
c                           If(MLk(jml,jK_a)*MLk(1,jK_a)
c     &                                         .ge. tau(jDen))then
c
c --- Here we use a non-exact bound for the exchange matrix to achieve
c --- linear scaling. The positive definiteness of the exchange matrix
c --- combined with the structure of the density matrix makes this
c --- bound acceptable and likely to be almost exact for what concerns
c --- the exchange energy
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
C------------------------------------------------------------------
                     ENDIF    ! Screening setup



C --- Transform vectors for shells in the list ML[k]
C
C --- Screening based on the Frobenius norm: sqrt(sum_ij  A(i,j)^2)
C
C ---   || La,J[k] ||  .le.  || Lab,J || * || Cb[k] ||

                     CALL CWTIME(TCT1,TWT1)

                     Do iSh=1,Indx(0,jk_a)

                        iaSh = Indx(iSh,jk_a)

                        iOffSha = kOffSh(iaSh,lSym)

                        Lab%Keep(iaSh,1)=.True.

                        ibcount=0

                        Do ibSh=1,nShell

                           iOffShb = kOffSh(ibSh,kSym)

                           iShp = iTri(iaSh,ibSh)

                           If ( iShp_rs(iShp)<=0) Cycle

                           If ( nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                       nBasSh(lSym,iaSh)*
     &                       nBasSh(kSym,ibSh) .gt. 0
     &                       .and. sqrt(abs(SumAClk(ibSh,jK_a)*
     &                                  SvShp(iShp_rs(iShp),1) ))
     &                       .ge. thrv(jDen) )Then

                             ibcount = ibcount + 1

                             IF (lSym.ge.kSym) Then

                                l1=1
                                If (iaSh<ibSh) l1=2

C ---  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
C ---------------------------------------
                                Mode(1:1)='N'
                                n1=nBasSh(lSym,iaSh)*JNUM
                                n2=nBasSh(kSym,ibSh)

                                If (JDen.eq.2) Then
                                  CALL DGEMV_(Mode(1:1),n1,n2,
     &                     One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1,
     &                              Ash(2)%SB(kSym)%A2(1+ioffShb:,jK),1,
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                Else
                                  CALL DGEMV_(Mode(1:1),n1,n2,
     &                     One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1,
     &                                 Work(ipMO+ioffShb),1,
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                End If

                             Else   ! lSym < kSym

                                l1=1
                                If (ibSh<iaSh) l1=2


C ---  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
C ---------------------------------------
                                Mode(1:1)='T'
                                n1=nBasSh(kSym,ibSh)
                                n2=JNUM*nBasSh(lSym,iaSh)

                                If (JDen.eq.2) Then
                                  CALL DGEMV_(Mode(1:1),n1,n2,
     &              One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1,
     &                              Ash(2)%SB(kSym)%A2(1+ioffShb:,jK),1,
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                Else
                                  CALL DGEMV_(Mode(1:1),n1,n2,
     &              One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1,
     &                                 Work(ipMO+ioffShb),1,
     &                     One,Lab%SB(iaSh,lSym,1)%A,1)
                                End If
                             EndIf


                           End If

                         End Do

c --- The following re-assignement is used later on to check if the
c --- iaSh vector LaJ[k] can be neglected because identically zero

                         If (ibcount==0) Lab%Keep(iaSh,1) = .False.

                      End Do

                     CALL CWTIME(TCT2,TWT2)
                     tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                     tmotr(2) = tmotr(2) + (TWT2 - TWT1)

C --- Prepare the J-screening

                      CALL CWTIME(TCS1,TWS1)

                      Do iSh=1,Indx(0,jk_a)

                         iaSh = Indx(iSh,jk_a)

                         If (.NOT.Lab%Keep(iaSh,1)) Cycle

                         IF (lSym.ge.kSym) Then

C ---  Faa,[k] = sum_J  (LaJ[k])**2
C ----------------------------------
                            Inc = nBasSh(lSym,iaSh)
                            n1 = 1

                         Else   ! lSym < kSym

C ---  Faa,[k] = sum_J  (LJa[k])**2
C ----------------------------------
                            Inc = 1
                            n1 = JNUM

                         End If

                         Tmp=Zero
                         Do ia=1,nBasSh(lSym,iaSh)
                            Fia(ia)=DDot_(JNUM,
     &                        Lab%SB(iaSh,lSym,1)%A(1+n1*(ia-1)),Inc,
     &                        Lab%SB(iaSh,lSym,1)%A(1+n1*(ia-1)),Inc)
                            Tmp=Max(Abs(Fia(ia)),Tmp)
                         End Do

                         Faa(iaSh)=Tmp

                      End Do

                      CALL CWTIME(TCS2,TWS2)
                      tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                      tscrn(2) = tscrn(2) + (TWS2 - TWS1)


C------------------------------------------------------------
C --- Compute exchange matrix for the interacting shell pairs
C------------------------------------------------------------
                     CALL CWTIME(TCX1,TWX1)


                     Do lSh=1,Indx(0,jk_a)

                        iaSh = Indx(lSh,jk_a)

                        iaSkip=Merge(1,0,Lab%Keep(iaSh,1))

                        iOffSha = kOffSh(iaSh,lSym)

                        mSh = 1

                        Do while (mSh.le.Indx(0,jk_a))

                           ibSh = Indx(mSh,jk_a)

                           ibSkip=Merge(1,0,Lab%Keep(ibSh,1))

                           iShp = iTri(iaSh,ibSh)

                           iOffShb = kOffSh(ibSh,lSym)

                           iOffAB = nnBfShp(iShp,lSym)

                           ipKI = ipKLT(jDen) + ISTLT(lSym) + iOffAB

                           xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                           If (MLk(mSh,jK_a)*MLk(lSh,jK_a)
     &                         .lt. tau(jDen))Then

                                mSh = Indx(0,jk_a) !skip rest

                           ElseIf(iaSh.eq.ibSh
     &                            .and.xFab.ge.tau(jDen)/MaxRedT
     &                            .and.iaSkip.eq.1)Then

                               IF (lSym.ge.kSym) Then

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(a,J)[k] * L(b,J)[k]
C -------------------------------------------------------------------
                               nBs = nBasSh(lSym,iaSh)

                               CALL DGEMM_Tri('N','T',
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                          FActX(jDen),Lab%SB(iaSh,lSym,1)%A,nBs,
     &                                      Lab%SB(iaSh,lSym,1)%A,nBs,
     &                                        ONE,Work(ipKI),nBs)

                               ELSE   ! lSym < kSym

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(J,a)[k] * L(J,b)[k]
C -------------------------------------------------------------------

                               nBs = nBasSh(lSym,iaSh)

                               CALL DGEMM_Tri('T','N',
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                           FActX(jDen),Lab%SB(iaSh,lSym,1)%A,JNUM,
     &                                       Lab%SB(iaSh,lSym,1)%A,JNUM,
     &                                       ONE,Work(ipKI),nBs)

                               End If


                           ElseIf (iaSh.gt.ibSh
     &                             .and.xFab.ge.tau(jDen)/MaxRedT
     &                             .and. iaSkip*ibSkip.eq.1) Then

                                 IF (lSym.ge.kSym) Then
C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(a,J)[k] * L(b,J)[k]
C -------------------------------------------------------------------
                                  nBsa = nBasSh(lSym,iaSh)
                                  nBsb = nBasSh(lSym,ibSh)

                                  CALL DGEMM_('N','T',
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                           FActX(jDen),Lab%SB(iaSh,lSym,1)%A,nBsa,
     &                                       Lab%SB(ibSh,lSym,1)%A,nBsb,
     &                                       ONE,Work(ipKI),nBsa)

                               ELSE   ! lSym < kSym

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(J,a)[k] * L(J,b)[k]
C -------------------------------------------------------------------

                                  nBs = nBasSh(lSym,iaSh)

                                  CALL DGEMM_('T','N',
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                           FActX(jDen),Lab%SB(iaSh,lSym,1)%A,JNUM,
     &                                       Lab%SB(ibSh,lSym,1)%A,JNUM,
     &                                       ONE,Work(ipKI),nBs)

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
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
               DoScreen=.false. ! avoid redo screening inside batch loop


C ************  END EXCHANGE CONTRIBUTION  ****************

C --- Diagonals updating. It only makes sense if Nscreen > 0

               If (Update .and. Nscreen .gt. 0) Then

                  CALL CWTIME(TCS1,TWS1)
C ---------------------------------------------------------------------
C --- update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2
C
C --- subtraction is done in the 1st reduced set
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
C --------------------------------------------------------------------
C --- First half Active transformation  Lvb,J = sum_a  C(v,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCINT1,TWINT1)

C --- Set up the skipping flags
C --- The memory used before for the full-dimension AO-vectors
C ---     is now re-used to store half and full transformed
C ---     vectors in the active space
C -------------------------------------------------------------

               kMOs = 1  !
               nMOs = 1  ! Active MOs (1st set)

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                          JSYM,iSwap,IREDC,nMOs,kMOs,Ash,
     &                          Laq(1),DoRead)

               if (irc.ne.0) then
                  RETURN
               endif


C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
C --------------------------------------------------------------------
               If (JSYM.eq.1) Then

                  Do iSymb=1,nSym

                     NAv = nAorb(iSymb)

                     If(NAv.gt.0)Then

                      Do JVC=1,JNUM


                       CALL DGEMM_Tri('N','T',NAv,NAv,NBAS(iSymb),
     &                            One,Laq(1)%SB(iSymb)%A3(:,:,JVC),NAv,
     &                                Ash(1)%SB(iSymb)%A2,NAv,
     &                           Zero,Lxy%SB(iSymb)%A2(:,JVC),NAv)

                      End Do

                     EndIf

                  End Do

               Else

C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
C --------------------------------------------------------------------
                  Do iSymb=1,nSym

                     iSymv = MulD2h(JSYM,iSymb)
                     NAv = nAorb(iSymv)
                     NAw = nAorb(iSymb) ! iSymb=iSymw

                     If(NAv*NAw.ne.0 .and. iSymv.lt.iSymb)Then

                      Do JVC=1,JNUM


                       CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),
     &                            One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                                Ash(1)%SB(iSymb)%A2,NAw,
     &                           Zero,Lxy%SB(iSymv)%A2(:,JVC),NAv)

                      End Do

                     EndIf

                  End Do

               EndIf
C
C
C *************** EVALUATION OF THE (WA|XY) INTEGRALS ***********

               DoTraInt = JRED.eq.myJRED2.and.iBatch.eq.nBatch

               CALL CHO_eval_waxy(irc,Scr,Laq(1),Lxy,ipInt,nAorb,
     &                            JSYM,JNUM,DoTraInt)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  RETURN
               endif

               Call Deallocate_SBA(Lxy)
               Call Deallocate_SBA(Laq(1))

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               add  = .true.
               Call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,
     &                           [ipFLT],Frs,mode,add)
            EndIf

C --- free memory
            Call mma_deallocate(Lrs)

            If(JSYM.eq.1)Then
              Call mma_deallocate(Frs)
              Call mma_deallocate(Drs)
            EndIf

999         Continue

C --- Screening control section
            DoScreen = kscreen.eq.Nscreen

            if (.not.DoScreen) then
                kscreen = kscreen + 1
            else
                kscreen = 1
            endif

#if defined (_MOLCAS_MPP_)
            If (nProcs.gt.1 .and. Update .and. DoScreen
     &          .and. Is_Real_Par()) Then
               Call GaDsum(DiagJ,nnBSTR(JSYM,1))
               Call Daxpy_(nnBSTR(JSYM,1),xone,DiagJ,1,
     &                    Diag(1+iiBstR(JSYM,1)),1)
               Call Fzero(DiagJ,nnBSTR(JSYM,1))
            EndIf
C--- Need to activate the screening to setup the contributing shell
C--- indeces the first time the loop is entered .OR. whenever other nodes
C--- have performed screening in the meanwhile
            If (nProcs.gt.1 .and. .not.DoScreen .and. nVrs.eq.0
     &          .and. Is_Real_Par()) Then
               ntv0=ntv0+1
               DoScreen = (JRED.lt.myJRED1 .or. ntv0.ge.Nscreen)
               if (DoScreen) ntv0=0
            EndIf
#endif

         END DO   ! loop over red sets

         Call Deallocate_twxy(Scr)

      END DO  ! loop over JSYM

* --- Accumulate Coulomb and Exchange contributions
      Do jDen=1,nDen

        Do iSym=1,nSym

         ipFX = ipFLT(jDen) + ISTLT(iSym)
         ipKX = ipKLT(jDen) + ISTLT(iSym)

         Do iaSh=1,nShell

            ioffa = kOffSh(iaSh,iSym)

c --- ibSh < iaSh
c ---------------
            Do ibSh=1,iaSh-1

               iShp = iaSh*(iaSh-1)/2 + ibSh

               iOffAB = nnBfShp(iShp,iSym)

               ioffb = kOffSh(ibSh,iSym)

               Do ib=1,nBasSh(iSym,ibSh)

                Do ia=1,nBasSh(iSym,iaSh)

                  iab = nBasSh(iSym,iaSh)*(ib-1) + ia

                  jKX = ipKX - 1 + iOffAB + iab

                  iag = ioffa + ia
                  ibg = ioffb + ib

                  jFX = ipFX - 1 + iTri(iag,ibg)

                  Work(jFX) = Work(jFX) + Work(jKX)

                End Do

               End Do

            End Do

c --- ibSh = iaSh
c ---------------
            iShp = iaSh*(iaSh+1)/2

            iOffAB = nnBfShp(iShp,iSym)

            Do ib=1,nBasSh(iSym,iaSh)

             Do ia=ib,nBasSh(iSym,iaSh)

               iab = ia*(ia-1)/2 + ib

               jKX = ipKX - 1 + iOffAB + iab

               iag = ioffa + ia
               ibg = ioffa + ib

               jFX = ipFX - 1 + iag*(iag-1)/2 + ibg

               Work(jFX) = Work(jFX) + Work(jKX)

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
#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par())
     &    Call mma_deallocate(DiagJ)
#endif
      Call mma_deallocate(Diag)


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1


*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky RASSCF timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'SCREENING OVERHEAD               '
     &                           //'         ',tscrn(1),tscrn(2)
         Write(6,'(2x,A26,2f10.2)')'MO HALF-TRANSFORM VECTORS        '
     &                           //'         ',tmotr(1),tmotr(2)
         Write(6,'(2x,A26,2f10.2)')'EXCHANGE                         '
     &                           //'         ',texch(1),texch(2)
         Write(6,'(2x,A26,2f10.2)')'(PU|VX) INTEGRALS                '
     &                           //'         ',tintg(1),tintg(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

c Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in CASSCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      WRITE(6,'(6X,A)')'***** INACTIVE FOCK MATRIX ***** '
      DO ISYM=1,NSYM
        ISFI=ipFI+ISTLT(ISYM)
        IF( NBAS(ISYM).GT.0 ) THEN
          WRITE(6,'(6X,A)')
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call TRIPRT('','',Work(ISFI),NBAS(ISYM))
        ENDIF
      END DO
      IF(DoActive)THEN
        WRITE(6,'(6X,A)')
        WRITE(6,'(6X,A)')'***** ACTIVE FOCK MATRIX ***** '
        DO ISYM=1,NSYM
         IF( NBAS(ISYM).GT.0 ) THEN
           ISFA=ipFA+ISTLT(ISYM)
           WRITE(6,'(6X,A)')
           WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
           call TRIPRT('','',Work(ISFA),NBAS(ISYM))
         ENDIF
        END DO
      END IF

      endif

#endif

      Return
      END
