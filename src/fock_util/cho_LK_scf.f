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
      SUBROUTINE CHO_LK_SCF(rc,nDen,ipFLT,ipKLT,nForb,nIorb,
     &                         ipPorb,ipPLT,FactXI,nSCReen,dmpk,dFmat)

**********************************************************************
*  Author : F. Aquilante
*
C *************** INACTIVE AO-BASIS FOCK MATRIX **********************
C
C   F(ab) = sum_J  Lab,J * V(J)  -  sum_Jk  Lka,J * Lkb,J
C
C   Exchange term computed using the Local-exchange (LK) algorithm
C
**********************************************************************
C
C      V(J) = sum_gd  Lgd,J * Ptot(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index   belonging to (Frozen+Inactive)
C
**********************************************************************
      use ChoArr, only: nBasSh, nDimRS
      use ChoSwp, only: nnBstRSh, iiBstRSh, InfVec, IndRed
      use Data_Structures, only: NDSBA_Type, Allocate_NDSBA,
     &                           Deallocate_NDSBA
      use Data_Structures, only: Allocate_L_Full, Deallocate_L_Full,
     &                           L_Full_Type

#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit Real*8 (a-h,o-z)

      Type (NDSBA_type) DiaH
      Type (L_Full_Type) L_Full

      Integer   rc,nDen
      Integer   ipOrb(8,2),nOrb(8,2)
      Integer   ISTLT(8),ISTSQ(8),kOff(8,2)
      Real*8    tread(2),tcoul(2),texch(2)
      Real*8    tscrn(2),tmotr(2)
      Real*8    FactXI,dmpk,dFmat,tau(2),thrv(2)
      Integer   ipPLT(nDen),ipFLT(nDen),ipKLT(nDen)
      Integer   ipPorb(nDen)
      Integer   nForb(8,nDen),nIorb(8,nDen)
#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   DoScreen
      Character*50 CFmt
      Character(LEN=10), Parameter :: SECNAM = 'CHO_LK_SCF'
#include "chotime.fh"
#include "choscreen.fh"
#include "real.fh"

      Real*8, Parameter :: xone = -One, FactCI = one

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "warnings.fh"
      Logical add
      Character(LEN=6) mode

      Real*8 LKThr
      Real*8, External:: Cho_LK_ScreeningThreshold
      Integer, External:: Cho_LK_MaxVecPerBatch
      Integer, External:: ip_of_Work

      Real*8, Allocatable:: Lrs(:,:), Drs(:), Frs(:), VJ(:)
      Integer, Allocatable:: iShp_rs(:), Indx(:,:)

      Integer, Allocatable:: nnBfShp(:,:), ipLab(:), kOffSh(:,:)
      Real*8, Allocatable:: SvShp(:,:), Diag(:), AbsC(:), SumAClk(:,:),
     &                      Ylk(:,:), MLk(:,:), Faa(:), Fia(:)
#if defined (_MOLCAS_MPP_)
      Real*8, Allocatable:: DiagJ(:)
#endif
************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
****** next is a trick to save memory. Memory in "location 2" is used
******      to store this offset array defined later on
      iOffShp(i,j) = iiBstRSh(i,j,2)
************************************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in SCF-debug
#endif

      IREDC= -1  ! unknwn reduced set

      iLoc = 3 ! use scratch location in reduced index arrays

      If (nDen.ne.1 .and. nDen.ne.2) then
         write(6,*)SECNAM//'Invalid parameter nDen= ',nDen
         call abend()
      EndIf


        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        ! 1 --> CPU   2 --> Wall
        tread(:) = zero  !time read/transform vectors
        tcoul(:) = zero  !time for computing Coulomb
        texch(:) = zero  !time for computing Exchange
        tmotr(:) = zero  !time for the half-transf of vectors
        tscrn(:) = zero  !time for screening overhead

C ==================================================================

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
      nT1=0
      DO jDen=1,nDen

         DO ISYM=1,NSYM

            ipOrb(iSym,jDen) = ipPorb(jDen) + ISTSQ(iSym)

            nOrb(iSym,jDen)  = nForb(iSym,jDen)+nIorb(iSym,jDen)

            kOff(iSym,jDen) = nnO

            nnO = nnO + nOrb(iSym,jDen)

         END DO

         nT1 = nT1 + nnO*(2-jDen) ! tot # alpha orbitals

      END DO

      nT2 = nnO - nT1

C --- Define max number of vectors to be treated in core at once

      MaxVecPerBatch=Cho_LK_MaxVecPerBatch()

C --- Define the basic screening threshold

      LKThr=Cho_LK_ScreeningThreshold(dFmat)

C --- Adjust the damping according to the Abs(Max offDiag FMOmat)
Ctbp, may 2013: adjustment moved to Cho_LK_ScreeningThreshold
      fcorr = dmpk
Ctbp  If (dFmat.gt.zero ) Then
Ctbp     If ( dFmat.lt.1.0d3*LKThr ) Then
Ctbp        fcorr=dmpk*1.0d-2
Ctbp     EndIf
Ctbp     If (dFmat.le.LKThr) fcorr=fcorr*1.0d-2
Ctbp  EndIf

      tau(1) = (LKThr/Max(1,nT1))*fcorr ! screening alpha Fock matrix
      tau(2) = (LKThr/Max(1,nT2))*fcorr ! screening beta Fock matrix

      MaxRedT=MaxRed
      Call GAIGOP_SCAL(MaxRedT,'+')
      If (Estimate) Then
         Do i=1,nDen
            tau(i)=tau(i)/MaxRedT
         End Do
      EndIf

c     xtau(1) = sqrt(tau(1))
c     xtau(2) = sqrt(tau(2))

C --- Vector MO transformation screening thresholds
      NumVT=NumChT
      Call GAIGOP_SCAL(NumVT,'+')
      thrv(1) = ( sqrt(LKThr/(Max(1,nT1)*NumVT)) )*fcorr
      thrv(2) = ( sqrt(LKThr/(Max(1,nT2)*NumVT)) )*fcorr

      CALL mma_allocate(DIAG,NNBSTRT(1),Label='DIAG')
      DIAG(:)=0.0D0

#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par()) Then
         NNBSTMX=0
         Do i=1,nSym
            NNBSTMX = Max(NNBSTMX,NNBSTR(i,1))
         End Do
         Call mma_allocate(DiagJ,NNBSTMX,Label='DiagJ')
         DiagJ(:)=0.0D0
      EndIf
#endif

C *************** Read the diagonal integrals (stored as 1st red set)
      If (Update) CALL CHO_IODIAG(DIAG,2) ! 2 means "read"

c --- allocate memory for sqrt(D(a,b)) stored in full (squared) dim
      Call Allocate_NDSBA(DiaH,nBas,nBas,nSym)
      DiaH%A0(:)=Zero

c --- allocate memory for the abs(C(l)[k])
      Call mma_allocate(AbsC,MaxB,Label='AbsC')
      ipAbs=ip_of_Work(AbsC(1))

c --- allocate memory for the Y(l)[k] vectors
      Call mma_allocate(Ylk,MaxB,nnO,Label='Ylk')

c --- allocate memory for the ML[k] list of largest elements
c --- in significant shells
      Call mma_allocate(MLk,nShell,nnO,Label='MLk')

c --- allocate memory for the list of  S:= sum_l abs(C(l)[k])
c --- for each shell
      Call mma_allocate(SumAClk,nShell,nnO,Label='SumAClk')

c --- allocate memory for the Index array
      Call mma_allocate(Indx,[0,nShell],[1,nnO],Label='Indx')

c --- allocate memory for ipLab
      Call mma_allocate(ipLab,nShell,Label='ipLab')

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
      Fia(:)=Zero
      Faa(:)=Zero

C *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
      Do jDen=1,nDen
         Do kSym=1,nSym
            Do jK=1,nOrb(kSym,jDen)
               jK_a = jK + kOff(kSym,jDen)

               ipMO = ipOrb(kSym,jDen)
     &              + nBas(kSym)*(jK-1)

               Do iaSh=1,nShell

                  ipMsh = ipMO + kOffSh(iaSh,kSym)

                  SKsh=zero
                  Do ik=0,nBasSh(kSym,iaSh)-1
                     SKsh = SKsh + Work(ipMsh+ik)**2
                  End Do

                  SumAClk(iaSh,jK_a) = SKsh

               End Do

            End Do
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

C *** Compute Shell pair Offsets   iOffShp(iSyma,iShp)
         JNUM=1
         Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,
     &                        Memory=LFULL)

C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- half-transformed for a given MO-index k  --> La,[k]
C ------------------------------------------------------------------
         mTvec = 0  ! mem for storing the half-transformed vec

         do l=1,nSym
            k=Muld2h(l,JSYM)
            Mmax = 0
            do jDen=1,nDen
               Mmax = Max(Mmax,nForb(k,jDen)+nIorb(k,jDen))
            end do
            if (Mmax.gt.0) mTvec = Max(mTvec,nBas(l))
         end do

         mTvec=Max(mTvec,1)

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
#if defined (_MOLCAS_MPP_)
         myJRED1=JRED1 ! first red set present on this node
         ntv0=0
#endif

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
               Write(6,*) 'nVrs=',nVrs
               call Abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED)
c           !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                        '   rc= ',irc
              call Abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then

               Call mma_allocate(Drs,nRS,Label='Drs')
               Call mma_allocate(Frs,nRS,Label='Frs')
               Drs(:)=Zero
               Frs(:)=Zero

            EndIf

            Call mma_maxDBLE(LWORK)

            nVec = min(LWORK/(nRS+mTvec+LFULL),min(nVrs,MaxVecPerBatch))

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) ' LWORK= ',LWORK
               WRITE(6,*) ' jsym= ',jsym
               WRITE(6,*) ' min. mem. need= ',nRS+mTvec+LFULL
               WRITE(6,*) ' nRS = ',nRS
               WRITE(6,*) ' mTvec = ',mTvec
               WRITE(6,*) ' LFULL = ',LFULL
               CALL Quit(_RC_MEMORY_ERROR_)
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               add  = .false.
               nMat=1
               Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                           [ipPLT],Drs,mode,add)
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
                  rc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               If(JSYM.eq.1)Then
C ************ (alpha+beta) COULOMB CONTRIBUTION  ****************
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * P(rs)
C==========================================================
C
                  CALL CWTIME(TCC1,TWC1)

                  Call mma_allocate(VJ,JNUM,Label='VJ')

                  CALL DGEMV_('T',nRS,JNUM,
     &                       ONE,Lrs,nRS,
     &                       Drs,1,ZERO,VJ,1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  Fact = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                       FactCI,Lrs,nRS,
     &                       VJ,1,Fact,Frs,1)


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

                  Call Fzero(DIAG(1+iiBstR(jSym,1)),NNBSTR(jSym,1))

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
               Call GetMem('ChoT','Allo','Real',ipChoT,mTvec*nVec)
               Call Allocate_L_Full(L_Full,nShell,iShp_rs,JNUM,JSYM,
     &                              nSym)
               ipLF = ip_of_Work(L_Full%A0(1))

               CALL CWTIME(TCX1,TWX1)

C *** Reorder vectors to Full-dimensions
C ***
C *** Vectors are returned in the storage LaJ,b with the restriction:
C ***
C ***    Sym(a).ge.Sym(b)
C ***
C *** and blocked in shell pairs

               L_Full%A0(:)=Zero

               CALL CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,ipLF,
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
     &                               DIAH,DIAG)


                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               ENDIF


               Do jDen=1,nDen

                  Do kSym=1,nSym

                     lSym=MulD2h(JSYM,kSym)

                     Do jK=1,nOrb(kSym,jDen)
                        jK_a = jK + kOff(kSym,jDen)

                      CALL FZero(Work(ipChoT),nBas(lSym)*JNUM)

                      ipMO = ipOrb(kSym,jDen)
     &                     + nBas(kSym)*(jK-1)

                      IF (DoScreen) THEN

                        CALL CWTIME(TCS1,TWS1)
C------------------------------------------------------------------
C --- Setup the screening
C------------------------------------------------------------------
                        Do ik=0,nBas(kSym)-1
                           AbsC(1+ik) = abs(Work(ipMO+ik))
                        End Do

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
     &                             ONE,DiaH%SB(lSym,kSym)%A2,n1,
     &                                 AbsC,1,
     &                            ZERO,Ylk(1,jK_a),1)

C --- List the shells present in Y(l)[k] by the largest element
                        Do ish=1,nShell
                           YshMax=zero
                           Do ibs=1,nBasSh(lSym,ish)
                              YshMax = Max(YshMax,
     &                               Ylk(koffSh(ish,lSym)+ibs,jK_a))
                           End Do
                           MLk(ish,jK_a) = YshMax
                        End Do

C --- Sort the list ML[k]
                        numSh=0  ! # of significant shells
                        jml=1

                        Do ish=1,nShell
                           Indx(ish,jK_a) = ish
                        End Do

                        Do while (jml.le.nShell)

                           YMax=MLk(jml,jK_a)
                           jmlmax=jml

                           Do iml=jml+1,nShell  ! get the max
                              If (MLk(iml,jK_a).gt.YMax) then
                                 YMax = MLk(iml,jK_a)
                                 jmlmax = iml
                              Endif
                           End Do

                           If (jmlmax.ne.jml) then  ! swap positions
                              xTmp = MLk(jml,jK_a)
                              iTmp = Indx(jml,jK_a)
                              MLk(jml,jK_a) = YMax
                              Indx(jml,jK_a) = Indx(jmlmax,jK_a)
                              MLk(jmlmax,jK_a) = xTmp
                              Indx(jmlmax,jK_a) = iTmp
                           Endif

c --- Exact bounds (quadratic scaling of the MO transformation)
ctbp, may 2013: reintroduce exact bound
                            If(MLk(jml,jK_a)*MLk(1,jK_a)
     &                                          .ge. tau(jDen))then
c
c --- Here we use a non-exact bound for the exchange matrix
c --- The positive definiteness of the exchange matrix
c --- combined with the structure of the density matrix makes this
c --- bound acceptable and likely to be almost exact for what concerns
c --- the exchange energy
ctbp                       If ( MLk(jml,jK_a) .ge. xtau(jDen) ) then
                              numSh = numSh + 1
                           else
                              jml=nShell  ! exit the loop
                           endif

                           jml=jml+1

                        End Do

                        Indx(0,jk_a)=numSh

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

                         iaSh = Indx(iSh,jK_a)

                         iOffSha = kOffSh(iaSh,lSym)

                         ipLab(iaSh) = ipChoT + iOffSha*JNUM

                         ibcount=0

                         Do ibSh=1,nShell

                            iOffShb = kOffSh(ibSh,kSym)

                            iShp = iTri(iaSh,ibSh)

                            If ( iShp_rs(iShp)<=0) Cycle

                            If ( nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                         nBasSh(lSym,iaSh)*
     &                         nBasSh(kSym,ibSh) .gt. 0
     &                         .and. sqrt(abs(SumAClk(ibSh,jK_a)*
     &                                    SvShp(iShp_rs(iShp),1) ))
     &                         .ge. thrv(jDen) )Then

                               ibcount = ibcount + 1

                               IF (lSym.ge.kSym) Then

                                  jOff = iOffShp(lSym,iShp_rs(iShp))
                                  If (iaSh<ibSh) jOff = jOff +
     &                            nBasSh(lSym,ibSh)*nBasSh(kSym,iaSh)

C ---  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
C ---------------------------------------
                                  Mode(1:1)='N'
                                  n1 = nBasSh(lSym,iaSh)*JNUM
                                  n2 = nBasSh(kSym,ibSh)

                               Else   ! lSym < kSym

                                  jOff = iOffShp(kSym,iShp_rs(iShp))
                                  If (ibSh<iaSh) jOff = jOff +
     &                            nBasSh(kSym,iaSh)*nBasSh(lSym,ibSh)

C ---  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
C ---------------------------------------
                                  Mode(1:1)='T'
                                  n1 = nBasSh(kSym,ibSh)
                                  n2 = JNUM*nBasSh(lSym,iaSh)

                               Endif

                               CALL DGEMV_(Mode(1:1),n1,n2,
     &                                  ONE,Work(ipLF+jOff*JNUM),n1,
     &                                      Work(ipMO+ioffShb),1,
     &                                  ONE,Work(ipLab(iaSh)),1)

                            EndIf

                          End Do

c --- The following re-assignement is used later on to check if the
c --- iaSh vector LaJ[k] can be neglected because identically zero

                          If (ibcount==0) ipLab(iaSh) = ipAbs

                      End Do

                      CALL CWTIME(TCT2,TWT2)
                      tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                      tmotr(2) = tmotr(2) + (TWT2 - TWT1)

C --- Prepare the J-screening

                      CALL CWTIME(TCS1,TWS1)


                      Do iSh=1,Indx(0,jk_a)

                         iaSh = Indx(iSh,jK_a)

                         iaSkip= Min(1,Max(0,
     &                           abs(ipLab(iaSh)-ipAbs))) ! = 1 or 0

                         If (iaSkip==0) Cycle

                         IF (lSym.ge.kSym) Then

C ---  Faa,[k] = sum_J  (LaJ[k])**2
C ----------------------------------
                            Inc=nBasSh(lSym,iaSh)
                            n1 = 1

                         Else   ! lSym < kSym

C ---  Faa,[k] = sum_J  (LJa[k])**2
C ----------------------------------
                            Inc=1
                            n1 = JNUM

                         End If

                         Tmp=Zero
                         Do ia=1,nBasSh(lSym,iaSh)
                            ipLaa = ipLab(iaSh) + n1*(ia-1)
                            Fia(ia)=DDot_(JNUM,Work(ipLaa),Inc,
     &                                         Work(ipLaa),Inc)
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

                      Do lSh=1,Indx(0,jK_a)

                         iaSh = Indx(lSh,jK_a)

                         iaSkip= Min(1,Max(0,
     &                           abs(ipLab(iaSh)-ipAbs))) ! = 1 or 0

                         iOffSha = kOffSh(iaSh,lSym)

                         mSh = 1

                         Do while (mSh.le.Indx(0,jK_a))

                            ibSh = Indx(mSh,jK_a)

                            ibSkip = Min(1,Max(0,
     &                                     abs(ipLab(ibSh)-ipAbs)))

                            iShp = iTri(iaSh,ibSh)

                            iOffShb = kOffSh(ibSh,lSym)

                            iOffAB = nnBfShp(iShp,lSym)

                            ipKI = ipKLT(jDen) + ISTLT(lSym) + iOffAB

                            xFab = sqrt(abs(Faa(iaSh)*Faa(ibSh)))

                            If ( MLk(lSh,jK_a)*MLk(mSh,jK_a)
     &                           .lt. tau(jDen) ) Then


                                mSh = Indx(0,jK_a)  ! skip the rest


                            ElseIf(iaSh.eq.ibSh
     &                             .and.xFab.ge.tau(jDen)/MaxRedT
     &                             .and.iaSkip.eq.1)Then

                              IF (lSym.ge.kSym) Then


C ---  F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(a,J)[k] * L(b,J)[k]
C -------------------------------------------------------------------
                                 Mode(1:1)='N'
                                 Mode(2:2)='T'
                                 n1 = nBasSh(lSym,iaSh)
                                 nBs = nBasSh(lSym,iaSh)

                              ELSE   ! lSym < kSym

C ---  F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(J,a)[k] * L(J,b)[k]
C -------------------------------------------------------------------
                                 Mode(1:1)='T'
                                 Mode(2:2)='N'
                                 n1 = JNUM
                                 nBs = nBasSh(lSym,iaSh)

                              End If

                              CALL DGEMM_Tri(Mode(1:1),Mode(2:2),
     &                        nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                              -FActXI,Work(ipLab(iaSh)),n1,
     &                                      Work(ipLab(ibsh)),n1,
     &                                  ONE,Work(ipKI),nBs)

                            ElseIf (iaSh.gt.ibSh
     &                                .and.xFab.ge.tau(jDen)/MaxRedT
     &                                 .and. iaSkip*ibSkip.eq.1) Then


                              IF (lSym.ge.kSym) Then

C ---  F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(a,J)[k] * L(b,J)[k]
C -------------------------------------------------------------------
                                 Mode(1:1)='N'
                                 Mode(2:2)='T'
                                 nBs  = nBasSh(lSym,iaSh)
                                 n1 = nBasSh(lSym,iaSh)
                                 n2 = nBasSh(lSym,ibSh)

                              ELSE   ! lSym < kSym

C ---  F(a,b)[k] = F(a,b)[k] - FactXI * sum_J  L(J,a)[k] * L(J,b)[k]
C -------------------------------------------------------------------
                                 Mode(1:1)='T'
                                 Mode(2:2)='N'
                                 n1 = JNUM
                                 n2 = JNUM
                                 nBs = nBasSh(lSym,iaSh)

                              End If

                              CALL DGEMM_(Mode(1:1),Mode(2:2),
     &                         nBasSh(lSym,iaSh),nBasSh(lSym,ibSh),JNUM,
     &                               -FactXI,Work(ipLab(iaSh)),n1,
     &                                       Work(ipLab(ibsh)),n2,
     &                                   ONE,Work(ipKI),nBs)

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

               Call Deallocate_L_Full(L_Full)
               Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
               DoScreen=.false. ! avoid redo screening inside batch loop

C --- Diagonals updating. It only makes sense if Nscreen > 0

               If (Update .and. Nscreen.gt.0) Then

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

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               add  = .true.
               nMat = 1
               Do jDen = 1, nDen
                  Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                           [ipFLT(jDen)],Frs,mode,add)
               End Do
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

      END DO   !loop over JSYM

* --- Accumulate Coulomb and Exchange contributions
      Do jDen=1,nDen

         Do iSym=1,nSym

            ipFI = ipFLT(jDen) + ISTLT(iSym)
            ipKI = ipKLT(jDen) + ISTLT(iSym)

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

                     jK = ipKI - 1 + iOffAB + iab

                     iag = ioffa + ia
                     ibg = ioffb + ib

                     jF = ipFI - 1 + iTri(iag,ibg)

                     Work(jF) = Work(jF) + Work(jK)

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

                  jK = ipKI - 1 + iOffAB + iab

                  iag = ioffa + ia
                  ibg = ioffa + ib

                  jF = ipFI - 1 + iag*(iag-1)/2 + ibg

                  Work(jF) = Work(jF) + Work(jK)

                End Do

               End Do


            End Do

         End Do

      End Do


      Call mma_deallocate(Fia)
      Call mma_deallocate(Faa)
      Call mma_deallocate(SvShp)
      Call mma_deallocate(iShp_rs)
      Call mma_deallocate(nnBfShp)
      Call mma_deallocate(kOffSh)
      Call mma_deallocate(ipLab)
      Call mma_deallocate(Indx)
      Call mma_deallocate(SumAClk)
      Call mma_deallocate(MLk)
      Call mma_deallocate(Ylk)
      Call mma_deallocate(AbsC)
      Call Deallocate_NDSBA(DiaH)
#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par())
     &    CALL mma_deallocate(DiagJ)
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
      Write(6,CFmt)'Cholesky SCF timing from '//SECNAM
      Write(6,CFmt)'------------------------------------'
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
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


c Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in SCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      WRITE(6,'(6X,A)')'***** FOCK MATRIX AO-BASIS ***** '
      Do jDen=1,nDen
        if(nDen.eq.2) Then
          if(jden.eq.1) WRITE(6,'(6X,A)')'******** ALPHA SPIN ******** '
          if(jden.eq.2) WRITE(6,'(6X,A)')'******** BETA SPIN ********* '
        endif
        DO ISYM=1,NSYM
           ISFI=ipFLT(jDen)+ISTLT(ISYM)
           IF( NBAS(ISYM).GT.0 ) THEN
             WRITE(6,'(6X,A)')
             WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES: ',ISYM
             call TRIPRT('','',Work(ISFI),NBAS(ISYM))
           ENDIF
        END DO
      END DO

      endif

#endif

      rc  = 0


      Return
      END
