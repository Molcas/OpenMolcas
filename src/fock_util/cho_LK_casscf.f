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
     &             FactXI,nFIorb,nAorb,nChM,ipAsh,DoActive,
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


      Implicit Real*8 (a-h,o-z)

      Integer   rc,ipLxy(8),ipScr(8,8),ipDIAH(1)
      Integer   ipLpq(8),ipDLT(2),ipFLT(2),ipKLT(2)
      Integer   iSkip(8),kOff(8,2),nnA(8,8)
      Integer   ISTLT(8),ISTSQ(8),ISSQ(8,8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2)
      Real*8    tmotr(2),tscrn(2)
      Integer   ipAsh(2),ipAorb(8,2),ipDab(2),ipFab(2),ipDD(2)
      Integer   nFIorb(8),nAorb(8),nChM(8)
      Logical   Debug,timings,DoRead,DoTraInt,DoActive,DoScreen
      Logical   Estimate,Update
      Real*8    FactC(2),FactX(2),tau(2),xtau(2),thrv(2)
      Real*8    dmpK, dFmat,ExFac
      Character*50 CFmt
      Character*13 SECNAM
      Parameter (SECNAM = 'CHO_LK_CASSCF')
      COMMON    /CHOTIME /timings
      COMMON    /CHOSCREEN/ Estimate,Update

      parameter (DoRead = .false. )
      parameter (zero = 0.0D0, one = 1.0D0, xone = -1.0D0)

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "para_info.fh"

      parameter ( N2 = InfVec_N2 )
      Logical add
      Character*6 mode
      Integer  Cho_F2SP
      External Cho_F2SP
      Real*8   LKThr
      Real*8   Cho_LK_ScreeningThreshold
      External Cho_LK_ScreeningThreshold
      Integer  Cho_LK_MaxVecPerBatch
      External Cho_LK_MaxVecPerBatch

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
******
      NBASSH(I,J)=IWORK(ip_NBASSH-1+NSYM*(J-1)+I)
******
      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
******
      nnBfShp(j,i) = iWork(ip_nnBfShp-1+nnShl_tot*(i-1)+j)
******
      ipLab(i) = iWork(ip_Lab+i-1)
******
      kOffSh(i,j) = iWork(ip_kOffSh+nShell*(j-1)+i-1)
******
      iShp_rs(i) = iWork(ip_iShp_rs+i-1)
******
      SvShp(i) = Work(ip_SvShp+i-1)
****** next is a trick to save memory. Memory in "location 2" is used
******      to store this offset array defined later on
      iOffShp(i,j) = iWork(ip_iiBstRSh+nSym*nnShl-1+nSym*(j-1)+i)
************************************************************************


#ifdef _DEBUG_
c      Debug=.true.
      Debug=.false.! to avoid double printing in CASSCF-debug
#else
      Debug=.false.
#endif


      DoTraInt = .false.
      IREDC = -1  ! unknown reduced set in core

      ipDLT(1) = ipDI    ! some definitions
      ipDLT(2) = ipDA1
      ipFLT(1) = ipFI
      ipFLT(2) = ipFA

      FactC(1) = one
      FactX(1) = FactXI*ExFac
      FactC(2) = one
      FactX(2) = -0.5D0*ExFac

      nDen=2  ! inactive and active density, respectively
      If (.not.DoActive) nDen=1

        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time read/transform vectors
           tcoul(i) = zero  !time for computing Coulomb
           texch(i) = zero  !time for computing Exchange
           tintg(i) = zero  !time for computing (tw|xy) integrals
           tmotr(i) = zero  !time for the half-transf of vectors
           tscrn(i) = zero  !time for screening overhead
        end do

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

      nnBSQ=0
      DO LSYM=1,NSYM
         DO KSYM=LSYM,NSYM
            ISSQ(KSYM,LSYM) = nnBSQ
            ISSQ(LSYM,KSYM) = nnBSQ ! symmetrization
            nnBSQ = nnBSQ + nBas(kSym)*nBas(lSym)
         END DO
      END DO

**************************************************

      nnO=0
      nIt=0
      DO jDen=1,nDen

         ipAorb(1,jDen)= ipAsh(jDen)

         kOff(1,jDen) = nnO

         DO ISYM=2,NSYM

            ipAorb(iSym,jDen) = ipAorb(iSym-1,jDen)
     &                        + nAorb(iSym-1)*nBas(iSym-1)*(2-jDen)
     &                        + (jDen-1)*nBas(iSym-1)**2

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

      CALL GETMEM('diagI','Allo','Real',ipDIAG,NNBSTRT(1))
      ipDD(1)=ipDIAG

#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par()) Then
         NNBSTMX=0
         Do i=1,nSym
            NNBSTMX = Max(NNBSTMX,NNBSTR(i,1))
         End Do
         CALL GETMEM('diagJ','Allo','Real',ipjDIAG,NNBSTMX)
         Call FZero(Work(ipjDIAG),NNBSTMX)
      EndIf
#endif

C *************** Read the diagonal integrals (stored as 1st red set)
      If (Update) CALL CHO_IODIAG(Work(ipDIAG),2) ! 2 means "read"

c --- allocate memory for sqrt(D(a,b)) stored in full (squared) dim
      CALL GETMEM('diahI','Allo','Real',ipDIAH(1),NNBSQ)
      CALL FZERO(Work(ipDIAH(1)),NNBSQ)

c --- allocate memory for the abs(C(l)[k])
      Call GetMem('absc','Allo','Real',ipAbs,MaxB)

c --- allocate memory for the Y(l)[k] vectors
      Call GetMem('yc','Allo','Real',ipY,MaxB*nnO)

c --- allocate memory for the ML[k] lists of largest elements
c --- in significant shells
      Call GetMem('MLk1','Allo','Real',ipML,nShell*nnO)

c --- allocate memory for the list of  S:= sum_l abs(C(l)[k])
c --- for each shell
      Call GetMem('SKsh','Allo','Real',ipSKsh,nShell*nnO)

c --- allocate memory for the Index arrays
      Call GetMem('Indx','Allo','Inte',ipIndx,(nShell+1)*nnO)

c --- allocate memory for ipLab
      Call GetMem('ip_Lab','Allo','Inte',ip_Lab,nShell)

c --- allocate memory for kOffSh
      Call GetMem('ip_kOffSh','Allo','Inte',ip_kOffSh,nShell*nSym)

c --- allocate memory for nnBfShp
      Call GetMem('ip_nnBfShp','Allo','Inte',ip_nnBfShp,nnShl_tot*nSym)

c --- allocate memory for iShp_rs
      Call GetMem('ip_iShp_rs','Allo','Inte',ip_iShp_rs,nnShl_tot)

c --- allocate memory for the shell-pair Frobenius norm of the vectors
      Call GetMem('ip_SvShp','Allo','Real',ip_SvShp,2*nnShl)


C *** Compute Shell Offsets ( MOs and transformed vectors)

      MxBasSh = 0

      Do iSyma=1,nSym

         LKsh=0

         Do iaSh=1,nShell    ! kOffSh(iSh,iSym)

            iWork(ip_kOffSh+nShell*(iSyma-1)+iaSh-1) = LKsh

            LKsh = LKsh + nBasSh(iSyma,iaSh)

            MxBasSh = Max(MxBasSh,nBasSh(iSyma,iaSh))

         End Do

      End Do

c --- allocate memory for the Diagonal of the Fock matrix
      Call GetMem('F(k)ss','Allo','Real',ipFk,MxBasSh+nShell)
      Call Fzero(Work(ipFk),MxBasSh+nShell)

C *** Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)
      Do jDen=1,nDen
         Do kSym=1,nSym

            nkOrb = nFIorb(kSym)*(2-jDen)
     &            + nChM(kSym)*(jDen-1)

            Do jK=1,nkOrb

               ipMO = (2-jDen)*(ipMSQ + ISTSQ(kSym)
     &              + nBas(kSym)*(jK-1))
     &              + (ipAorb(kSym,jDen)
     &              + nBas(kSym)*(jK-1))*(jDen-1)

               ipSk = ipSKsh + nShell*(kOff(kSym,jDen) + jK - 1)

               Do iaSh=1,nShell

                  ipMsh = ipMO + kOffSh(iaSh,kSym)

                  SKsh=zero
                  Do ik=0,nBasSh(kSym,iaSh)-1
                     SKsh = SKsh + Work(ipMsh+ik)**2
                  End Do

                  Work(ipSk+iaSh-1) = SKsh

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

               iWork(ip_nnBfShp - 1 + nnShl_tot*(iSyma-1)
     &        + iShp) = LKShp   ! nnBfShp(iShp,iSyma)

               LKShp = LKShp + nBasSh(iSyma,iaSh)*nBasSh(iSyma,ibSh)
     &                       - (1-Min((iaSh-ibSh),1))*nBasSh(iSyma,iaSh)
     &                       * (nBasSh(iSyma,iaSh) - 1)/2

             End Do

            End Do

         End Do

C *** Mapping shell pairs from the full to the reduced set

      Do iaSh=1,nShell
         Do ibSh=1,iaSh
            iShp = iaSh*(iaSh-1)/2 + ibSh
            iWork(ip_iShp_rs+iShp-1) = Cho_F2SP(iShp)
         End Do
      End Do


C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

        NumCV=NumCho(jSym)
        Call GAIGOP_SCAL(NumCV,'max')
        If (NumCV .lt. 1) GOTO 1000

C *** Compute Shell pair Offsets   iOffShp(iSyma,iShp)

         LFULL=0

         Do iaSh=1,nShell

          Do ibSh=1,iaSh

           iShp = iaSh*(iaSh-1)/2 + ibSh

           If (iShp_rs(iShp).gt.0) Then

            If (nnBstRSh(Jsym,iShp_rs(iShp),1).gt.0) Then

             Do iSymb=1,nSym

              iSyma=MulD2h(iSymb,Jsym)

              If (iSyma.ge.iSymb) Then

               iWork(ip_iiBstRSh + nSym*nnShl - 1
     &         + nSym*(iShp_rs(iShp)-1) + iSyma) = LFULL

                 LFULL = LFULL + nBasSh(iSyma,iaSh)*nBasSh(iSymb,ibSh)
     &        + Min(1,(iaSh-ibSh))*nBasSh(iSyma,ibSh)*nBasSh(iSymb,iaSh)

              EndIf

             End Do

            EndIf

           EndIf

          End Do

         End Do


C *** memory for the (wa|xy) integrals --- temporary array
         Mwaxy = 0
         Do iSymy=1,nSym
            iSymx=MulD2h(iSymy,JSYM)
            If (iSymx.le.iSymy) then
               Do iSyma=1,nSym
                  iSymw=MulD2h(iSyma,JSYM)
                  Mwaxy = Mwaxy + nAorb(iSymw)*nBas(iSyma)
     &                                        *nnA(iSymx,iSymy)
               End Do
             End If
         End Do

         Call GetMem('Mtmp','ALLO','REAL',ipItmp,Mwaxy)
         Call Fzero(Work(ipItmp),Mwaxy)

C *** setup pointers to the symmetry blocks of (wa|xy)
         Do i=1,nSym
            Do j=1,nSym
               ipScr(j,i) = ipItmp
            End Do
         End Do

         kScr=ipItmp
         Do iSymy=1,nSym
            iSymx=MulD2h(iSymy,JSYM)
            If (iSymx.le.iSymy) then
               Do iSyma=1,nSym
                  iSymw=MulD2h(iSyma,JSYM)
                  ipScr(iSymw,iSymx) = kScr
                  kScr = kScr + nAorb(iSymw)*nBas(iSyma)
     &                                      *nnA(iSymx,iSymy)
               End Do
             End If
         End Do

C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec = 0
         MxB=0
         do l=1,nSym
            k=Muld2h(l,JSYM)
            If ((nFIorb(k)+nChM(k)).gt.0) MxB=Max(MxB,nBas(l))
            mTvec = mTvec + nAorb(k)*nBas(l)
            If (k.le.l) mTvec = mTvec + nnA(k,l)
         end do

         LFMAX = Max(mTvec,LFULL) ! re-use memory for the active vec
         mTvec = Max(MxB,1) ! mem for storing half-transformed vec

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         iLoc = 3 ! use scratch location in reduced index arrays

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
         myJRED1=JRED1 ! first red set present on this node
         myJRED2=JRED2 ! last  red set present on this node

c --- entire red sets range for parallel run
         Call GAIGOP_SCAL(JRED1,'min')
         Call GAIGOP_SCAL(JRED2,'max')

         ntv0=0
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
             Do jden=1,nden
               Call GetMem('rsD','Allo','Real',ipDab(jden),nRS)
               Call GetMem('rsFC','Allo','Real',ipFab(jden),nRS)
               Call Fzero(Work(ipDab(jden)),nRS)
               Call Fzero(Work(ipFab(jden)),nRS)
             End Do
            EndIf

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec=min(LWORK/(nRS+mTvec+LFMAX),min(nVrs,MaxVecPerBatch))

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec+LFMAX
               WRITE(6,*) 'jsym= ',jsym
               WRITE(6,*) ' nRS = ',nRS
               WRITE(6,*) ' mTvec = ',mTvec
               WRITE(6,*) ' LFMAX = ',LFMAX
               rc = 33
               CALL QTrace()
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
            Call GetMem('ChoT','Allo','Real',ipChoT,mTvec*nVec)
            CALL GETMEM('FullV','Allo','Real',ipLF,LFMAX*nVec)

            If(JSYM.eq.1)Then
C --- Transform the densities to reduced set storage
               mode = 'toreds'
               add  = .false.
               Call play_casscf_sto(irc,iLoc,nDen,JSYM,ISTLT,ISSQ,
     &                                  ipDLT,ipDab,mode,add)
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

               CALL CHO_VECRD(Work(ipLrs),LREAD,JVEC,IVEC2,JSYM,
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  rc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               If(JSYM.eq.1)Then

                 CALL CWTIME(TCC1,TWC1)

                 ipVJ = ipChoT

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
     &                  ONE,Work(ipLrs),nRS,
     &                  Work(ipDab(jDen)),1,ZERO,Work(ipVJ),1)

C --- F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
C===============================================================

                   CALL DGEMV_('N',nRS,JNUM,
     &                 FactC(jDen),Work(ipLrs),nRS,
     &                 Work(ipVJ),1,Fact,Work(ipFab(jDen)),1)


                 End Do

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

                  Call Fzero(Work(ipDiag+iiBstR(jSym,1)),NNBSTR(jSym,1))

                  Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        ipL = ipLrs + nRS*(jvc-1)

                        Work(ipDiag+jrs-1) = Work(ipDiag+jrs-1)
     &                                  + Work(ipL+krs-1)**2

                     End Do

                  End Do

               EndIf

               CALL CWTIME(TCS2,TWS2)
               tscrn(1) = tscrn(1) + (TCS2 - TCS1)
               tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               CALL CWTIME(TCX1,TWX1)

C *** Reorder vectors to Full-dimensions
C ***
C *** Vectors are returned in the storage LaJ,b with the restriction:
C ***
C ***    Sym(a).ge.Sym(b)
C ***
C *** and blocked in shell pairs

               CALL FZero(Work(ipLF),LFULL*JNUM)
               CALL FZero(Work(ip_SvShp),2*nnShl)

               CALL CHO_getShFull(Work(ipLrs),lread,JNUM,JSYM,
     &                            IREDC,ipLF,Work(ip_SvShp),
     &                            iWork(ip_iShp_rs))

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)


               IF (DoScreen) THEN

                  CALL CWTIME(TCS1,TWS1)

c --- Compute DH(a,b)=sqrt(D(a,b)) from the updated diagonals.
c ---                              Only the symmetry blocks with
c ---                              compound symmetry JSYM are computed
c --------------------------------------------------------------------
                   mode = 'tosqrt'
                   ired1 = 1 ! location of the 1st red set
                   add  = .false.
                   nMat = 1
                   Call play_casscf_sto(irc,ired1,nMat,JSYM,ISTLT,ISSQ,
     &                                     ipDIAH,ipDD,mode,add)

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               ENDIF


               Do jDen=1,nDen

                 Do kSym=1,nSym

                    lSym=MulD2h(JSYM,kSym)

                    nkOrb = nFIorb(kSym)*(2-jDen)
     &                    + nChM(kSym)*(jDen-1)

                    Do jK=1,nkOrb

                     CALL FZero(Work(ipChoT),nBas(lSym)*JNUM)

                     ipMO = (2-jDen)*(ipMSQ + ISTSQ(kSym)
     &                    + nBas(kSym)*(jK-1))
     &                    + (ipAorb(kSym,jDen)
     &                    + nBas(kSym)*(jK-1))*(jDen-1)

                     ipYk = ipY + MaxB*(kOff(kSym,jDen) + jK - 1)

                     ipMLk = ipML + nShell*(kOff(kSym,jDen)
     &                            + jK - 1)

                     ipIndSh = ipIndx + (nShell+1)*(kOff(kSym,jDen)
     &                                + jK - 1)

                     ipSk = ipSKsh + nShell*(kOff(kSym,jDen) + jK - 1)


                     IF (DoScreen) THEN

                        CALL CWTIME(TCS1,TWS1)
C------------------------------------------------------------------
C --- Setup the screening
C------------------------------------------------------------------
                        ipDIH = ipDIAH(1) + ISSQ(lSym,kSym)

                        Do ik=0,nBas(kSym)-1
                           Work(ipAbs+ik) = abs(Work(ipMO+ik))
                        End Do

                        If (lSym.ge.kSym) Then
c --------------------------------------------------------------
C --- Y(l)[k] = sum_n  DH(l,n) * |C(n)[k]|
C===============================================================
                           nBs = Max(1,nBas(lSym))

                           CALL DGEMV_('N',nBas(lSym),nBas(kSym),
     &                                ONE,Work(ipDIH),nBs,
     &                                    Work(ipAbs),1,
     &                               ZERO,Work(ipYk),1)

                        Else
c --------------------------------------------------------------
C --- Y(l)[k] = sum_n  DH(n,l) * |C(n)[k]|
C===============================================================
                           nBs = Max(1,nBas(kSym))

                           CALL DGEMV_('T',nBas(kSym),nBas(lSym),
     &                                ONE,Work(ipDIH),nBs,
     &                                    Work(ipAbs),1,
     &                               ZERO,Work(ipYk),1)

                        EndIf

c            write(6,*)'Y(k)= ',(Work(ipYk+i-1),i=1,nBas(lSym))
c            write(6,*)'|C(k)|= ',(Work(ipAbs+i-1),i=1,nBas(kSym))

c            Call recprt('DH','',Work(ipDIH),nBas(lSym),nBas(kSym))

C --- List the shells present in Y(l)[k] by the largest element
                        Do ish=1,nShell
                           YshMax=zero
                           Do ibs=1,nBasSh(lSym,ish)
                              YshMax = Max(YshMax,
     &                          Work(ipYk+koffSh(ish,lSym)+ibs-1))
                           End Do
                           Work(ipMLk+ish-1) = YshMax
                        End Do

C --- Sort the lists ML[k]
                        Do ish=1,nShell
                           iWork(ipIndSh+ish) = ish
                        End Do

C **** Sort the list
                        numSh=0  ! # of significant shells
                        jml=1
                        Do while (jml.le.nShell)

                           YMax=Work(ipMLk+jml-1)
                           jmlmax=jml

                           Do iml=jml+1,nShell  ! get the max
                              If (Work(ipMLk+iml-1).gt.YMax) then
                                 YMax = Work(ipMLk+iml-1)
                                 jmlmax = iml
                              Endif
                           End Do

                           If(jmlmax.ne.jml) then  ! swap positions
                             xTmp = Work(ipMLk+jml-1)
                             iTmp = iWork(ipIndSh+jml)
                             Work(ipMLk+jml-1) = YMax
                             iWork(ipIndSh+jml) = iWork(ipIndSh+jmlmax)
                             Work(ipMLk+jmlmax-1) = xTmp
                             iWork(ipIndSh+jmlmax) = iTmp
                           Endif

c --- Exact bounds (quadratic scaling of the MO transformation)
c
c                           If(Work(ipMLk+jml-1)*Work(ipMLk)
c     &                                         .ge. tau(jDen))then
c
c --- Here we use a non-exact bound for the exchange matrix to achieve
c --- linear scaling. The positive definiteness of the exchange matrix
c --- combined with the structure of the density matrix makes this
c --- bound acceptable and likely to be almost exact for what concerns
c --- the exchange energy
                           If ( Work(ipMLk+jml-1) .ge. xtau(jDen) ) then
                             numSh = numSh + 1
                           else
                             jml=nShell  ! exit the loop
                           endif

                           jml=jml+1

                        End Do

                        iWork(ipIndSh) = numSh

c         write(6,*)'ord-ML(k)= ',(Work(ipMLk+i-1),i=1,nShell)
c         write(6,*)'Ind-ML(k)= ',(iWork(ipIndSh+i-1),i=1,nShell+1)
c         write(6,*)'lSym,kSym,jSym,jk,nShell,numSh= ',lSym,
c     &              kSym,jSym,jk,nShell,numSh

                        CALL CWTIME(TCS2,TWS2)
                        tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                        tscrn(2) = tscrn(2) + (TWS2 - TWS1)
C------------------------------------------------------------------
                     ENDIF    ! Screening setup


                     CALL CWTIME(TCT1,TWT1)

C --- Transform vectors for shells in the list ML[k]
C
C --- Screening based on the Frobenius norm: sqrt(sum_ij  A(i,j)^2)
C
C ---   || La,J[k] ||  .le.  || Lab,J || * || Cb[k] ||

                     IF (lSym.ge.kSym) Then


                         Do iSh=1,iWork(ipIndSh)

                            iaSh = iWork(ipIndSh+iSh)

                            iOffSha = kOffSh(iaSh,lSym)

                            iWork(ip_Lab+iaSh-1) = ipChoT + iOffSha*JNUM

                            ibcount=0

                            Do ibSh=1,nShell

                               iOffShb = kOffSh(ibSh,kSym)

                               iShp = iTri(iaSh,ibSh)

                               If ( iShp_rs(iShp) .gt. 0 ) Then

                                 If ( nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                              nBasSh(lSym,iaSh)*
     &                              nBasSh(kSym,ibSh) .gt. 0
     &                              .and. sqrt(abs(Work(ipSk+ibSh-1)*
     &                                         SvShp(iShp_rs(iShp)) ))
     &                              .ge. thrv(jDen) )Then

                                    ibcount = ibcount + 1

                                    jOff = iOffShp(lSym,iShp_rs(iShp)) -
     &                                     nBasSh(lSym,ibSh)*
     &                                     nBasSh(kSym,iaSh)*
     &                             Min(0,(iaSh-ibSh))/Max(1,(ibSh-iaSh))


C ---  LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
C ---------------------------------------

                                 CALL DGEMV_('N',nBasSh(lSym,iaSh)*JNUM,
     &                                        nBasSh(kSym,ibSh),
     &                                    ONE,Work(ipLF+jOff*JNUM),
     &                                        nBasSh(lSym,iaSh)*JNUM,
     &                                     Work(ipMO+ioffShb),1,
     &                                    ONE,Work(ipLab(iaSh)),1)


                                 EndIf

                               EndIf

                            End Do

c --- The following re-assignement is used later on to check if the
c --- iaSh vector LaJ[k] can be neglected because identically zero


                            iWork(ip_Lab+iaSh-1) = ipLab(iaSh)*
     &                                             Min(1,ibcount)
     &                                  + ipAbs*(1-Min(1,ibcount))

                         End Do


                     Else   ! lSym < kSym


                         Do iSh=1,iWork(ipIndSh)

                            iaSh = iWork(ipIndSh+iSh)

                            iOffSha = kOffSh(iaSh,lSym)

                            iWork(ip_Lab+iaSh-1) = ipChoT + iOffSha*JNUM

                            ibcount=0

                            Do ibSh=1,nShell

                               iOffShb = kOffSh(ibSh,kSym)

                               iShp = iTri(iaSh,ibSh)

                               If ( iShp_rs(iShp) .gt. 0 ) Then

                                 If (nnBstRSh(JSym,iShp_rs(iShp),iLoc)*
     &                             nBasSh(lSym,iaSh)*
     &                             nBasSh(kSym,ibSh) .gt. 0
     &                             .and. sqrt(abs(Work(ipSk+ibSh-1)*
     &                                        SvShp(iShp_rs(iShp)) ))
     &                             .ge. thrv(jDen) ) Then

                                   ibcount = ibcount + 1

                                   jOff = iOffShp(kSym,iShp_rs(iShp)) -
     &                                    nBasSh(kSym,iaSh)*
     &                                    nBasSh(lSym,ibSh)*
     &                             Min(0,(ibSh-iaSh))/Max(1,(iaSh-ibSh))


C ---  LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
C ---------------------------------------

                                 CALL DGEMV_('T',nBasSh(kSym,ibSh),
     &                                       JNUM*nBasSh(lSym,iaSh),
     &                                    ONE,Work(ipLF+jOff*JNUM),
     &                                        nBasSh(kSym,ibSh),
     &                                     Work(ipMO+ioffShb),1,
     &                                    ONE,Work(ipLab(iaSh)),1)

                                 EndIf

                               Endif

                            End Do

c --- The following re-assignement is used later on to check if the
c --- iaSh vector LaJ[k] can be neglected because identically zero

                            iWork(ip_Lab+iaSh-1) = ipLab(iaSh)*
     &                                             Min(1,ibcount)
     &                                  + ipAbs*(1-Min(1,ibcount))

                         End Do

                     EndIf

                     CALL CWTIME(TCT2,TWT2)
                     tmotr(1) = tmotr(1) + (TCT2 - TCT1)
                     tmotr(2) = tmotr(2) + (TWT2 - TWT1)

C --- Prepare the J-screening

                      CALL CWTIME(TCS1,TWS1)

                      IF (lSym.ge.kSym) Then


                         Do iSh=1,iWork(ipIndSh)

                            iaSh = iWork(ipIndSh+iSh)

                            ipFaa = ipFk + MxBasSh + iaSh - 1

                            iaSkip= Min(1,Max(0,
     &                              abs(ipLab(iaSh)-ipAbs))) ! = 1 or 0

                            Do i=1,iaSkip

C ---  Faa,[k] = sum_J  (LaJ[k])**2
C ----------------------------------
                               Do jv=1,JNUM

                                  xfjv = dble(min(1,jv-1))

                                  Do ia=1,nBasSh(lSym,iaSh)

                                     ipFia = ipFk + ia - 1

                                     ipLaa = ipLab(iaSh)
     &                                     + nBasSh(lSym,iaSh)*(jv-1)
     &                                     + ia - 1

                                     Work(ipFia) = xfjv*Work(ipFia)
     &                                           + Work(ipLaa)**2
                                  End Do
                               End Do

                               CALL FindMax(ipFk,'N',
     &                                      nBasSh(lSym,iaSh),
     &                                      1,ipFaa)

                            End Do

                         End Do


                      Else   ! lSym < kSym


                         Do iSh=1,iWork(ipIndSh)

                            iaSh = iWork(ipIndSh+iSh)

                            ipFaa = ipFk + MxBasSh + iaSh - 1

                            iaSkip= Min(1,Max(0,
     &                              abs(ipLab(iaSh)-ipAbs))) ! = 1 or 0

                            Do i=1,iaSkip

C ---  Faa,[k] = sum_J  (LJa[k])**2
C ----------------------------------
                               Do ia=1,nBasSh(lSym,iaSh)

                                  ipFia = ipFk + ia - 1

                                  Do jv=1,JNUM

                                     xfjv = dble(min(1,jv-1))

                                     ipLaa = ipLab(iaSh)
     &                                     + JNUM*(ia-1)
     &                                     + jv - 1

                                     Work(ipFia) = xfjv*Work(ipFia)
     &                                           + Work(ipLaa)**2
                                  End Do
                               End Do

                               CALL FindMax(ipFk,'N',
     &                                      nBasSh(lSym,iaSh),
     &                                      1,ipFaa)

                            End Do

                         End Do


                      EndIf

                      CALL CWTIME(TCS2,TWS2)
                      tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                      tscrn(2) = tscrn(2) + (TWS2 - TWS1)


C------------------------------------------------------------
C --- Compute exchange matrix for the interacting shell pairs
C------------------------------------------------------------
                     CALL CWTIME(TCX1,TWX1)

                     IF (lSym.ge.kSym) Then

                         Do lSh=1,iWork(ipIndSh)

                            iaSh = iWork(ipIndSh+lSh)

                            ipFaa = ipFk + MxBasSh + iaSh - 1

                            iaSkip= Min(1,Max(0,
     &                              abs(ipLab(iaSh)-ipAbs))) ! = 1 or 0

                            iOffSha = kOffSh(iaSh,lSym)

                            mSh = 1

                            Do while (mSh.le.iWork(ipIndSh))

                               ibSh = iWork(ipIndSh+mSh)

                               ipFbb = ipFk + MxBasSh + ibSh - 1

                               ibSkip = Min(1,Max(0,
     &                                        abs(ipLab(ibSh)-ipAbs)))

                               iShp = iTri(iaSh,ibSh)

                               iOffShb = kOffSh(ibSh,lSym)

                               iOffAB = nnBfShp(iShp,lSym)

                               ipKI = ipKLT(jDen) + ISTLT(lSym) + iOffAB

                               xFab = sqrt(abs(Work(ipFaa)*Work(ipFbb)))

                               If (Work(ipMLk+lSh-1)*
     &                             Work(ipMLk+mSh-1) .lt. tau(jDen))Then


                                    mSh = iWork(ipIndSh) !skip rest


                               ElseIf(iaSh.eq.ibSh
     &                                .and.xFab.ge.tau(jDen)/MaxRedT
     &                                .and.iaSkip.eq.1)Then

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(a,J)[k] * L(b,J)[k]
C -------------------------------------------------------------------
                               nBs = Max(1,nBasSh(lSym,iaSh))

                               CALL DGEMM_Tri('N','T',nBasSh(lSym,iaSh),
     &                                           nBasSh(lSym,ibSh),JNUM,
     &                                  FActX(jDen),Work(ipLab(iaSh)),
     &                                           nBs,
     &                                           Work(ipLab(iaSh)),
     &                                           nBs,
     &                                       ONE,Work(ipKI),
     &                                               nBs)

c                         write(6,*)'Symm= ',lsym,' Sh-pair: ',iaSh,ibSh
c                               CALL TRIPRT('F',' ',Work(ipKI),nBs)

                               ElseIf (iaSh.gt.ibSh
     &                                .and.xFab.ge.tau(jDen)/MaxRedT
     &                                 .and. iaSkip*ibSkip.eq.1) Then

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(a,J)[k] * L(b,J)[k]
C -------------------------------------------------------------------
                                  nBsa = Max(1,nBasSh(lSym,iaSh))
                                  nBsb = Max(1,nBasSh(lSym,ibSh))

                                  CALL DGEMM_('N','T',nBasSh(lSym,iaSh),
     &                                           nBasSh(lSym,ibSh),JNUM,
     &                                   FActX(jDen),Work(ipLab(iaSh)),
     &                                           nBsa,
     &                                           Work(ipLab(ibSh)),
     &                                           nBsb,
     &                                       ONE,Work(ipKI),
     &                                               nBsa)

c                         write(6,*)'Symm= ',lsym,' Sh-pair: ',iaSh,ibSh
c                               CALL TRIPRT('FI',' ',Work(ipKI),nBsa)

                               EndIf


                               mSh = mSh + 1  ! update shell counter


                            End Do

                         End Do


                     ELSE   ! lSym < kSym


                         Do lSh=1,iWork(ipIndSh)

                            iaSh = iWork(ipIndSh+lSh)

                            ipFaa = ipFk + MxBasSh + iaSh - 1

                            iaSkip= Min(1,Max(0,
     &                              abs(ipLab(iaSh)-ipAbs))) ! = 1 or 0

                            iOffSha = kOffSh(iaSh,lSym)

                            mSh = 1

                            Do while (mSh.le.iWork(ipIndSh))

                               ibSh = iWork(ipIndSh+mSh)

                               ipFbb = ipFk + MxBasSh + ibSh - 1

                               ibSkip = Min(1,Max(0,
     &                                        abs(ipLab(ibSh)-ipAbs)))

                               iShp = iTri(iaSh,ibSh)

                               iOffShb = kOffSh(ibSh,lSym)

                               iOffAB = nnBfShp(iShp,lSym)

                               ipKI = ipKLT(jDen) + ISTLT(lSym) + iOffAB

                               xFab = sqrt(abs(Work(ipFaa)*Work(ipFbb)))

                               If (Work(ipMLk+lSh-1)*
     &                             Work(ipMLk+mSh-1) .lt. tau(jDen))Then


                                    mSh = iWork(ipIndSh)  ! skip rest


                               ElseIf(iaSh.eq.ibSh
     &                                .and.xFab.ge.tau(jDen)/MaxRedT
     &                                .and.iaSkip.eq.1)Then

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(J,a)[k] * L(J,b)[k]
C -------------------------------------------------------------------

                               nBs = Max(1,nBasSh(lSym,iaSh))

                               CALL DGEMM_Tri('T','N',nBasSh(lSym,iaSh),
     &                                           nBasSh(lSym,ibSh),JNUM,
     &                                    FActX(jDen),Work(ipLab(iaSh)),
     &                                                JNUM,
     &                                              Work(ipLab(iaSh)),
     &                                                JNUM,
     &                                            ONE,Work(ipKI),
     &                                                nBs)


c                         write(6,*)'Symm= ',lsym,' Sh-pair: ',iaSh,ibSh
c                               CALL TRIPRT('F',' ',Work(ipKI),nBs)

                               ElseIf (iaSh.gt.ibSh
     &                                .and.xFab.ge.tau(jDen)/MaxRedT
     &                                 .and. iaSkip*ibSkip.eq.1) Then

C ---  F(a,b)[k] = F(a,b)[k] + FactX * sum_J  L(J,a)[k] * L(J,b)[k]
C -------------------------------------------------------------------

                                  nBs = Max(1,nBasSh(lSym,iaSh))

                                  CALL DGEMM_('T','N',nBasSh(lSym,iaSh),
     &                                           nBasSh(lSym,ibSh),JNUM,
     &                                    FactX(jDen),Work(ipLab(iaSh)),
     &                                           JNUM,
     &                                           Work(ipLab(ibSh)),
     &                                           JNUM,
     &                                       ONE,Work(ipKI),
     &                                               nBs)

c                         write(6,*)'Symm= ',lsym,' Sh-pair: ',iaSh,ibSh
c                               CALL TRIPRT('F',' ',Work(ipKI),nBs)

                               EndIf


                               mSh = mSh + 1  ! update shell counter


                            End Do

                         End Do


                     ENDIF

                     CALL CWTIME(TCX2,TWX2)
                     texch(1) = texch(1) + (TCX2 - TCX1)
                     texch(2) = texch(2) + (TWX2 - TWX1)


                    End Do  ! loop over k MOs

                 End Do   ! loop over MOs symmetry

               End Do   ! loop over densities


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

                        ipL = ipLrs + nRS*(jvc-1)
                        Work(ipjDiag+jrs-1) = Work(ipjDiag+jrs-1)
     &                                      + Work(ipL+krs-1)**2
                     End Do

                   End Do

                  Else

                   Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        ipL = ipLrs + nRS*(jvc-1)
                        Work(ipDiag+jrs-1) = Work(ipDiag+jrs-1)
     &                                     - Work(ipL+krs-1)**2
                     End Do

                   End Do

                  EndIf

#else
                  Do krs=1,nRS

                     mrs = iiBstR(JSYM,iLoc) + krs
                     jrs = IndRed(mrs,iLoc) ! address in 1st red set

                     Do jvc=1,JNUM

                        ipL = ipLrs + nRS*(jvc-1)
                        Work(ipDiag+jrs-1) = Work(ipDiag+jrs-1)
     &                                     - Work(ipL+krs-1)**2
                     End Do

                  End Do
#endif

                  CALL CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1) + (TCS2 - TCS1)
                  tscrn(2) = tscrn(2) + (TWS2 - TWS1)

               EndIf

C --------------------------------------------------------------------
C --- First half Active transformation  Lvb,J = sum_a  C(v,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCINT1,TWINT1)

C --- Set up the skipping flags and the pointers ipLpq
C --- The memory used before for the full-dimension AO-vectors
C ---     is now re-used to store half and full transformed
C ---     vectors in the active space
C -------------------------------------------------------------
               lChoa=0
               Do i=1,nSym

                  k = Muld2h(i,JSYM)
                  iSkip(k) = Min(1,
     &                 nAorb(k)*nBas(i))

                  ipLpq(k) = ipLF + lChoa       ! Lvb,J
                  ipLxy(k) = ipLF + lChoa       ! Lvw,J
     &                     + nAorb(k)*nBas(i)*JNUM

                  lChoa = lChoa + nAorb(k)*nBas(i)*JNUM
                  If (k.le.i) lChoa = lChoa + nnA(k,i)*JNUM

               End Do

               iSwap = 0  ! Lvb,J are returned
               kMOs = 1  !
               nMOs = 1  ! Active MOs (1st set)

               CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                          JSYM,iSwap,IREDC,nMOs,kMOs,ipAorb,nAorb,
     &                          ipLpq,iSkip,DoRead)


               if (irc.ne.0) then
                  rc = irc
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

                       ipLvb = ipLpq(iSymb) + NAv*NBAS(iSymb)*(JVC-1)
                       ipLvw = ipLxy(iSymb) + nnA(iSymb,iSymb)*(JVC-1)

                       CALL DGEMM_Tri('N','T',NAv,NAv,NBAS(iSymb),
     &                            One,Work(ipLvb),NAv,
     &                                Work(ipAorb(iSymb,1)),NAv,
     &                           Zero,Work(ipLvw),NAv)

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

                       ipLvb = ipLpq(iSymv) + NAv*NBAS(iSymb)*(JVC-1)
                       ipLvw = ipLxy(iSymv) + nnA(iSymv,iSymb)*(JVC-1)

                       CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),
     &                            One,Work(ipLvb),NAv,
     &                                Work(ipAorb(iSymb,1)),NAw,
     &                           Zero,Work(ipLvw),NAv)

                      End Do

                     EndIf

                  End Do

               EndIf
C
C
C *************** EVALUATION OF THE (WA|XY) INTEGRALS ***********

               DoTraInt = JRED.eq.myJRED2.and.iBatch.eq.nBatch

               CALL CHO_eval_waxy(irc,ipScr,ipLpq,ipLxy,ipInt,nAorb,
     &                            JSYM,JNUM,DoTraInt)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C --------------------------------------------------------------------
C --------------------------------------------------------------------


            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               add  = .true.
               Call play_casscf_sto(irc,iLoc,nDen,JSYM,ISTLT,ISSQ,
     &                                  ipFLT,ipFab,mode,add)
            EndIf

C --- free memory
            CALL GETMEM('FullV','Free','Real',ipLF,LFMAX*nVec)
            Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

            If(JSYM.eq.1)Then
              do jden=nden,1,-1
               Call GetMem('rsFC','Free','Real',ipFab(jden),nRS)
               Call GetMem('rsD','Free','Real',ipDab(jden),nRS)
              end do
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
               Call GaDsum(Work(ipjDiag),nnBSTR(JSYM,1))
               Call Daxpy_(nnBSTR(JSYM,1),xone,Work(ipjDiag),1,
     &                    Work(ipDiag+iiBstR(JSYM,1)),1)
               Call Fzero(Work(ipjDiag),nnBSTR(JSYM,1))
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


         Call GetMem('Mtmp','Free','REAL',ipItmp,Mwaxy)

1000  CONTINUE


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


      CALL GETMEM('F(k)ss','Free','Real',ipFk,MxBasSh+nShell)
      Call GetMem('ip_SvShp','Free','Real',ip_SvShp,2*nnShl)
      Call GetMem('ip_iShp_rs','Free','Inte',ip_iShp_rs,nnShl_tot)
      Call GetMem('ip_nnBfShp','Free','Inte',ip_nnBfShp,nnShl_tot*nSym)
      Call GetMem('ip_kOffSh','Free','Inte',ip_kOffSh,nShell*nSym)
      Call GetMem('ip_Lab','Free','Inte',ip_Lab,nShell)
      Call GetMem('Indx','Free','Inte',ipIndx,(nShell+1)*nnO)
      Call GetMem('SKsh','Free','Real',ipSKsh,nShell*nnO)
      Call GetMem('MLk','Free','Real',ipML,nShell*nnO)
      Call GetMem('yc','Free','Real',ipY,MaxB*nnO)
      Call GetMem('absc','Free','Real',ipAbs,MaxB)
      CALL GETMEM('diahI','Free','Real',ipDIAH(1),NNBSQ)
#if defined (_MOLCAS_MPP_)
      If (nProcs.gt.1 .and. Update .and. Is_Real_Par())
     &    CALL GETMEM('diagJ','Free','Real',ipjDIAG,NNBSTMX)
#endif
      CALL GETMEM('diagI','Free','Real',ipDIAG,NNBSTRT(1))


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
#ifdef _DEBUG_

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


      rc  = 0


      Return
      END

**************************************************************

      SUBROUTINE play_casscf_sto(irc,iLoc,nDen,JSYM,ISLT,ISSQ,
     &                        ipXLT,ipXab,mode,add)

      Implicit Real*8 (a-h,o-z)
      Integer  ISLT(8),ISSQ(8,8),cho_isao,nDen
      External cho_isao
      Integer ipXLT(nDen),ipXab(nDen)
      Logical add
      Character*6 mode

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)
******
      iRS2F(i,j)  = iWork(ip_iRS2F-1+2*(j-1)+i)
************************************************************************



      xf=0.0d0
      if (add) xf=1.0d0 !accumulate contributions

      If (mode.eq.'toreds'.and.JSYM.eq.1) then ! TOTAL SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            Do jDen=1,nDen

               kfrom = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(ipXab(jDen)+jRab-1) = xf*Work(ipXab(jDen)+jRab-1)
     &                                  +    Work(kfrom)

            End Do

         End Do  ! jRab loop


      ElseIf (mode.eq.'tofull'.and.JSYM.eq.1) then  ! TOTAL SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            Do jDen=1,nDen

               kto = ipXLT(jDen) + isLT(iSyma) + iab - 1

               Work(kto) = xf*Work(kto)
     &                   +    Work(ipXab(jDen)+jRab-1)

            End Do

         End Do  ! jRab loop


      ElseIf (mode.eq.'tosqrt'.and.JSYM.ne.1) then
c      ! NON TOTAL-SYMMETRIC

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set

            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = MulD2h(jSym,iSyma) ! sym(a) .gt. sym(b)

            ias   = iag - ibas(iSyma)
            ibs   = ibg - ibas(iSymb)

            iab   = nBas(iSyma)*(ibs-1) + ias

            Do jDen=1,nDen

               kto = ipXLT(jDen) - 1 + isSQ(iSyma,iSymb) + iab

               Work(kto) = sqrt(abs(Work(ipXab(jDen)+kRab-1)))

            End Do

         End Do  ! jRab loop


      ElseIf (mode.eq.'tosqrt'.and.JSYM.eq.1) then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab ! already in 1st red set

            iag   = iRS2F(1,kRab)  !global address
            ibg   = iRS2F(2,kRab)

            iSyma = cho_isao(iag)  ! sym(a)=sym(b)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)

            iab   = nBas(iSyma)*(ibs-1) + ias
            iba   = nBas(iSyma)*(ias-1) + ibs

            Do jDen=1,nDen

               kto = ipXLT(jDen) - 1 + isSQ(iSyma,iSyma)

               Work(kto+iab) = sqrt(abs(Work(ipXab(jDen)+kRab-1)))

               Work(kto+iba) = sqrt(abs(Work(ipXab(jDen)+kRab-1)))

            End Do

         End Do  ! jRab loop


      Else

         write(6,*)'Wrong input parameters. JSYM,mode = ',JSYM,mode
         irc = 66
         Call Qtrace()
         Call abend()

      EndIf

      irc = 0

      Return
      End
