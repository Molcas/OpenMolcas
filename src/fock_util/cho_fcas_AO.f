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

      SUBROUTINE CHO_FCAS_AO(rc,ipFA,ipFI,ipQmat,nForb,nIorb,nAorb,
     &                          FactXI,ipPorb,ipDI,ipDA1,ipDA2,DoActive,
     &                          DoQmat,ipChM,nChM,ipInt,ExFac)

**********************************************************************
*  Author : F. Aquilante
*
C *************** INACTIVE AO-BASIS FOCK MATRIX **********************
C
C   FI(ab) = 2 * sum_J  Lab,J * U(J)  -  sum_Jk  Lka,J * Lkb,J
C
C ***************   ACTIVE AO-BASIS FOCK MATRIX **********************
C
C   FA(ab) = sum_J  Lab,J * V(J)  -  0.5 * sum_Jw  Lwa,J * Lwb,J
C
C ***************   AUXILIARY  Q-MATRIX   ****************************
C
C   Q(av) = 2 * sum_Jw  L(w,aJ) * sum_xy P(vw,xy) * L(xy,J)
C
**********************************************************************
C
C      U(J) = sum_gd  Lgd,J * DI(gd)
C      V(J) = sum_gd  Lgd,J * DA(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index   belonging to (Frozen+Inactive)
C      v,w,x,y:  MO-indeces belonging to (Active)
C
**********************************************************************

      Implicit Real*8 (a-h,o-z)

      Integer   rc,ipLab(8,3),ipLxy(8),ipScr(8,8)
      Integer   ipOrb(8,3),nOrb(8,3)
      Integer   ISTAQ(8),ISTAV(8),iSkip(8)
      Integer   ISTLT(8),ISTCH(8),ISZW(8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2),tqmat(2)
      Real*8    ExFac
*      Integer   ipDA1,ipDA2,ipDI
      Integer   ipDA1,ipDA2(8,8,8),ipDI
      Integer   ipPorb,ipFA,ipFI
      Integer   ipDLT(2),ipFLT(2),ipDab(2),ipFab(2)
      Integer   nForb(8),nIorb(8),nAorb(8),nPorb(8),nnA(8,8),nChM(8)
      Logical   Debug,timings,DoRead,DoReord,DoActive,DoQmat
      Character*50 CFmt
      Character*11 SECNAM
      Parameter (SECNAM = 'CHO_FCAS_AO')
      COMMON    /CHOTIME /timings

      parameter (FactCI = 1.0D0)
      parameter (FactCA = 1.0D0, FactXA = -0.5D0)
      parameter (zero = 0.0D0, one = 1.0D0, two = 2.0d0)

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      parameter ( N2 = InfVec_N2 )
      Logical add
      Character*6 mode

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
************************************************************************

#ifdef _DEBUGPRINT_
c      Debug=.true.
      Debug=.false.! to avoid double printing in CASSCF-debug
#else
      Debug=.false.
#endif


      DoRead  = .false.
      DoReord = .false.
      IREDC = -1  ! unknown reduced set in core

      ipDLT(1) = ipDI    ! some definitions
      ipDLT(2) = ipDA1
      ipFLT(1) = ipFI
      ipFLT(2) = ipFA

      if(ExFac.ne.1.0d0) then
          write(6,*) 'WARNING: if MCPDFT is used and your code is'
          write(6,*) 'passing by here errors may occur.'
          write(6,*) 'Please check with:'
          write(6,*) 'Giovanni Li Manni'
          write(6,*) 'giovannilimanni@gmail.com'
      end if
      nDen = 1
      if (DoActive) nDen=2

        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time read/transform vectors
           tcoul(i) = zero  !time for computing Coulomb
           texch(i) = zero  !time for computing Exchange
           tintg(i) = zero  !time for computing (tw|xy) integrals
           tqmat(i) = zero  !time for computing Q-matrix
        end do

c --- Define MOs in the Primary space
c --- ( Frozen + Inactive + Active )
c -----------------------------------
        do i=1,nSym
           nPorb(i) = nForb(i) + nIorb(i) + nAorb(i)
        end do

      Call set_nnA(nSym,nAorb,nnA)

C ==================================================================

c --- Various offsets
c --------------------
        ISTAQ(1)=0
        ISTAV(1)=0
        ISTLT(1)=0
        ISTCH(1)=0
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        NP=NPORB(ISYM-1)
        NP2=NB*NP
        NV2=NB*NAORB(ISYM-1)
        NCH=NB*NCHM(ISYM-1)
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive and Active D and F mat
        ISTAQ(ISYM)=ISTAQ(ISYM-1)+NP2 ! MOs coefficients
        ISTAV(ISYM)=ISTAV(ISYM-1)+NV2 ! Q-matrix
        ISTCH(ISYM)=ISTCH(ISYM-1)+NCH ! "Cholesky MOs"
      END DO

      Do iSym=1,nSym        ! MOs to feed in cho_x_getvtra

         ipOrb(iSym,1) = ipPorb + ISTAQ(iSym)
         nOrb(iSym,1)  = nForb(iSym)+nIorb(iSym)

         ipOrb(iSym,2) = ipChM + ISTCH(iSym)
         nOrb(iSym,2)  = nChM(iSym)

         ipOrb(iSym,3) = ipPorb + ISTAQ(iSym)
     &                 + nOrb(iSym,1)*nBas(iSym)
         nOrb(iSym,3)  = nAorb(iSym)

      End Do

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

         If (NumCho(jSym).lt.1) GOTO 1000

C *** memory for the (tw|xy) integrals --- temporary array
      Mtwxy = 0
      Do iSymy=1,nSym
         iSymx=MulD2h(iSymy,JSYM)
         If (iSymx.ge.iSymy) then
            Do iSymw=iSymy,nSym
               iSymt=MulD2h(isymw,JSYM)
               If (iSymt.ge.iSymw) then
                  Mtwxy = Mtwxy + nnA(iSymt,iSymw)*nnA(iSymx,iSymy)
               End If
            End Do
          End If
      End Do

      Call GetMem('Mtmp','ALLO','REAL',ipItmp,Mtwxy)
      Call Fzero(Work(ipItmp),Mtwxy)

C *** setup pointers to the symmetry blocks of (tw|xy)
      Do i=1,nSym
         Do j=1,nSym
            ipScr(j,i) = ipItmp
         End Do
      End Do

      kScr=ipItmp
      Do iSymy=1,nSym
         iSymx=MulD2h(iSymy,JSYM)
         If (iSymx.ge.iSymy) then
            Do iSymw=iSymy,nSym   ! iSymw.ge.iSymy
               iSymt=MulD2h(isymw,JSYM)
               If (iSymt.ge.iSymw) then
                  ipScr(iSymw,iSymy) = kScr
                  ipScr(iSymy,iSymw) = kScr !symmetrization
                  kScr = kScr + nnA(iSymt,iSymw)*nnA(iSymx,iSymy)
               End If
            End Do
          End If
      End Do

C --- Set up the skipping flags + some initializations --------
C -------------------------------------------------------------
         Do i=1,nSym

            k = Muld2h(i,JSYM)
            iSkip(i) = Min(1,nBas(i)*nBas(k)) ! skip Lik vector
            iSkip(i) = iSkip(i)*(nPorb(i)+nChM(i))

            ipLab(i,1) = -6666  ! pointers to Lk,Jb
            ipLab(i,2) = -6666  ! pointers to "cholesky MOs vectors"
            ipLab(i,3) = -6666  ! pointers to Lvb,J
            ipLxy(i) = -6666  ! pointers to Lxy,J

         End Do
C -------------------------------------------------------------


C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec  = 0  ! mem for storing the half-transformed vec
         mTZvec = 0

         do l=1,nSym
            k=Muld2h(l,JSYM)
            mTvec = mTvec + nBas(l)*Max((nForb(k)+nIorb(k)),
     &              nAorb(k),nChM(k)) + nnA(k,l)
            mTZvec = mTZvec + nAorb(l)*nAorb(k)
         end do

         if (.not.DoQmat) mTZvec=0
         mTvec = Max(mTvec,1)

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
               call abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM //' cho_X_setred non-zero rc (',irc,').'
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then

               Call GetMem('rsDI','Allo','Real',ipDab(1),nRS)
               Call GetMem('rsFI','Allo','Real',ipFab(1),nRS)

               if(DoActive)then
                 Call GetMem('rsDA1','Allo','Real',ipDab(2),nRS)
                 Call GetMem('rsFA','Allo','Real',ipFab(2),nRS)
               endif

            EndIf

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            mNeed1 = Max(nRS,mTZvec)

            nVec  = Min(LWORK/(mNeed1+mTvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',mNeed1+mTvec
               WRITE(6,*) 'read in ',mNeed1,' and transform to ',mTvec
               WRITE(6,*) 'of jsym= ',jsym,' and JRED= ',JRED
               rc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,mNeed1*nVec)
            Call GetMem('ChoT','Allo','Real',ipChoT,mTvec*nVec)

            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               add  = .false.
               Call change_sto(irc,iLoc,nDen,ipDLT,ipDab,mode,add)
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

               If (NUMV.le.0 .or.NUMV.ne.JNUM) then
                  rc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

               IF(JSYM.eq.1)THEN

C ************ INACTIVE COULOMB CONTRIBUTION  ****************
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
C==========================================================

                  CALL CWTIME(TCC1,TWC1)

                  ipVJ = ipChoT

                  CALL DGEMV_('T',nRS,JNUM,
     &                 ONE,Work(ipLrs),nRS,
     &                 Work(ipDab(1)),1,ZERO,Work(ipVJ),1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  xfac = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                 FactCI,Work(ipLrs),nRS,
     &                 Work(ipVJ),1,xfac,Work(ipFab(1)),1)


                  If (DoActive) Then
C ************ ACTIVE COULOMB CONTRIBUTION  ****************
C
C---- Computing the intermediate vector V(J)
C ---------------------------------------
C --- U{#J} <- U{#J}  +  sum_rs  L(rs,{#J}) * DA(rs)
C==========================================================
C
                     CALL DGEMV_('T',nRS,JNUM,
     &                          ONE,Work(ipLrs),nRS,
     &                          Work(ipDab(2)),1,ZERO,Work(ipVJ),1)

C --- FA(rs){#J} <- FA(rs){#J} + FactCA * sum_J L(rs,{#J})*U{#J}
C===============================================================

                     xfac = dble(min(jVec-iVrs,1))

                     CALL DGEMV_('N',nRS,JNUM,
     &                          FactCA,Work(ipLrs),nRS,
     &                          Work(ipVJ),1,xfac,Work(ipFab(2)),1)

                  EndIf


                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)


               ENDIf  ! Coulomb (jsym=1)

C ************ BEGIN EXCHANGE CONTRIBUTIONS  ****************

C --- Set pointers to the half-transformed Cholesky vectors
               lChoT=0
               Do iSymb=1,nSym

                  iSymp = MulD2h(JSYM,iSymb)

                  ipLab(iSymp,1) = ipChoT + lChoT  ! LkJ,b
                  ipLab(iSymp,2) = ipLab(iSymp,1)  ! LxJ,b

                  ipLab(iSymp,3) = ipLab(iSymp,1)  ! Lbv,J
                  ipLxy(iSymp) = ipLab(iSymp,3)    ! Lvw,J
     &                         + nAorb(iSymp)*nBas(iSymb)*JNUM

                  lChoT = lChoT + nBas(iSymb)*
     &                    Max((nForb(iSymp)+nIorb(iSymp)),
     &                         nAorb(iSymp),nChM(iSymp))*JNUM
     &                  + nnA(iSymp,iSymb)*JNUM

#ifdef _DEBUGPRINT_
            write(6,*)'JRED,iBatch,nBatch= ',jred,iBatch,nBatch
            write(6,*)'JSYM,iSymp,iSymb,lChot= ',JSYM,iSymp,iSymb,lChot
            write(6,*)'Lxb starts in= ',ipLab(iSymp,3)
            write(6,*)'and occupies ',nAorb(iSymp)*nBas(iSymb)*JNUM
            write(6,*)'Lxy starts in= ',ipLxy(iSymp)
            write(6,*)'and occupies ',nnA(iSymp,iSymb)*JNUM
#endif

               End Do

               iSwap = 2  ! LpJ,b are returned

C *********************** INACTIVE HALF-TRANSFORMATION  ****************

               kMOs = 1  ! inactive MOs
               nMOs = 1

               CALL CWTIME(TCR3,TWR3)

               CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                         JSYM,iSwap,IREDC,nMOs,kMOs,ipOrb,nOrb,
     &                         ipLab,iSkip,DoRead)


               CALL CWTIME(TCR4,TWR4)
               tread(1) = tread(1) + (TCR4 - TCR3)
               tread(2) = tread(2) + (TWR4 - TWR3)

#ifdef _DEBUGPRINT_
       write(6,*) 'Half-transformation in the Inactive space'
       write(6,*) 'Total allocated :     ',mTvec*nVec,' at ',ipChoT
       write(6,*) 'Mem pointers ipLab :  ',(ipLab(i,1),i=1,nSym)
       write(6,*) 'ipLxy :  ',(ipLxy(i),i=1,nSym)
       write(6,*) 'iSkip :        ',(iSkip(i),i=1,nSym)
       write(6,*) 'LREAD: ',LREAD,' allocated at ',ipLrs
       write(6,*) 'JRED :        ',JRED
       write(6,*) 'JSYM :        ',JSYM
       write(6,*) 'JNUM :        ',JNUM
#endif

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif


               CALL CWTIME(TCX1,TWX1)

               Do iSyma=1,nSym

                  iSymk = MulD2h(JSYM,iSyma)

C ---------------------------------------------------------------------
c *** Compute only the LT part of the InActive exchange matrix ********
C
C     FI(ab) = FI(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
C ---------------------------------------------------------------------
                  NK = nForb(iSymk) + nIorb(iSymk)

                  If (iSkip(iSymk)*NK.ne.0) Then

                     ISFI = ipFI + ISTLT(iSyma)

                     CALL DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),
     &                         NK*JNUM,FactXI,Work(ipLab(iSymk,1)),
     &                         NK*JNUM,Work(ipLab(iSymk,1)),NK*JNUM,
     &                         One,Work(ISFI),nBas(iSyma))


c          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMA
c          CALL TRIPRT('FI: ',' ',Work(ISFI),nBas(iSyma))

                  EndIf

C --------------------------------------------------------------------
               End Do  !loop over MOs symmetries

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)


               IF(DoActive)THEN
C *********************** "CHOLESKY" HALF-TRANSFORMATION  ****************
C --------------------------------------------------------------------
c --- Using "Cholesky MOs" obtained by cholesky decomposing DA
C --------------------------------------------------------------------

                  CALL CWTIME(TCR5,TWR5)

                  kMOs = 2  ! Cholesky MOs
                  nMOs = 2

                  CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,ipOrb,nOrb,
     &                            ipLab,iSkip,DoRead)


                  CALL CWTIME(TCR6,TWR6)
                  tread(1) = tread(1) + (TCR6 - TCR5)
                  tread(2) = tread(2) + (TWR6 - TWR5)


                  if (irc.ne.0) then
                     rc = irc
                     RETURN
                  endif

                  CALL CWTIME(TCX3,TWX3)

C ---------------------------------------------------------------------
c *** Compute only the LT part of the Active exchange matrix **********
C
C     FA(ab) = FA(ab) + FactXA * sum_wJ  LwJ,a * LwJ,b
C ---------------------------------------------------------------------
                    Do iSyma=1,nSym

                     iSymw = MulD2h(JSYM,iSyma)

                     NAch= nChM(iSymw)

                     If (iSkip(iSymw)*NAch.ne.0) Then

                        ISFA = ipFA + ISTLT(iSyma)

                        CALL DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),
     &                         NAch*JNUM,FactXA,Work(ipLab(iSymw,2)),
     &                         NAch*JNUM,Work(ipLab(iSymw,2)),NAch*JNUM,
     &                         One,Work(ISFA),nBas(iSyma))


c          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMA
c          CALL TRIPRT('FA: ',' ',Work(ISFA),nBas(iSyma))

                     EndIf

C --------------------------------------------------------------------
                    End Do  !loop over MOs symmetries

                  CALL CWTIME(TCX4,TWX4)
                  texch(1) = texch(1) + (TCX4 - TCX3)
                  texch(2) = texch(2) + (TWX4 - TWX3)

               ENDIF   ! Do Active Exchange


C ************  END EXCHANGE CONTRIBUTIONS  ****************


C --------------------------------------------------------------------
C --- First half Active transformation  Lav,J = sum_b  Lab,J * C(v,a)
C --------------------------------------------------------------------

               CALL CWTIME(TCR7,TWR7)

               iSwap = 1  ! Lav,J are returned
               kMOs = 3  ! Active MOs
               nMOs = 3  ! Active MOs

               CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,ipOrb,nOrb,
     &                            ipLab,iSkip,DoRead)


               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C --------------------------------------------------------------------
C --- Active-Active transformation  Lwv,J = sum_a  C(w,a) * La,vJ
C --------------------------------------------------------------------
               IF (JSYM.eq.1) THEN   !  Lwv,J in LT-storage

                  Do iSyma=1,nSym

                     NAv = nAorb(iSyma)

                     If(NAv.ne.0)Then

                      NK   = nForb(iSyma) + nIorb(iSyma)
                      ISMO = ipPorb + ISTAQ(iSyma) + NK*nBas(iSyma)

                      Do JVC=1,JNUM

                       ipLvb = ipLab(iSyma,3) + nBas(iSyma)*NAv*(JVC-1)
                       ipLvw = ipLxy(iSyma) + nnA(iSyma,iSyma)*(JVC-1)

                       CALL DGEMM_Tri('N','N',NAv,NAv,nBas(iSyma),
     &                                One,Work(ISMO),NAv,
     &                                    Work(ipLvb),nBas(iSyma),
     &                               Zero,Work(ipLvw),NAv)

                      End Do

                     EndIf

                  End Do


               ELSE

C --------------------------------------------------------------------
C --- Active-Active transformation  Lw,vJ = sum_a  C(w,b) * Lb,vJ
C --------------------------------------------------------------------
                  Do iSymb=1,nSym

                     iSymv = MulD2h(JSYM,iSymb)
                     NAv = nAorb(iSymv)
                     NAw = nAorb(iSymb) ! iSymb=iSymw

                     If(NAv*NAw.ne.0.and.iSymv.gt.iSymb)Then

                      NK = nForb(iSymb) + nIorb(iSymb)
                      ISMO = ipPorb + ISTAQ(iSymb) + NK*nBas(iSymb)
                      ipLvb = ipLab(iSymv,3)
                      ipLvw = ipLxy(iSymv)

                      CALL DGEMM_('N','N',NAw,NAv*JNUM,nBas(iSymb),
     &                            One,Work(ISMO),NAw,
     &                                Work(ipLvb),nBas(iSymb),
     &                           Zero,Work(ipLvw),NAw)

                     EndIf

                  End Do


               ENDIF   ! jSym >=< 1

               CALL CWTIME(TCR8,TWR8)
               tread(1) = tread(1) + (TCR8 - TCR7)
               tread(2) = tread(2) + (TWR8 - TWR7)


C
C *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********

               CALL CWTIME(TCINT1,TWINT1)

               DoReord = JRED.eq.JRED2.and.iBatch.eq.nBatch

               CALL CHO_eval_twxy(irc,ipScr,ipLxy,ipInt,nAorb,
     &                      JSYM,JNUM,DoReord)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C ---------------- END (TW|XY) EVALUATION -----------------------

               kZvw = ipLrs ! mem used earlier for reading AO vectors

               IF (DoQmat) THEN
C ***************  Q-MATRIX SECTION *****************************

                  CALL CWTIME(TCQMT1,TWQMT1)

                  NTvw=0
                  do iSymw=1,nSym
                     iSymv=muld2h(iSymw,JSYM)
                     ISZW(iSymw)=NTvw     ! Offset to symm block of Zvw
                     NTvw = NTvw + nAorb(iSymv)*nAorb(iSymw)*JNUM
                  end do

C ---------------------------------------------------------------------
C --- P(vw,xy) : is totally symmetric stored by compound symmetry JSYM
C ---            the indeces {xy} are stored as packed (sym x.le.sym y)
C ---            the indeces [vw] are stored as squared
C======================================================================
C --- Z(vw){#J} <- Z(vw){#J} + sum_xy  P(vw,xy) * L(xy,{#J})
C======================================================================
                  CALL FZERO(Work(kZvw),NTvw)

                  Do iSymy=1,nSym

                    iSymx=MulD2h(iSymy,JSYM)

                    If (iSymx.le.iSymy .and. nnA(iSymx,iSymy).ne.0) Then

                      Do iSymw=1,nSym

                        iSymv = MulD2h(iSymw,JSYM)
                        Nvw = nAorb(iSymv)*nAorb(iSymw)

                        If (Nvw.gt.0) Then

                           ipZvw = kZvw + ISZW(iSymw)

                          CALL DGEMM_('N','N',Nvw,JNUM,nnA(iSymx,iSymy),
     &                         ONE,Work(ipDA2(iSymw,iSymy,JSYM)),Nvw,
     &                             Work(ipLxy(iSymy)),nnA(iSymx,iSymy),
     &                         ONE,Work(ipZvw),Nvw)

                        EndIf

                      End Do

                    Endif

                  End Do

C
C --------------------------------------------------------------
C --- Q(av) <- Q(av) + 2 * sum_Jw  L(a,wJ) * Z(v,wJ)
C===============================================================
C
                  DO iSymw=1,nSym

                     iSyma = MULD2H(iSymw,JSYM)
                     iSymv = iSyma

                     IF(iSkip(iSymw)*nAorb(iSymw)*nAorb(iSymv).ne.0)THEN

                        Nv   = nAorb(iSymv)
                        Nw   = nAorb(iSymw)
                        ISQ  = ipQmat + ISTAV(iSyma)
                        ipZvw= kZvw + ISZW(iSymw)

                        CALL DGEMM_('N','T',
     &                                 NBAS(ISYMA),Nv,Nw*NUMV,
     &                             TWO,Work(ipLab(iSymw,3)),nBas(iSyma),
     &                                 Work(ipZvw),Nv,
     &                             ONE,Work(ISQ),nBas(iSyma))

c        WRITE(6,'(6X,A,I2)')'Q-MATRIX symm=',ISYMA
c        CALL CHO_OUTPUT(Work(ISQ),1,nBas(iSyma),1,Nv,nBas(iSyma),Nv,1,6)

                     ENDIF

                  END DO

C --------------------------------------------------------

                  CALL CWTIME(TCQMT2,TWQMT2)
                  tqmat(1) = tqmat(1) + (TCQMT2 - TCQMT1)
                  tqmat(2) = tqmat(2) + (TWQMT2 - TWQMT1)

C ---------------- END Q-MATRIX SECTION -------------------------
               ENDIF
C
C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop


            if (JSYM.eq.1) then
c --- backtransform fock matrix in full storage
               mode = 'tofull'
               add  = .true.
               Call change_sto(irc,iLoc,nDen,ipFLT,ipFab,mode,add)
            endif

C --- free memory
            Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,mNeed1*nVec)

            if (JSYM.eq.1) then
               if(DoActive)then
                 Call GetMem('rsFA','Free','Real',ipFab(2),nRS)
                 Call GetMem('rsDA1','Free','Real',ipDab(2),nRS)
               endif
               Call GetMem('rsFI','Free','Real',ipFab(1),nRS)
               Call GetMem('rsDI','Free','Real',ipDab(1),nRS)
            endif


999         CONTINUE

         END DO   ! loop over red sets

         Call GetMem('Mtmp','Free','REAL',ipItmp,Mtwxy)

1000     CONTINUE

      END DO   !loop over JSYM


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1


*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky CASSCF timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ/TRANSFORM VECTORS           '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'EXCHANGE                         '
     &                           //'         ',texch(1),texch(2)
         Write(6,'(2x,A26,2f10.2)')'(TW|XY) INTEGRALS                '
     &                           //'         ',tintg(1),tintg(2)
         Write(6,'(2x,A26,2f10.2)')'Q-MATRIX                         '
     &                           //'         ',tqmat(1),tqmat(2)
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
      IF(DoQmat)THEN
        WRITE(6,'(6X,A)')
        WRITE(6,'(6X,A)')'****** Q-MATRIX ****** '
        DO ISYM=1,NSYM
         IF( NBAS(ISYM).GT.0 ) THEN
          Nv = nAorb(iSym)
          ISQ = ipQmat + ISTAV(iSym)
          WRITE(6,'(6X,A)')
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL CHO_OUTPUT(Work(ISQ),1,nBas(iSym),1,Nv,nBas(iSym),Nv,1,6)
         ENDIF
        END DO
      ENDIF

      endif

#endif

      rc  = 0


      Return
      END

**************************************************************
**************************************************************
