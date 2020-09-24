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
      SUBROUTINE CHO_FOCK_RASSI(ipDLT,ipMO1,ipMO2,ipFLT,ipInt)

**********************************************************************
*  Author : F. Aquilante
*
*  Note: this routine should be used only to compute the FI matrix
*        for cases in which the 2 sets of MOs are identical and
*        therefore FI is symmetric and stores as LT matrix.
*        For historical reasons, the routine is written as if it
*        could handle the more general case of different MOs.
*
*        Please, for that purpose, use CHO_FOCK_RASSI_X instead!
*
*
C *************** INACTIVE AO-BASIS FOCK MATRIX **********************
C
C   FI(ab) = 2 * sum_J  Lab,J * U(J)  -  sum_Jk  Yka,J * Xkb,J
C
C      U(J) = sum_gd  Lgd,J * DI(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index   belonging to (Inactive)
C      v,w,x,y:  MO-indeces belonging to (Active)
C
**********************************************************************

      Implicit Real*8 (a-h,o-z)

      Integer   rc,ipLxy(8),ipScr(8,8)
      Integer   ipLab(8,2),ipOrb(8,2),nOrb(8,2)
      Integer   iSkip(8)
      Integer   ISTLT(8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2)
      Integer   ipAorb(8,2)
      Logical   Debug,timings,DoRead,DoReord
      Character*50 CFmt
      Character*14 SECNAM
      Parameter (SECNAM = 'CHO_FOCK_RASSI')
      COMMON    /CHOTIME /timings

      parameter (DoRead = .false. )
      parameter (FactCI = 1.0D0, FactXI = -1.0D0)
      parameter (zero = 0.0D0, one = 1.0D0, two = 2.0D0)
      Character*6 mode
      Logical Fake_CMO2
      COMMON / CHO_JOBS / Fake_CMO2

#include "rassi.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      parameter ( N2 = InfVec_N2 )

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
**************************************************


#ifdef _DEBUGPRINT_
c      Debug=.true.
      Debug=.false.! to avoid double printing in CASSCF-debug
#else
      Debug=.false.
#endif

      Call QEnter(SECNAM)

      DoReord = .false.
      IREDC = -1  ! unknown reduced set in core

      nDen=2
      If (Fake_CMO2) nDen = 1  ! MO1 = MO2
      kDen=nDen

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero  !time read/transform vectors
         tcoul(i) = zero  !time for computing Coulomb
         texch(i) = zero  !time for computing Exchange
         tintg(i) = zero  !time for computing (tw|xy) integrals
      end do

C ==================================================================

c --- Various offsets
c --------------------
        ISTLT(1)=0
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive F matrix
      END DO


      ipOrb(1,1) = ipMO1
      ipOrb(1,2) = ipMO2

      DO jDen=1,nDen

         nOrb(1,jDen)  = nIsh(1)
         ipAorb(1,jDen)= ipOrb(1,jDen)
     &                 + nOrb(1,jDen)*NBAS(1)

         DO ISYM=2,NSYM

            ipOrb(iSym,jDen) = ipAorb(iSym-1,jDen)
     &                       + nAsh(iSym-1)*NBAS(iSym-1)

            nOrb(iSym,jDen)  = nIsh(iSym)

            ipAorb(iSym,jDen)= ipOrb(iSym,jDen)
     &                       + nOrb(iSym,jDen)*NBAS(iSym)

         END DO

      END DO

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

        If (NumCho(jSym).lt.1) GOTO 1000

C *** memory for the (tw|xy) integrals --- temporary array
      Mtwxy = 0
      Do iSymy=1,nSym
         iSymx=MulD2h(iSymy,JSYM)
         Do iSymw=iSymy,nSym    ! iSymw.ge.iSymy (particle symmetry)
            iSymt=MulD2h(isymw,JSYM)
            Mtwxy=Mtwxy+nAsh(iSymt)*nAsh(iSymw)*nAsh(iSymx)*nAsh(iSymy)
         End Do
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
         Do iSymw=iSymy,nSym   ! iSymw.ge.iSymy (particle symmetry)
            iSymt=MulD2h(isymw,JSYM)
            ipScr(iSymw,iSymy) = kScr
            kScr=kScr+nAsh(iSymt)*nAsh(iSymw)*nAsh(iSymx)*nAsh(iSymy)
         End Do
      End Do


      iLoc = 3 ! use scratch location in reduced index arrays

C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec = 0  ! mem for storing the half-transformed vec
         mTTvec= 0  ! mem for Lvb,J and Lvw,J
         do l=1,nSym
            k=Muld2h(l,JSYM)
            mTvec = mTvec + nDen*nBas(l)*nIsh(k)
            mTTvec = mTTvec + (nBas(l)+nAsh(l))*nAsh(k)
         end do

         mTvec=Max(mTvec,mTTvec,1)

C ------------------------------------------------------------------
C ------------------------------------------------------------------

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
               call qtrace()
               call abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                        '   rc= ',irc
              call qtrace()
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then

               Call GetMem('rsDtot','Allo','Real',ipDab,nRS)
               Call GetMem('rsFC','Allo','Real',ipFab,nRS)
               Call Fzero(Work(ipDab),nRS)
               Call Fzero(Work(ipFab),nRS)
            EndIf

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec  = Min(LWORK/(nRS+mTvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec
               WRITE(6,*) 'jsym= ',jsym
               rc = 33
               CALL QTrace()
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
            Call GetMem('ChoT','Allo','Real',ipChoT,mTvec*nVec)

            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               Call swap_sto(irc,iLoc,ipDLT,ISTLT,ipDab,mode)
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
                  write(6,*)'return code = ',rc
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
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
C==========================================================
C
                  CALL CWTIME(TCC1,TWC1)

                  ipVJ = ipChoT

                  CALL DGEMV_('T',nRS,JNUM,
     &                 ONE,Work(ipLrs),nRS,
     &                 Work(ipDab),1,ZERO,Work(ipVJ),1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  Fact = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                 FactCI,Work(ipLrs),nRS,
     &                 Work(ipVJ),1,Fact,Work(ipFab),1)


                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)

               EndIf  ! Coulomb contribution


C *************** EXCHANGE CONTRIBUTIONS  ***********************

C --- Set pointers to the half-transformed Cholesky vectors
               lChoT=0
               Do iSymb=1,nSym

                  iSymp = MulD2h(JSYM,iSymb)

                  Do jDen=1,nDen

                     ipLab(iSymp,jDen) = ipChoT + lChoT  ! LpJ,b

                     lChoT = lChoT + nIsh(iSymp)*nBas(iSymb)*JNUM

                  End Do

               End Do


               iSwap = 2  ! LpJ,b are returned

               CALL CWTIME(TCR3,TWR3)

               kMOs = 1
               nMOs = nDen

C --- Set up the skipping flags
C -------------------------------------------------------------
               Do i=1,nSym

                  k = Muld2h(i,JSYM)
                  iSkip(k) = Min(1,NBAS(i)*nIsh(k))

               End Do
C -------------------------------------------------------------


C *********************** HALF-TRANSFORMATION  ****************

               CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,ipOrb,nOrb,
     &                            ipLab,iSkip,DoRead)


               CALL CWTIME(TCR4,TWR4)
               tread(1) = tread(1) + (TCR4 - TCR3)
               tread(2) = tread(2) + (TWR4 - TWR3)

               if (irc.ne.0) then
                  rc = irc
                  write(6,*)'CHO_X_getVtra failed! '
                  RETURN
               endif

               CALL CWTIME(TCX1,TWX1)

               Do iSyma=1,nSym

                  iSymk = MulD2h(JSYM,iSyma)

C ---------------------------------------------------------------------
c *** Compute the InActive exchange matrix
C
C     FI(ab) = FI(ab) + FactXI * sum_Jk  X(2)kJ,a * X(1)kJ,b
C ---------------------------------------------------------------------
                  NK = nIsh(iSymk)

                  If (iSkip(iSymk).ne.0) Then

                     ISFI = ipFLT + ISTLT(iSyma)

                     CALL DGEMM_Tri('T','N',NBAS(iSyma),NBAS(iSyma),
     &                         NK*JNUM,FactXI,Work(ipLab(iSymk,kDen)),
     &                         NK*JNUM,Work(ipLab(iSymk,1)),NK*JNUM,
     &                             One,Work(ISFI),NBAS(iSyma))



                  EndIf

C --------------------------------------------------------------------
               End Do  !loop over MOs symmetries

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)


C ************  END EXCHANGE CONTRIBUTION  ****************


C --------------------------------------------------------------------
C --- First half Active transformation  Lvb,J = sum_a  C1(v,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCR7,TWR7)

C --- Set pointers to the half-transformed Cholesky vectors
               lChoa=0
               Do iSymb=1,nSym

                  iSymp = MulD2h(JSYM,iSymb)

                  ipLab(iSymp,1) = ipChoT + lChoa  ! LpJ,b
                  ipLab(iSymp,2) = ipLab(iSymp,1)
     &                           + nAsh(iSymp)*nBas(iSymb)*JNUM

                  lChoa = lChoa
     &                  + nAsh(iSymp)*(nBas(iSymb)+nAsh(iSymb))*JNUM

               End Do

C --- Set up the skipping flags
C -------------------------------------------------------------
               Do i=1,nSym

                  k = Muld2h(i,JSYM)
                  iSkip(k) = Min(1,
     &                 NBAS(i)*nAsh(k))

               End Do

               iSwap = 0  ! Lvb,J are returned
               kMOs = 1  !
               nMOs = 1  ! Active MOs (1st set)

               CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                           JSYM,iSwap,IREDC,nMOs,kMOs,ipAorb,nAsh,
     &                           ipLab,iSkip,DoRead)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif


C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
C --------------------------------------------------------------------
                  Do iSymb=1,nSym

                     iSymv = MulD2h(JSYM,iSymb)
                     NAv = nAsh(iSymv)
                     NAw = nAsh(iSymb) ! iSymb=iSymw

                     If(NAv*NAw.ne.0)Then

                      Do JVC=1,JNUM

                       ipLvb = ipLab(iSymv,1) + NAv*NBAS(iSymb)*(JVC-1)
                       ipLvw = ipLab(iSymv,2)
     &                       + nAsh(iSymv)*nAsh(iSymb)*(JVC-1)

                       CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),
     &                            One,Work(ipLvb),NAv,
     &                                Work(ipAorb(iSymb,kDen)),NAw,
     &                           Zero,Work(ipLvw),NAv)

                      End Do

                     EndIf

                  End Do


C
C
C *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********

               CALL CWTIME(TCINT1,TWINT1)

               DoReord = JRED.eq.JRED2.and.iBatch.eq.nBatch

               Do iSym=1,nSym
                  ksym = mulD2h(iSym,JSYM)
                  ipLxy(kSym) = ipLab(iSym,2) ! switch to column storage
               End Do

               CALL CHO_rassi_twxy(irc,ipScr,ipLxy,ipInt,nAsh,
     &                                 JSYM,JNUM,DoReord)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C ---------------- END (TW|XY) EVALUATION -----------------------


            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               Call swap_sto(irc,iLoc,ipFLT,ISTLT,ipFab,mode)
            EndIf

C --- free memory
            Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

            If(JSYM.eq.1)Then
              Call GetMem('rsFC','Free','Real',ipFab,nRS)
              Call GetMem('rsDtot','Free','Real',ipDab,nRS)
            EndIf


999         Continue

         END DO   ! loop over red sets

         Call GetMem('Mtmp','Free','REAL',ipItmp,Mtwxy)

1000  CONTINUE


      END DO  ! loop over JSYM


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1


*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky RASSI timing from '//SECNAM
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
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


c Print the Fock-matrix
#ifdef _DEBUGPRINT_

      if(Debug) then !to avoid double printing in RASSI-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      DO ISYM=1,NSYM
        ISLFI=ipFLT+ISTLT(ISYM)
        IF( NBAS(ISYM).GT.0 ) THEN
          WRITE(6,'(6X,A)')'***** INACTIVE FOCK MATRIX ***** '
          WRITE(6,'(6X,A)')
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call TRIPRT('','',Work(ISLFI),NBAS(ISYM))
        ENDIF
      END DO

      endif

#endif

      rc  = 0

      CAll QExit(SECNAM)

      Return
      END

**************************************************************

      SUBROUTINE swap_sto(irc,iLoc,ipXLT,ISLT,ipXab,mode)

      Implicit Real*8 (a-h,o-z)
      Integer  ISLT(8),cho_isao
      External cho_isao
      Integer ipXLT,ipXab
      Character*6 mode

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)
******
      iRS2F(i,j)  = iWork(ip_iRS2F-1+2*(j-1)+i)
************************************************************************


      jSym = 1 ! only total symmetric density

      If (mode.eq.'toreds') then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            kfrom = ipXLT + isLT(iSyma) + iab - 1

            Work(ipXab+jRab-1) = Work(kfrom)

         End Do  ! jRab loop

      ElseIf (mode.eq.'tofull') then

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            iab   = iTri(ias,ibs)

            kto = ipXLT + isLT(iSyma) + iab - 1

            Work(kto) = Work(kto)
     &                + Work(ipXab+jRab-1)


         End Do  ! jRab loop

      Else

         write(6,*)'Wrong input parameter. mode = ',mode
         irc = 66
         Call Qtrace()
         Call abend()

      EndIf

      irc = 0

      Return
      End
