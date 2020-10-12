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
      SUBROUTINE CHO_Fmat(rc,nDen,ipFLT,nMOs,ipMOs,
     &                       ExFac,ipDLT)

********************************************************************
*  Author : F. Aquilante
*
C ******************* AO-BASIS FOCK MATRIX *************************
C
C   F(ab) = sum_J  Lab,J * V(J)  -  ExFac * sum_Jk  Lka,J * Lkb,J
C
********************************************************************
C
C      V(J) = sum_gd  Lgd,J * D(gd)
C
C      a,b,g,d:  AO-index
C      k:        MO-index such that D(g,d) = sum_k  X(g,k) * X(d,k)
C
C      Lkb,J = sum_a  Lab,J * X(a,k)
C
********************************************************************

      Implicit Real*8 (a-h,o-z)


      Integer   rc,nDen,ipDLT(nDen),ipFLT(nDen),ipMOs(nDen)
      Integer   nMOs(8,nDen)
      Real*8    ExFac(nDen)

      Integer   iSkip(8),ISTLT(8)
      Real*8    tread(2),tcoul(2),texch(2)

      Character*50 CFmt
      Character*8 SECNAM
      Parameter (SECNAM = 'Cho_Fmat')

      Logical   Debug,timings,DoRead
      COMMON    /CHOTIME /timings

      parameter (zero = 0.0D0, one = 1.0D0, xone = -1.0D0)

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
******
      ipOrb(i,j) = iWork(ip_Orb-1+8*(j-1)+i)
******
      ipLab(i,j) = iWork(ip_Lab-1+8*(j-1)+i)
******
      ipDr(i) = iWork(ipDab+i-1)
******
      ipFr(i) = iWork(ipFab+i-1)
************************************************************************

#ifdef _DEBUGPRINT_
c      Debug=.true.
      Debug=.false.! to avoid double printing in CASSCF-debug
#else
      Debug=.false.
#endif


      FactC = one

      DoRead  = .false.
      IREDC= -1  ! unknwn reduced set

      If (nDen.lt.1) then
         write(6,*)SECNAM//'Invalid parameter nDen= ',nDen
         call abend()
      EndIf


        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time read/transform vectors
           tcoul(i) = zero  !time for computing Coulomb
           texch(i) = zero  !time for computing Exchange
        end do

C ==================================================================

c --- Various offsets
c --------------------
        ISTLT(1)=0
      DO ISYM=2,NSYM
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive D and F matrices
      END DO

      Call GetMem('ipRedM','Allo','Inte',ipDab,2*nDen)
      ipFab = ipDab + nDen - 1

      Call GetMem('ip_Lab','Allo','Inte',ip_Lab,8*nDen)
      Call GetMem('ip_Orb','Allo','Inte',ip_Orb,8*nDen)

      DO jDen=1,nDen

         iWork(ip_Orb+8*(jDen-1)) = ipMOs(jDen)

         DO ISYM=2,NSYM

            iWork(ip_Orb-1+8*(jDen-1)+iSym) = ipOrb(iSym-1,jDen)
     &                       + nMOs(iSym-1,jDen)*nBas(iSym-1)
         END DO

      END DO

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

         If (NumCho(jSym).lt.1) GOTO 1000


C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec = 0  ! mem for storing the half-transformed vec

         do l=1,nSym
            k=Muld2h(l,JSYM)
            Mmax = 0
            do jDen=1,nDen
               if (ExFac(jDen).ne.zero) Mmax = Max(Mmax,nMOs(k,jDen))
            end do
            mTvec = mTvec + nBas(l)*Mmax
         end do

         mTvec=Max(mTvec,1)

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
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                       '    rc= ',irc
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            If(JSYM.eq.1)Then

              Do jDen=1,nDen
                Call GetMem('rsD','Allo','Real',iWork(ipDab+jDen-1),nRS)
                Call GetMem('rsF','Allo','Real',iWork(ipFab+jDen-1),nRS)
                Call Fzero(Work(ipDr(jDen)),nRS)
                Call Fzero(Work(ipFr(jDen)),nRS)
             End Do

            EndIf

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec  = Min(LWORK/(nRS+mTvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec
               WRITE(6,*) 'jsym= ',jsym
               rc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
            Call GetMem('ChoT','Allo','Real',ipChoT,mTvec*nVec)

            If(JSYM.eq.1)Then
C --- Transform the densities to reduced storage
              mode = 'toreds'
              add  = .false.
              Call swap_rs2full(irc,iLoc,nDen,JSYM,ISTLT,
     &                              ipDLT,iWork(ipDab),mode,add)
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
C ************ (alpha+beta) COULOMB CONTRIBUTION  ****************
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * DI(rs)
C==========================================================
C
                  CALL CWTIME(TCC1,TWC1)

                  ipVJ = ipChoT

                  Do jDen=1,nDen

                     CALL DGEMV_('T',nRS,JNUM,
     &                          ONE,Work(ipLrs),nRS,
     &                 Work(ipDr(jDen)),1,ZERO,Work(ipVJ),1)

C --- F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
C===============================================================

                     Fact = dble(min(jVec-iVrs,1))

                     CALL DGEMV_('N',nRS,JNUM,
     &                          FactC,Work(ipLrs),nRS,
     &                 Work(ipVJ),1,Fact,Work(ipFr(jDen)),1)


                  End Do

                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)

               EndIf  ! Coulomb contribution


C *************** EXCHANGE CONTRIBUTIONS  ***********************

C --- Set pointers (ipLab) to the half-transformed Cholesky vectors
               Do jDen=1,nDen

                  lChoT=0
                  Do iSymb=1,nSym

                     iSymp = MulD2h(JSYM,iSymb)

                     iWork(ip_Lab - 1 + 8*(jDen-1)
     &                                + iSymp) = ipChoT + lChoT ! LpJ,b

                     lChoT = lChoT + nMOs(iSymp,jDen)*JNUM*nBas(iSymb)

                  End Do

               End Do


               iSwap = 2  ! LpJ,b are returned

               Do jDen=1,nDen

                 CALL CWTIME(TCR3,TWR3)

                 IF (ExFac(jDen).ne.zero) THEN

                     FactX = xone*ExFac(jDen)

C --- Set up the skipping flags
C -------------------------------------------------------------
                  Do i=1,nSym

                     k = Muld2h(i,JSYM)
                     iSkip(k) = Min(1,nBas(i)*nMOs(k,jDen))

                  End Do
C -------------------------------------------------------------


C *********************** HALF-TRANSFORMATION  ****************

                  CALL CHO_X_getVtra(irc,Work(ipLrs),LREAD,jVEC,JNUM,
     &                         JSYM,iSwap,IREDC,jDen,jDen,iWork(ip_Orb),
     &                         nMOs,iWork(ip_Lab),iSkip,DoRead)

                  CALL CWTIME(TCR4,TWR4)
                  tread(1) = tread(1) + (TCR4 - TCR3)
                  tread(2) = tread(2) + (TWR4 - TWR3)

#ifdef _DEBUGPRINT_
      write(6,*) 'Half-transformation in the MO space: ',jDen
      write(6,*) 'Total allocated :     ',mTvec*nVec,' at ',ipChoT
      write(6,*) 'Mem pointers ipLab: ',((ipLab(i,j),i=1,nSym),j=1,jDen)
      write(6,*) 'iSkip :             ',(iSkip(k),k=1,nSym)
      write(6,*) 'LREAD: ',LREAD,' allocated at ',ipLrs
      write(6,*) 'JSYM :                ',JSYM
      write(6,*) 'JRED :                ',JRED
      write(6,*) 'JNUM :                ',JNUM
#endif

                  if (irc.ne.0) then
                     rc = irc
                     RETURN
                  endif


                  CALL CWTIME(TCX1,TWX1)

                  Do iSyma=1,nSym

                     iSymk = MulD2h(JSYM,iSyma)

C ---------------------------------------------------------------------
c *** Compute only the LT part of the Exchange matrix ********
C
C     F(ab) = F(ab) + FactX * sum_Jk  LkJ,a * LkJ,b
C ---------------------------------------------------------------------
                     NK = nMOs(iSymk,jDen)

                     If (iSkip(iSymk).ne.0) Then

                        ISFI = ipFLT(jDen) + ISTLT(iSyma)

                        CALL DGEMM_TRI('T','N',nBas(iSyma),nBas(iSyma),
     &                         NK*JNUM,FactX,Work(ipLab(iSymk,jDen)),
     &                         NK*JNUM,Work(ipLab(iSymk,jDen)),NK*JNUM,
     &                         One,Work(ISFI),nBas(iSyma))


c          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMA
c          CALL TRIPRT('FI: jSym>1',' ',Work(ISFI),nBas(iSyma))

                     EndIf

C --------------------------------------------------------------------
                  End Do  !loop over MOs symmetries

                  CALL CWTIME(TCX2,TWX2)
                  texch(1) = texch(1) + (TCX2 - TCX1)
                  texch(2) = texch(2) + (TWX2 - TWX1)


                 ENDIF  ! ExFac(jDen)=0.0d0

               End Do   ! loop over densities

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop


            If(JSYM.eq.1)Then

             if (nVrs.gt.0) then
c --- backtransform fock matrix in full storage
               mode = 'tofull'
               add  = JRED.gt.JRED1
               Call swap_rs2full(irc,iLoc,nDen,JSYM,ISTLT,
     &                               ipFLT,iWork(ipFab),mode,add)
             endif

            EndIf


C --- free memory
            Call GetMem('ChoT','Free','Real',ipChoT,mTvec*nVec)
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

            If(JSYM.eq.1)Then
              Do jDen=nDen,1,-1
                Call GetMem('rsD','Free','Real',iWork(ipDab+jDen-1),nRS)
                Call GetMem('rsF','Free','Real',iWork(ipFab+jDen-1),nRS)
              End Do
            EndIf


999         Continue

         END DO   ! loop over red sets

1000     CONTINUE

      END DO   !loop over JSYM

      Call GetMem('ip_Orb','Free','Inte',ip_Orb,8*nDen)
      Call GetMem('ip_Lab','Free','Inte',ip_Lab,8*nDen)
      Call GetMem('ipRedM','Free','Inte',ipDab,2*nDen)

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

         Write(6,'(2x,A26,2f10.2)')'READ/TRANSFORM VECTORS           '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
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
      Do jDen=1,nDen
       WRITE(6,'(6X,A)')'***** FOCK MATRIX AO-BASIS ***** '
       WRITE(6,'(6X,A)')'***** computed from Density #',jDen
        DO ISYM=1,NSYM
           ISFI=ipFLT(jDen)+ISTLT(ISYM)
           IF( NBAS(ISYM).GT.0 ) THEN
             WRITE(6,'(6X,A)')
             WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
             call TRIPRT('','',Work(ISFI),NBAS(ISYM))
           ENDIF
        END DO
      END DO

      endif

#endif

      rc  = 0


      Return
      END

**************************************************************
**************************************************************
