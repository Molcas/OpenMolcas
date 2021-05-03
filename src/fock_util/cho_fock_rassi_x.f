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
      SUBROUTINE CHO_FOCK_RASSI_X(DLT,MO1,MO2,FLT,FSQ,ipInt)

**********************************************************************
*  Author : F. Aquilante
*
*
*  Note:  this routine differs from CHO_FOCK_RASSI because it can
*         handle ALSO the case where the 2 sets of MOs are different!
*         The Exchange contribution is non-symmetric and so is FI
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
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use Data_Structures, only: DSBA_Type, SBA_Type
      use Data_Structures, only: Allocate_SBA, Deallocate_SBA
      use Data_Structures, only: twxy_Type
      use Data_Structures, only: Allocate_twxy, Deallocate_twxy
      Implicit Real*8 (a-h,o-z)

      Type (DSBA_Type) DLT, MO1(2), MO2(2), FLT, FSQ
      Type (SBA_Type), Target:: Laq(2)
      Type (Twxy_Type) Scr

      Integer   rc
      Integer   iSkip(8)
      Integer   ISTLT(8), ISTSQ(8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2)
#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   DoReord, add
      Character*50 CFmt
      Character(LEN=16), Parameter:: SECNAM = 'CHO_FOCK_RASSI_X'
#include "chotime.fh"
#include "real.fh"

      Logical, Parameter :: DoRead = .false.
      Real*8, Parameter ::FactCI = One, FactXI = -One
      Character*6 mode
#include "cho_jobs.fh"

#include "rassi.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"

      Real*8, Allocatable:: Lrs(:,:), Drs(:), Frs(:)

      Real*8, Pointer:: VJ(:)

**************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
**************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in CASSCF-debug
#endif
      DoReord = .false.
      IREDC = -1  ! unknown reduced set in core

      nDen=2
      If (Fake_CMO2) nDen = 1  ! MO1 = MO2
      kDen=nDen
      ipFLT = ip_of_Work(FLT%A0(1))
      ipK   = ip_of_Work(FSQ%A0(1))

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      ! 1 --> CPU   2 --> Wall
      tread(:) = zero  !time read/transform vectors
      tcoul(:) = zero  !time for computing Coulomb
      texch(:) = zero  !time for computing Exchange
      tintg(:) = zero  !time for computing (tw|xy) integrals

C ==================================================================

c --- Various offsets
c --------------------
      ISTLT(1)=0
      ISTSQ(1)=0
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB ! Inactive Coul matrix
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NB**2 ! Inactive Exch matrix
      END DO

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

        If (NumCho(jSym).lt.1) GOTO 1000

        iCase = 0 ! twxy
        Call Allocate_twxy(Scr,nAsh,nAsh,JSYM,nSym,iCase)

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
               call abend()
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                        '   rc= ',irc
              call abend()
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

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               add = .False.
               mDen=1
               ipDLT = ip_of_Work(DLT%A0(1))
               Call swap_rs2full(irc,iLoc,nRS,mDen,JSYM,[ipDLT],Drs,
     &                           mode,add)
            EndIf

C --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch

               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               else
                  JNUM = nVec
               endif

               iSwap = 2  ! LpJ,b are returned
               Do jDen = 1, nDen
                  Call Allocate_SBA(Laq(jDen),nIsh,nBas,nVec,JSYM,nSym,
     &                              iSwap)
               End Do

               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CWTIME(TCR1,TWR1)

               CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,
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

                  VJ(1:JNUM) => Laq(1)%A0(1:JNUM)

                  CALL DGEMV_('T',nRS,JNUM,
     &                 ONE,Lrs,nRS,
     &                 Drs,1,ZERO,VJ,1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  Fact = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                 FactCI,Lrs,nRS,
     &                 VJ,1,Fact,Frs,1)


                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)

                  VJ=>Null()

               EndIf  ! Coulomb contribution


C *************** EXCHANGE CONTRIBUTIONS  ***********************

               CALL CWTIME(TCR3,TWR3)

               kMOs = 1
               nMOs = nDen

C --- Set up the skipping flags
C -------------------------------------------------------------
               Do i=1,nSym

                  k = Muld2h(i,JSYM)
                  iSkip(k) = Min(1,nIsh(k)*NBAS(i))

               End Do
C -------------------------------------------------------------


C *********************** HALF-TRANSFORMATION  ****************

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,MO1,
     &                            Laq,DoRead)


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

                     ISFI = ipK + ISTSQ(iSyma)

                     CALL DGEMM_('T','N',NBAS(iSyma),NBAS(iSyma),
     &                         NK*JNUM,FactXI,Laq(kDen)%SB(iSymk)%A3,
     &                         NK*JNUM,Laq(1)%SB(iSymk)%A3,NK*JNUM,
     &                             One,Work(ISFI),NBAS(iSyma))

                  EndIf

C --------------------------------------------------------------------
               End Do  !loop over MOs symmetries

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)

               Do jDen = 1, nDen
                  Call Deallocate_SBA(Laq(jDen))
               End Do

C ************  END EXCHANGE CONTRIBUTION  ****************

               iSwap = 0  ! Lvb,J are returned
               Call Allocate_SBA(Laq(1),nAsh,nBas,nVec,JSYM,nSym,iSwap)
               Call Allocate_SBA(Laq(2),nAsh,nAsh,nVec,JSYM,nSym,iSwap)

C --------------------------------------------------------------------
C --- First half Active transformation  Lvb,J = sum_a  C1(v,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCR7,TWR7)

C --- Set up the skipping flags
C -------------------------------------------------------------
               Do i=1,nSym

                  k = Muld2h(i,JSYM)
                  iSkip(k) = Min(1,NBAS(i)*nAsh(k))

               End Do

               kMOs = 1  !
               nMOs = 1  ! Active MOs (1st set)

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                           JSYM,iSwap,IREDC,nMOs,kMOs,MO2,
     &                           Laq,DoRead)

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

                       CALL DGEMM_('N','T',NAv,NAw,NBAS(iSymb),
     &                            One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                                MO2(kDen)%SB(iSymb)%A2,NAw,
     &                           Zero,Laq(2)%SB(iSymv)%A3(:,:,JVC),NAv)

                      End Do

                     EndIf

                  End Do


C
C
C *************** EVALUATION OF THE (TW|XY) INTEGRALS ***********

               CALL CWTIME(TCINT1,TWINT1)

               DoReord = JRED.eq.JRED2.and.iBatch.eq.nBatch

               CALL CHO_rassi_twxy(irc,Scr,Laq(2),ipInt,nAsh,
     &                                 JSYM,JNUM,DoReord)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C ---------------- END (TW|XY) EVALUATION -----------------------


               Call Deallocate_SBA(Laq(2))
               Call Deallocate_SBA(Laq(1))
            END DO  ! end batch loop


            If(JSYM.eq.1)Then
c --- backtransform fock matrix to full storage
               mode = 'tofull'
               add = .True.
               mDen=1
               Call swap_rs2full(irc,iLoc,nRS,mDen,JSYM,[ipFLT],Frs,
     &                           mode,add)
            EndIf

C --- free memory
            Call mma_deallocate(Lrs)

            If(JSYM.eq.1)Then
              Call mma_deallocate(Frs)
              Call mma_deallocate(Drs)
            EndIf


999         Continue

         END DO   ! loop over red sets

         Call Deallocate_twxy(Scr)

1000  CONTINUE


      END DO  ! loop over JSYM

* --- Accumulate Coulomb and Exchange contributions
      Do iSym=1,nSym

         ipFI = ipFLT - 1 + ISTLT(iSym)
         ipKI = ipK -1 + ISTSQ(iSym)

         Do ia=1,nBas(iSym)
            Do ib=1,ia-1
               iabt = ipFI + ia*(ia-1)/2 + ib
               iabq = ipKI + nBas(iSym)*(ia-1) + ib
               Work(iabq)=Work(iabq)+Work(iabt)
               iabq = ipKI + nBas(iSym)*(ib-1) + ia
               Work(iabq)=Work(iabq)+Work(iabt)
            End Do
            iabt = ipFI + ia*(ia+1)/2
            iabq = ipKI + nBas(iSym)*(ia-1) + ia
            Work(iabq)=Work(iabq)+Work(iabt)
         End Do

      End Do



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
        ISFI=ipK+ISTSQ(ISYM)
        IF( NBAS(ISYM).GT.0 ) THEN
          WRITE(6,'(6X,A)')'***** INACTIVE FOCK MATRIX ***** '
          WRITE(6,'(6X,A)')
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          call CHO_OUTPUT(Work(ISFI),1,NBAS(ISYM),1,NBAS(ISYM),
     &                    NBAS(ISYM),NBAS(ISYM),1,6)
        ENDIF
      END DO

      endif

#endif

      rc  = 0


      Return
      END

**************************************************************
