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

      SUBROUTINE CHO_FMCSCF(rc,ipFA,ipFI,nForb,nIorb,nAorb,FactXI,
     &                      ipDI,ipDA1,DoActive,POrb,nChM,ipInt,ExFac)

**********************************************************************
*  Author : F. Aquilante
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
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use Data_structures, only: DSBA_Type, SBA_Type
      use Data_structures, only: Allocate_SBA, Deallocate_SBA
      use Data_structures, only: twxy_Type
      use Data_structures, only: Allocate_twxy, Deallocate_twxy
      Implicit Real*8 (a-h,o-z)

      Type (DSBA_Type) POrb(3)
      Type (SBA_Type), Target:: Laq(3), Lxy
      Type (twxy_type) Scr

      Integer   rc
      Integer   iSkip(8)
      Integer   ISTLT(8)
      Real*8    tread(2),tcoul(2),texch(2),tintg(2), ExFac
      Integer   ipDA1,ipDI
      Integer   ipFA,ipFI
      Integer   ipDLT(2),ipFLT(2),ipDab(2),ipFab(2)
      Integer   nForb(8),nIorb(8),nAorb(8),nPorb(8),nnA(8,8),nChM(8)
#ifdef _DEBUGPRINT_
      Logical   Debug
#endif
      Logical   DoRead,DoTraInt,DoActive
      Character*50 CFmt
      Character(LEN=10), Parameter:: SECNAM = 'CHO_FMCSCF'
#include "chotime.fh"

#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Real*8, Parameter:: FactCI = One, FactCA = One, FactXA = -Half

      Logical add
      Character*6 mode

      Integer nAux(8)
      Real*8, Allocatable:: Lrs(:,:), Drs(:,:), Frs(:,:)
      Real*8, Pointer:: VJ(:)=>Null()


************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
************************************************************************

#ifdef _DEBUGPRINT_
      Debug=.false.! to avoid double printing in CASSCF-debug
#endif
      if(ExFac.ne.1.0d0) then
          write(6,*) 'WARNING: if you are running MCPDFT calculations'
          write(6,*) 'and end up with this message, you are in trouble.'
          write(6,*) 'Tweaks are needed!! Please contact'
          write(6,*) 'Giovanni Li Manni'
          write(6,*) 'giovannilimanni@gmail.com'
      end if
      DoRead  = .false.
      DoTraInt = .false.
      IREDC = -1  ! unknown reduced set in core

      ipDLT(1) = ipDI    ! some definitions
      ipDLT(2) = ipDA1
      ipFLT(1) = ipFI
      ipFLT(2) = ipFA

      nDen = 1
      if (DoActive) nDen=2

        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time read/transform vectors
           tcoul(i) = zero  !time for computing Coulomb
           texch(i) = zero  !time for computing Exchange
           tintg(i) = zero  !time for computing (pu|vx) integrals
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
      ISTLT(1)=0
      DO ISYM=2,NSYM
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB
      END DO

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

         If (NumCho(jSym).lt.1) GOTO 1000

         iCase = 1 ! (wa|xy)
         Call Allocate_twxy(Scr,nAorb,nBas,JSYM,nSym,iCase)

C --- Set up the skipping flags --------
C -------------------------------------------------------------
         Do i=1,nSym
            k = Muld2h(i,JSYM)
            iSkip(i) = Min(1,nBas(i)*nBas(k)) ! skip Lik vector
            iSkip(i) = iSkip(i)*(nPorb(i)+nChM(i))
         End Do
C -------------------------------------------------------------


C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec1 = 0  ! mem for storing the half-transformed vec
         mTvec2 = 0  ! mem for storing the half-transformed vec
         mTvec3 = 0  ! mem for storing the half-transformed vec
         mTvec4 = 0  ! mem for storing the half-transformed vec

         nAux(:) = nForb(:) + nIorb(:)
         do l=1,nSym
            k=Muld2h(l,JSYM)
            mTvec1= mTvec1+ nBas(k)*nAux(l)
            mTvec2= mTvec2+ nBas(k)*nChM(l)
            mTvec3= mTvec3+ nBas(k)*nAorb(l)
            mTvec4= mTvec4+ nnA(k,l)
         end do

         mTvec = Max(mTvec1,mTvec2,mTvec3+mTvec4,1)

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

               if(DoActive)then
                 Call mma_allocate(Drs,nRS,2,Label='Drs')
                 ipDab(1) = ip_of_Work(Drs(1,1))
                 ipDab(2) = ip_of_Work(Drs(1,2))
                 Call mma_allocate(Frs,nRS,2,Label='Frs')
                 ipFab(1) = ip_of_Work(Frs(1,1))
                 ipFab(2) = ip_of_Work(Frs(1,2))
               Else
                 Call mma_allocate(Drs,nRS,1,Label='Drs')
                 ipDab(1) = ip_of_Work(Drs(1,1))
                 Call mma_allocate(Frs,nRS,1,Label='Frs')
                 ipFab(1) = ip_of_Work(Frs(1,1))
               endif

            EndIf

            Call mma_maxDBLE(LWORK)

            nVec  = Min(LWORK/(nRS+mTvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec
               WRITE(6,*) 'reading ',nRS,' and transforming to ',mTvec
               WRITE(6,*) 'of jsym= ',jsym,' and JRED= ',JRED
               rc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')


            If(JSYM.eq.1)Then
C --- Transform the density to reduced storage
               mode = 'toreds'
               add  = .false.
               Call swap_rs2full(irc,iLoc,nDen,ipDLT,ipDab,mode,add)
            EndIf

C --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch
               iSwap = 2  ! LpJ,b are returned
               Call Allocate_SBA(Laq(1),nAux,nBas,nVec,JSYM,nSym,iSwap)
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

                  VJ(1:JNUM) => Laq(1)%A0(1:JNUM)

                  CALL DGEMV_('T',nRS,JNUM,
     &                 ONE,Lrs,nRS,
     &                 Drs(:,1),1,ZERO,VJ,1)

C --- FI(rs){#J} <- FI(rs){#J} + FactCI * sum_J L(rs,{#J})*V{#J}
C===============================================================

                  xfac = dble(min(jVec-iVrs,1))

                  CALL DGEMV_('N',nRS,JNUM,
     &                 FactCI,Lrs,nRS,
     &                 VJ,1,xfac,Frs(:,1),1)


                  If (DoActive) Then
C ************ ACTIVE COULOMB CONTRIBUTION  ****************
C
C---- Computing the intermediate vector V(J)
C ---------------------------------------
C --- U{#J} <- U{#J}  +  sum_rs  L(rs,{#J}) * DA(rs)
C==========================================================
C
                     CALL DGEMV_('T',nRS,JNUM,
     &                          ONE,Lrs,nRS,
     &                          Drs(:,2),1,ZERO,VJ,1)

C --- FA(rs){#J} <- FA(rs){#J} + FactCA * sum_J L(rs,{#J})*U{#J}
C===============================================================

                     xfac = dble(min(jVec-iVrs,1))

                     CALL DGEMV_('N',nRS,JNUM,
     &                          FactCA,Lrs,nRS,
     &                          VJ,1,xfac,Frs(:,2),1)

                  EndIf


                  CALL CWTIME(TCC2,TWC2)
                  tcoul(1) = tcoul(1) + (TCC2 - TCC1)
                  tcoul(2) = tcoul(2) + (TWC2 - TWC1)

                  VJ=>Null()

               ENDIf  ! Coulomb (jsym=1)

C ************ BEGIN EXCHANGE CONTRIBUTIONS  ****************

C *********************** INACTIVE HALF-TRANSFORMATION  ****************

               kMOs = 1  ! inactive MOs
               nMOs = 1

               CALL CWTIME(TCR3,TWR3)

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                         JSYM,iSwap,IREDC,nMOs,kMOs,POrb,
     &                         Laq,DoRead)


               CALL CWTIME(TCR4,TWR4)
               tread(1) = tread(1) + (TCR4 - TCR3)
               tread(2) = tread(2) + (TWR4 - TWR3)


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
     &                         NK*JNUM,FactXI,Laq(1)%SB(iSymk)%A3,
     &                         NK*JNUM,Laq(1)%SB(iSymk)%A3,NK*JNUM,
     &                         One,Work(ISFI),nBas(iSyma))


c          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMA
c          CALL TRIPRT('FI: ',' ',Work(ISFI),nBas(iSyma))

                  EndIf

C --------------------------------------------------------------------
               End Do  !loop over MOs symmetries

               CALL CWTIME(TCX2,TWX2)
               texch(1) = texch(1) + (TCX2 - TCX1)
               texch(2) = texch(2) + (TWX2 - TWX1)

               Call Deallocate_SBA(Laq(1))

               IF(DoActive)THEN
                  iSwap = 2  ! LxJ,b are returned
                  Call Allocate_SBA(Laq(2),nChM,nBas,nVec,JSYM,nSym,
     &                              iSwap)
C *********************** "CHOLESKY" HALF-TRANSFORMATION  ****************
C --------------------------------------------------------------------
c --- Using "Cholesky MOs" obtained by cholesky decomposing DA
C --------------------------------------------------------------------

                  CALL CWTIME(TCR5,TWR5)

                  kMOs = 2  ! Cholesky MOs
                  nMOs = 2

                  CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,POrb,
     &                            Laq,DoRead)


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
     &                         NAch*JNUM,FactXA,Laq(2)%SB(iSymw)%A3,
     &                         NAch*JNUM,Laq(2)%SB(iSymw)%A3,NAch*JNUM,
     &                         One,Work(ISFA),nBas(iSyma))


c          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYMA
c          CALL TRIPRT('FA: ',' ',Work(ISFA),nBas(iSyma))

                     EndIf

C --------------------------------------------------------------------
                    End Do  !loop over MOs symmetries

                  CALL CWTIME(TCX4,TWX4)
                  texch(1) = texch(1) + (TCX4 - TCX3)
                  texch(2) = texch(2) + (TWX4 - TWX3)

                  Call Deallocate_SBA(Laq(2))
               ENDIF   ! Do Active Exchange
C ************  END EXCHANGE CONTRIBUTIONS  ****************


C --------------------------------------------------------------------
C --- First half Active transformation  Lvb,J = sum_a  C(v,a) * Lab,J
C --------------------------------------------------------------------
               ! Lvw,J, LT-storage for the diagonal symmetry blocks
               iSwap = 4
               Call Allocate_SBA(Lxy,nAorb,nAorb,nVec,JSYM,nSym,iSwap)
*              Call mma_allocate(Lxy%Lxy_full,mTvec4*nVec,Label='Lxy')

               iSwap = 0  ! Lvb,J are returned
               Call Allocate_SBA(Laq(3),nAorb,nBas,nVec,JSYM,nSym,iSwap)

               CALL CWTIME(TCR7,TWR7)

               kMOs = 3  ! Active MOs
               nMOs = 3  ! Active MOs

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                            JSYM,iSwap,IREDC,nMOs,kMOs,POrb,
     &                            Laq,DoRead)


               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
C --------------------------------------------------------------------
               IF (JSYM.eq.1) THEN   !  Lvw,J in LT-storage

                  Do iSyma=1,nSym

                     NAv = nAorb(iSyma)

                     If(NAv.ne.0)Then

                      Do JVC=1,JNUM

                       CALL DGEMM_Tri('N','T',NAv,NAv,nBas(iSyma),
     &                             One,Laq(3)%SB(iSyma)%A3(:,:,JVC),NAv,
     &                                    POrb(3)%SB(iSyma)%A2,NAv,
     &                               Zero,Lxy%SB(iSyma)%A2(:,JVC),NAv)

                      End Do

                     EndIf

                  End Do


               ELSE

C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C(w,b)
C --------------------------------------------------------------------
                  Do iSymb=1,nSym

                     iSymv = MulD2h(JSYM,iSymb)
                     NAv = nAorb(iSymv)
                     NAw = nAorb(iSymb) ! iSymb=iSymw

                     If(NAv*NAw.ne.0.and.iSymv.lt.iSymb)Then

                      Do JVC=1,JNUM

                       CALL DGEMM_('N','T',NAv,NAw,nBas(iSymb),
     &                            One,Laq(3)%SB(iSymv)%A3(:,:,JVC),NAv,
     &                                POrb(3)%SB(iSymb)%A2,NAw,
     &                           Zero,Lxy%SB(iSymv)%A2(:,JVC),NAv)

                      End Do

                     EndIf

                  End Do


               ENDIF   ! jSym >=< 1

               CALL CWTIME(TCR8,TWR8)
               tread(1) = tread(1) + (TCR8 - TCR7)
               tread(2) = tread(2) + (TWR8 - TWR7)


C *************** EVALUATION OF THE (WA|XY) INTEGRALS ***********

               CALL CWTIME(TCINT1,TWINT1)

               DoTraInt = JRED.eq.JRED2.and.iBatch.eq.nBatch

               CALL CHO_eval_waxy(irc,Scr,Laq(3),Lxy,ipInt,
     &                            nAorb,JSYM,JNUM,DoTraInt)

               CALL CWTIME(TCINT2,TWINT2)
               tintg(1) = tintg(1) + (TCINT2 - TCINT1)
               tintg(2) = tintg(2) + (TWINT2 - TWINT1)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

C --------------------------------------------------------------------
C --------------------------------------------------------------------
               Call Deallocate_SBA(Lxy)
               Call Deallocate_SBA(Laq(3))

            END DO  ! end batch loop


            if (JSYM.eq.1) then
c --- backtransform fock matrix in full storage
               mode = 'tofull'
               add  = .true.
               Call swap_rs2full(irc,iLoc,nDen,ipFLT,ipFab,mode,add)
            endif

C --- free memory
            Call mma_deallocate(Lrs)

            if (JSYM.eq.1) then
               Call mma_deallocate(Drs)
               Call mma_deallocate(Frs)
            endif



999         CONTINUE

         END DO   ! loop over red sets

         Call Deallocate_twxy(Scr)

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
         Write(6,'(2x,A26,2f10.2)')'(WA|XY) INTEGRALS                '
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

      rc  = 0


      Return
      END
