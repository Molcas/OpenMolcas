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
      SUBROUTINE CHO_FOCK_DFT_RED(irc,DLT,FLT)

!********************************************************
!
! Author:  F. Aquilante
!
! Coulomb fock matrix  (level 2 BLAS)
!
! --- F(ab) = 2 * sum_J  Lab,J * sum_gd  D(gd) * Lgd,J
!
!********************************************************
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use Data_Structures, only: DSBA_Type
      Implicit Real*8 (a-h,o-z)
#ifdef _DEBUGPRINT_
      Logical Debug
#endif
      Logical add
      Type (DSBA_Type) DLT, FLT
      Real*8  tread(2),tcoul(2)
      Character*16  SECNAM
      Character*6   mode
      Character*50 CFmt
#include "chotime.fh"

      parameter (SECNAM = 'CHO_FOCK_DFT_RED')

#include "real.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Real*8, Allocatable :: Lrs(:,:), Drs(:), Frs(:), VJ(:)

#ifdef _DEBUGPRINT_
      Debug=.true.
#endif

      FactC = one

! --- For Coulomb only, the vectors symmetry is restricted to 1
      JSYM=1
      If (NumCho(JSYM).lt.1) Return

        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time for reading the vectors
           tcoul(i) = zero  !time for computing Coulomb
        end do

      iLoc = 3 ! use scratch location in reduced index arrays

      JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
      JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

      Do JRED=JRED1,JRED2

! --- Memory management section -----------------------------
! ---
      CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

      if (nVrs.eq.0) goto 999

      if (nVrs.lt.0) then
         Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs < 0. STOP!!'
         call abend()
      endif

      Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
      if(irc.ne.0)then
        Write(6,*)SECNAM//'cho_X_setred non-zero return code. rc= ',irc
        call abend()
      endif

      nRS = nDimRS(JSYM,JRED)

      Call mma_allocate(Drs,nRS,Label='Drs')
      Call mma_allocate(Frs,nRS,Label='Frs')
      Drs(:)=Zero
      Frs(:)=Zero

      Call mma_maxDBLE(LWork)

      nVec  = Min(LWORK/(nRS+1),nVrs)

      If (nVec.lt.1) Then
         WRITE(6,*) SECNAM//': Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'min. mem. need= ',nRS+1
         irc = 33
         CALL Abend()
         nBatch = -9999  ! dummy assignment
      End If

      LREAD = nRS*nVec

      Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')
      Call mma_allocate(VJ,nVec,Label='VJ')

! --- Transform the density to reduced storage
      mode = 'toreds'
      add  = .false.
      nDen=1
      Call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,[DLT],Drs,mode,add)

! --- BATCH over the vectors in JSYM=1 ----------------------------

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

         CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,JRED,MUSED)

         If (NUMV.le.0 .or. NUMV.ne.JNUM) then
            irc=77
            RETURN
         End If

         CALL CWTIME(TCR2,TWR2)
         tread(1) = tread(1) + (TCR2 - TCR1)
         tread(2) = tread(2) + (TWR2 - TWR1)

! ************ BEGIN COULOMB CONTRIBUTION  ****************
!
!---- Computing the intermediate vector V(J)
!
! --- Contraction with the density matrix
! ---------------------------------------
! --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
!==========================================================
!
         CALL CWTIME(TCC1,TWC1)

         CALL DGEMV_('T',nRS,JNUM,                                      &
     &              ONE,Lrs,nRS,                                        &
     &              Drs,1,ZERO,VJ,1)

! --- Frs{#J} <- Frs{#J} + sum_J L(rs,{#J})*V{#J}
!==========================================================

         xfac = dble(min(jVec-iVrs,1))

         CALL DGEMV_('N',nRS,JNUM,                                      &
     &              FactC,Lrs,nRS,                                      &
     &              VJ,1,xfac,Frs,1)


         CALL CWTIME(TCC2,TWC2)
         tcoul(1) = tcoul(1) + (TCC2 - TCC1)
         tcoul(2) = tcoul(2) + (TWC2 - TWC1)


      END DO  !end batch loop

      if (nVrs.gt.0) then
! --- backtransform fock matrix in full storage
         mode = 'tofull'
         add  = JRED.gt.JRED1
         Call swap_rs2full(irc,iLoc,nRS,nDen,JSYM,[FLT],Frs,mode,add)
      endif

! --- free memory
      Call mma_deallocate(VJ)
      Call mma_deallocate(Lrs)
      Call mma_deallocate(Frs)
      Call mma_deallocate(Drs)


999   Continue

      END DO   ! loop over red sets


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

!
!---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky SCF timing from '//SECNAM
      Write(6,CFmt)'-----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '  &
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '  &
     &                           //'         ',tcoul(1),tcoul(2)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '  &
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


! Print the Fock-matrix
#ifdef _DEBUGPRINT_
      if(Debug) then !to avoid double printing in SCF-debug

      WRITE(6,'(6X,A)')'TEST PRINT FROM '//SECNAM
      WRITE(6,'(6X,A)')
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT('Coulomb Fmat',' ',FLT%SB(ISYM)%A1,NB)
        END IF
      END DO

      endif

#endif


      irc=0

      Return
      End
