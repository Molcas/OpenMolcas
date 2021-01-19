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
      SUBROUTINE CHO_FOCK_DFT_RED(irc,DLT,FLT)

C********************************************************
C
C Author:  F. Aquilante
C
C Coulomb fock matrix  (level 2 BLAS)
C
C --- F(ab) = 2 * sum_J  Lab,J * sum_gd  D(gd) * Lgd,J
C
C********************************************************
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
#ifdef _DEBUGPRINT_
      Logical Debug
#endif
      Logical add,timings
      Real*8  DLT(*),FLT(*)
      Real*8  tread(2),tcoul(2)
      Character*16  SECNAM
      Character*6   mode
      Character*50 CFmt
      COMMON    /CHOTIME /timings

      parameter (SECNAM = 'CHO_FOCK_DFT_RED')
      parameter (zero = 0.0d0, one = 1.0d0)

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

#ifdef _DEBUGPRINT_
      Debug=.true.
#endif


      FactC = one

C --- For Coulomb only, the vectors symmetry is restricted to 1
      JSYM=1
      If (NumCho(JSYM).lt.1) Return

        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero  !time for reading the vectors
           tcoul(i) = zero  !time for computing Coulomb
        end do


      ipDLT = ip_of_Work(DLT(1))
      ipFLT = ip_of_Work(FLT(1))

      iLoc = 3 ! use scratch location in reduced index arrays

      JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
      JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

      Do JRED=JRED1,JRED2

C --- Memory management section -----------------------------
C ---
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

      Call GetMem('rsD','Allo','Real',ipDab,nRS)
      Call GetMem('rsF','Allo','Real',ipFab,nRS)
      Call Fzero(Work(ipDab),nRS)
      Call Fzero(Work(ipFab),nRS)

      Call GetMem('MaxM','Max','Real',KDUM,LWORK)

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

      Call GetMem('rsL','Allo','Real',ipLab,LREAD)
      Call GetMem('VJ','Allo','Real',ipVJ,nVec)

C --- Transform the density to reduced storage
      mode = 'toreds'
      add  = .false.
      Call switch_sto(irc,iLoc,ipDLT,ipDab,mode,add)

C --- BATCH over the vectors in JSYM=1 ----------------------------

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

         CALL CHO_VECRD(Work(ipLab),LREAD,JVEC,IVEC2,JSYM,
     &                  NUMV,JRED,MUSED)

         If (NUMV.le.0 .or. NUMV.ne.JNUM) then
            irc=77
            RETURN
         End If

         CALL CWTIME(TCR2,TWR2)
         tread(1) = tread(1) + (TCR2 - TCR1)
         tread(2) = tread(2) + (TWR2 - TWR1)

C ************ BEGIN COULOMB CONTRIBUTION  ****************
C
C---- Computing the intermediate vector V(J)
C
C --- Contraction with the density matrix
C ---------------------------------------
C --- V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
C==========================================================
C
         CALL CWTIME(TCC1,TWC1)

         CALL DGEMV_('T',nRS,JNUM,
     &              ONE,Work(ipLab),nRS,
     &              Work(ipDab),1,ZERO,Work(ipVJ),1)

C --- Frs{#J} <- Frs{#J} + sum_J L(rs,{#J})*V{#J}
C==========================================================

         xfac = dble(min(jVec-iVrs,1))

         CALL DGEMV_('N',nRS,JNUM,
     &              FactC,Work(ipLab),nRS,
     &              Work(ipVJ),1,xfac,Work(ipFab),1)


         CALL CWTIME(TCC2,TWC2)
         tcoul(1) = tcoul(1) + (TCC2 - TCC1)
         tcoul(2) = tcoul(2) + (TWC2 - TWC1)


      END DO  !end batch loop

      if (nVrs.gt.0) then
c --- backtransform fock matrix in full storage
         mode = 'tofull'
         add  = JRED.gt.JRED1
         Call switch_sto(irc,iLoc,ipFLT,ipFab,mode,add)
      endif

C --- free memory
      Call GetMem('VJ','Free','Real',ipVJ,nVec)
      Call GetMem('rsL','Free','Real',ipLab,LREAD)
      Call GetMem('rsF','Free','Real',ipFab,nRS)
      Call GetMem('rsD','Free','Real',ipDab,nRS)


999   Continue

      END DO   ! loop over red sets


      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky SCF timing from '//SECNAM
      Write(6,CFmt)'-----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COULOMB                          '
     &                           //'         ',tcoul(1),tcoul(2)
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
      ioff=0
      DO ISYM=1,NSYM
        ISFI= ioff + 1
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT('Coulomb Fmat',' ',FLT(ISFI),NB)
        END IF
        ioff = ioff + NB*(NB+1)/2
      END DO

      endif

#endif


      irc=0

      Return
      End



      SUBROUTINE switch_sto(irc,iLoc,ipXLT,ipXab,mode,add)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer  ISLT(8),cho_isao
      External cho_isao
      Logical  add
      Character*6 mode

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

************************************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************


c Offsets to symmetry block in the LT matrix
      ISLT(1)=0
      DO ISYM=2,NSYM
         ISLT(ISYM) = ISLT(ISYM-1)
     &              + NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
      END DO

**************************************************

      jSym = 1 ! only total symmetric density

      xf=0.0d0
      if (add) xf=1.0d0 !accumulate contributions

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

            Work(ipXab+jRab-1) = xf*Work(ipXab+jRab-1)
     &                         +    Work(kfrom)

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

            Work(kto) = xf*Work(kto)
     &                +    Work(ipXab+jRab-1)

         End Do  ! jRab loop

      Else

         write(6,*)'Wrong input parameter. mode = ',mode
         irc = 66
         Call abend()

      EndIf

      irc = 0

      Return
      End
