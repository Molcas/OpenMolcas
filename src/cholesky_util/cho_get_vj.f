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
*  CHO_get_VJ
*
*> @author F. Aquilante
*>
*> @details
*> Computes the Coulomb intermediate obtained
*> contracting the Cholesky vectors with the AO density matrix:
*>
*> \f[ \mathit{VJ}_k = \sum_{ab} L_{ab,k} D_{ab} \f]
*>
*> Reduced sets storage is used throughout the calculation.
*>
*> @note
*> Requires initialization of the Cholesky information.
*>
*> @param[out] irc     Return code
*> @param[in]  ipDLT   Pointer to density matrix, stored as packed LT, consecutive symmetry blocks
*> @param[out] VJ      Coulomb intermediate
*> @param[in]  Mvec    Number of Cholesky vectors (min. length of \p VJ)
*> @param[in]  timings Switch on/off timings printout
************************************************************************
      SUBROUTINE CHO_get_VJ(irc,ipDLT,VJ,Mvec,timings)
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
      Logical timings
      Integer ipDLT, Mvec
      Real*8  VJ(Mvec)
      Real*8  tread(2), tcalv(2)

      Logical  add
      Character*6  mode
      Character*10  SECNAM

      parameter (SECNAM = 'CHO_get_VJ')
      parameter (zero = 0.0d0, one = 1.0d0)

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      JSYM=1

      If (Mvec.ne.NumCho(jSym)) Then
         Write(6,*)SECNAM//': Mvec must be equal to NumCho(1) !!'
         irc = 77
         Return
      EndIf

      If (Mvec.lt.1) Then
         irc = 0
         Return  ! could happen in a parallel run
      EndIf

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero  !time for reading the vectors
         tcalv(i) = zero  !time for computing VJ
      end do

      iLoc = 3 ! use scratch location in reduced index arrays

C ------------------------------------------------------------------

      JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
      JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

      Do JRED=JRED1,JRED2

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
       Call FZero(Work(ipDab),nRS)

       Call GetMem('MaxM','Max','Real',KDUM,LWORK)

       nVec  = Min(LWORK/nRS,nVrs)

       If (nVec.lt.1) Then
         WRITE(6,*) SECNAM//': Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'min. mem. need= ',nRS+1
         irc = 33
         CALL Abend()
         nBatch = -9669  ! dummy assignment
       End If

       LREAD = nRS*nVec

       Call GetMem('rsL','Allo','Real',ipLab,LREAD)

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

         CALL CWTIME(TCI1,TWI1)

C --- Transform the density to reduced storage
         mode = 'toreds'
         add  = .false.
         Call switch_sto(irc,iLoc,ipDLT,ipDab,mode,add)

C ------------------------------------------------------------
C --- V{#J} <- V{#J} + sum_ab  L(ab,{#J}) * D(ab)
C=============================================================

         CALL DGEMV_('T',nRS,JNUM,
     &              One,Work(ipLab),nRS,
     &                  Work(ipDab),1,
     &              zero,VJ(jVec),1)


         CALL CWTIME(TCI2,TWI2)
         tcalv(1) = tcalv(1) + (TCI2 - TCI1)
         tcalv(2) = tcalv(2) + (TWI2 - TWI1)


       END DO  !end batch loop

C --- free memory
       Call GetMem('rsL','Free','Real',ipLab,LREAD)
       Call GetMem('rsD','Free','Real',ipDab,nRS)


999    Continue

      END DO   ! loop over red sets

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

*
*---- Write out timing information
      if(timings)then

      Write(6,*)
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)'Timing from ', SECNAM,'            CPU      WALL '
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                  '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COMPUTE V^J                   '
     &                           //'         ',tcalv(1),tcalv(2)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                         '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      irc=0

      Return
      End
