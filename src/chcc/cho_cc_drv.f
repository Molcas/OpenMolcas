************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CHO_CC_drv(rc,CMO)

************************************************************************
C
C      a,b,g,d:  AO-index
C      p,q,r,s:  MO-indices belonging to (probably frozen excluded ?)
C
************************************************************************
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use Data_Structures, only: DSBA_Type, SBA_Type
      use Data_Structures, only: Allocate_DT, Deallocate_DT
      Implicit Real*8 (a-h,o-z)

      Integer   rc
      Type (DSBA_Type) CMO

      Real*8    tread(2),tmotr1(2),tmotr2(2)
      Logical   DoRead
      Integer   nPorb(8)

      Character*50 CFmt
      Character(LEN=10), Parameter:: SECNAM = 'CHO_CC_drv'

#include "real.fh"
#include "chotime.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Type (SBA_Type) Laq(1)
      Real*8, Allocatable:: Lrs(:,:)
      Real*8, Allocatable,Target:: Lpq(:)
      Real*8, Pointer:: pLpq(:,:,:)=>Null()

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
************************************************************************

      LunChVF = 80
      LunChVF = isfreeunit(LunChVF)
      call DaName_mf_wa (LunChVF,'CD1tmp')
      idisk=1
      DoRead  = .false.
      IREDC = -1  ! unknown reduced set in core

      iSwap = 0  ! Lpb,J are returned by cho_x_getVtra
      kMOs = 1
      nMOs = 1


        CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        do i=1,2            ! 1 --> CPU   2 --> Wall
           tread(i) = zero   !time read/write vectors
           tmotr1(i) = zero  !time 1st MO half-transf.
           tmotr2(i) = zero  !time 2nd MO half-transf.
        end do

c --- Define MOs used in CC
c -----------------------------------
        do i=1,nSym
           nPorb(i) = SIZE(CMO%SB(i)%A2,1)
        end do


C ==================================================================

c --- Various offsets & pointers
c ------------------------------

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
c
c
      DO jSym=1,nSym

         If (NumCho(jSym).lt.1) GOTO 1000

C -------------------------------------------------------------


C ****************     MEMORY MANAGEMENT SECTION    *****************
C ------------------------------------------------------------------
C --- compute memory needed to store at least 1 vector of JSYM
C --- and do all the subsequent calculations
C ------------------------------------------------------------------
         mTvec  = 0  ! mem for storing half-transformed vec Laq,J
         mTTvec = 0  ! mem for storing transformed vec Lpq,J

         do l=1,nSym
            k=Muld2h(l,JSYM)
            mTvec = mTvec + nPorb(l)*nBas(k)
            mTTvec = Max(mTTvec,nPorb(l)*nPorb(k))
         end do

         mvec = mTvec + mTTvec

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
     &                         ' rc= ',irc
              call abend()
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            Call mma_maxDBLE(LWORK)

            nVec  = Min(LWORK/(nRS+mvec),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) 'LWORK= ',LWORK
               WRITE(6,*) 'min. mem. need= ',nRS+mTvec
               WRITE(6,*) 'reading ',nRS,' and transforming to ',mvec
               WRITE(6,*) 'of jsym= ',jsym,' and JRED= ',JRED
               rc = 33
               CALL Abend()
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call mma_allocate(Lrs,nRS,nVec,Label='Lrs')
            Call Allocate_DT(Laq(1),nPorb,nBas,nVec,jSym,nSym,iSwap)
            Call mma_allocate(Lpq,mTTVec*nVec,Label='Lpq')

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

               If (NUMV.le.0 .or.NUMV.ne.JNUM) then
                  rc=77
                  RETURN
               End If

               CALL CWTIME(TCR2,TWR2)
               tread(1) = tread(1) + (TCR2 - TCR1)
               tread(2) = tread(2) + (TWR2 - TWR1)

C --------------------------------------------------------------------
C --- First half MO transformation  Lpb,J = sum_a  C(p,a) * Lab,J
C --------------------------------------------------------------------

               CALL CWTIME(TCM1,TWM1)

               CALL CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,
     &                           JSYM,iSwap,IREDC,nMOs,kMOs,[CMO],
     &                           Laq(1),DoRead)

               if (irc.ne.0) then
                  rc = irc
                  RETURN
               endif

               CALL CWTIME(TCM2,TWM2)
               tmotr1(1) = tmotr1(1) + (TCM2 - TCM1)
               tmotr1(2) = tmotr1(2) + (TWM2 - TWM1)

C --------------------------------------------------------------------
C --- 2nd half of MO transformation  Lpq,J = sum_b  Lpb,J * C(q,b)
C --------------------------------------------------------------------
               Do iSymb=1,nSym

                  iSymp = MulD2h(JSYM,iSymb)
                  NAp = nPorb(iSymp)
                  NAq = nPorb(iSymb) ! iSymb=iSymq
                  iS = 1
                  iE = NAp * NAq * JNUM

                  pLpq(1:NAp,1:NAq,1:JNUM) => Lpq(iS:iE)

                  CALL CWTIME(TCM3,TWM3)

                  If(NAp*NAq.ne.0)Then

                    Do JVC=1,JNUM

                      CALL DGEMM_('N','T',NAp,NAq,nBas(iSymb),
     &                           One,Laq(1)%SB(iSymb)%A3(:,:,JVC),NAp,
     &                               CMO%SB(iSymb)%A2,NAq,
     &                          Zero,pLpq(:,:,JVC),NAp)

                      End Do

                  EndIf

                  CALL CWTIME(TCM4,TWM4)
                  tmotr2(1) = tmotr2(1) + (TCM4 - TCM3)
                  tmotr2(2) = tmotr2(2) + (TWM4 - TWM3)

C if u need to compute fock matrix elements this should be done probably here
C     I can help you with that

                  CALL CWTIME(TCR3,TWR3)
C --- WRITE transformed vectors to disk (each Jsym on a separate file!)
c
                  call ddafile (LunChVF,1,Lpq,NAp*NAq*JNUM,idisk)
c
C --- remember that this is inside a batch over J, the vector index

                  CALL CWTIME(TCR4,TWR4)
                  tread(1) = tread(1) + (TCR4 - TCR3)
                  tread(2) = tread(2) + (TWR4 - TWR3)

                  pLpq => Null()

               End Do

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop

C --- free memory
            Call mma_deallocate(Lpq)
            Call Deallocate_DT(Laq(1))
            Call mma_deallocate(Lrs)

999         CONTINUE

         END DO   ! loop over red sets


1000     CONTINUE

      END DO   !loop over JSYM
      call daclos(LunChVF)

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1
*
*---- Write out timing information
      if(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky-CC timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'MO transf. Cholesky vectors     CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

         Write(6,'(2x,A26,2f10.2)')'READ/WRITE VECTORS               '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'1st half-transf.                 '
     &                           //'         ',tmotr1(1),tmotr1(2)
         Write(6,'(2x,A26,2f10.2)')'2nd half-transf.                 '
     &                           //'         ',tmotr2(1),tmotr2(2)
         Write(6,*)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      rc  = 0

      Return
      END
