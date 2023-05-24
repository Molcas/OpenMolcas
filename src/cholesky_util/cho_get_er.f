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
!  CHO_get_ER
!
!> @brief
!>   Compute the Edmiston--Ruedenberg functional for a given set of occupied MOs
!> @author F. Aquilante
!>
!> @details
!> Computes the Edmiston--Ruedenberg functional
!>
!> \f[ W = \sum_i \mathit{ER}[i] = \sum_i (ii|ii) \f]
!>
!> for a given set of occupied MOs.
!>
!> The functional orbital components \p ER(i) are
!> computed by using the Cholesky representation \f$ L_{ab,J} \f$
!> of the AO two-electron integrals, namely
!>
!> \f[ \mathit{ER}[i] = \sum_J V[i]_J V[i]_J \f]
!>
!> where
!>
!> \f[ V[i]_J = \sum_{ab} D[i]_{ab} L_{ab,J} \f]
!>
!> and
!>
!> \f[ D[i]_{a,b} = C[i]_a C[i]_b \f]
!>
!> @note
!> Requires initialization of the Cholesky information.
!>
!> @param[out] irc     return code
!> @param[in]  CMO     MOs matrix, stored as \p C(a,k)
!> @param[in]  nOcc    number of occupied orbitals in each symmetry
!> @param[out] ER      orbital components of the ER functional}
!> @param[out] W       value of the Edmiston--Ruedenberg functional
!> @param[in]  timings switch on/off timings printout
!***********************************************************************
      SUBROUTINE CHO_get_ER(irc,CMO,nOcc,ER,W,timings)
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      use Constants
      Implicit Real*8 (a-h,o-z)
      Integer irc
      Integer nOcc(*)
      Real*8  CMO(*),ER(*),W
      Logical timings

      Integer iOcc(8),isMO(8)
      Real*8  tread(2),tintg(2)
      Character(LEN=10), Parameter:: SECNAM = 'CHO_get_ER'

#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      Real*8, Allocatable:: DLT(:), Lab(:), Dab(:), VJ(:)

      JSYM=1
      If (NumCho(JSYM).lt.1) Then
         Write(6,*)SECNAM//'No total symmetric vectors present'
         irc = 77
         Return
      EndIf

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero  !time for reading the vectors
         tintg(i) = zero  !time for computing the functional
      end do

! --- compute some offsets and other quantities
      MaxB=nBas(1)
      iOcc(1)=0
      isMO(1)=0
      Do kSym=2,nSym
         iOcc(kSym) = iOcc(kSym-1) + nOcc(kSym-1)
         isMO(kSym) = isMO(kSym-1) + nBas(kSym-1)**2
         MaxB = Max(MaxB,nBas(kSym))
      End Do

      MaxBB = MaxB*(MaxB+1)/2
      nOccT = iOcc(nSym) + nOcc(nSym)

      W = zero ! initialization of the ER-functional value
      Call Fzero(ER(1),nOccT) ! and its orbital components

      Call mma_allocate(DLT,MaxBB,Label='DLT')
      DLT(:)=Zero

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

      Call mma_allocate(Dab,nRS,Label='Dab')

      Call mma_maxDBLE(LWORK)

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

      Call mma_allocate(Lab,LREAD,Label='Lab')
      Call mma_allocate(VJ,nVec,Label='VJ')

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

         CALL CHO_VECRD(Lab,LREAD,JVEC,IVEC2,JSYM,NUMV,JRED,MUSED)

         If (NUMV.le.0 .or. NUMV.ne.JNUM) then
            irc=77
            RETURN
         End If

         CALL CWTIME(TCR2,TWR2)
         tread(1) = tread(1) + (TCR2 - TCR1)
         tread(2) = tread(2) + (TWR2 - TWR1)

         CALL CWTIME(TCI1,TWI1)

         Do kSym=1,nSym

            Do ik=1,nOcc(kSym)

! --- Compute the i-th orbital component D[i](a,b) of the density
! --- D[i](a,b) is stored as upper-triangular
               Do ib=1,nBas(kSym)

                  ipb = isMO(kSym) + nBas(kSym)*(ik-1) + ib

                  Do ia=1,ib-1

                     ipa = isMO(kSym) + nBas(kSym)*(ik-1) + ia
                     ipab = ib*(ib-1)/2 + ia
                     DLT(ipab) = CMO(ipa)*CMO(ipb)

                  End Do

                  ipab = ib*(ib+1)/2
                  DLT(ipab) = half*CMO(ipb)**2  !diagonal scaled

               End Do

! --- Transform the density to reduced storage
               Call switch_density(iLoc,DLT,Dab,kSym)

!  ( r .ge. s )
! ------------------------------------------------------------
! --- V[i]{#J} <- V[i]{#J} + 2 * sum_rs  L(rs,{#J}) * D[i](rs)
!=============================================================

               CALL DGEMV_('T',nRS,JNUM,
     &                    TWO,Lab,nRS,
     &                    Dab,1,ZERO,VJ,1)


! ----------------------------------------------------------
! --- ER[i] <- ER[i]  +  sum_J V[i](J)^2
!===========================================================
               Do jv=1,JNUM

                  ER(iOcc(kSym)+ik) = ER(iOcc(kSym)+ik) + VJ(jv)**2
               End Do


            End Do


         End Do

         CALL CWTIME(TCI2,TWI2)
         tintg(1) = tintg(1) + (TCI2 - TCI1)
         tintg(2) = tintg(2) + (TWI2 - TWI1)


      END DO  !end batch loop

! --- free memory
      Call mma_deallocate(VJ)
      Call mma_deallocate(Lab)
      Call mma_deallocate(Dab)


999   Continue

      END DO   ! loop over red sets

      Call mma_deallocate(DLT)

! --- Sync components
      Call GAdGOp(ER,nOccT,'+')

! --- Compute the ER-functional from its orbital components
      W = Zero
      Do ik=1,nOccT
         W = W + ER(ik)
      End Do



      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

!
!---- Write out timing information
      if(timings)then

      Write(6,*)
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)'Timing from ', SECNAM,'            CPU      WALL '
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'COMPUTE (ii|ii)                  '
     &                           //'         ',tintg(1),tintg(2)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


      irc=0

      Return
      End



      SUBROUTINE switch_density(iLoc,XLT,Xab,kSym)
      use ChoArr, only: iRS2F
      use ChoSwp, only: IndRed
      Implicit Real*8 (a-h,o-z)
      Integer, External:: cho_isao
      Integer iLoc, kSym
      Real*8 XLT(*), Xab(*)

#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"


!***********************************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
!***********************************************************************

      jSym = 1 ! only total symmetric density

      Do jRab=1,nnBstR(jSym,iLoc)

         kRab = iiBstr(jSym,iLoc) + jRab
         iRab = IndRed(kRab,iLoc)

         iag   = iRS2F(1,iRab)  !global address
         ibg   = iRS2F(2,iRab)

         iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

         xf = dble(1-min(1,abs(iSyma-kSym)))

         ias   = iag - ibas(iSyma)  !address within that symm block
         ibs   = ibg - ibas(iSyma)
         iab   = iTri(ias,ibs)

         Xab(jRab) = xf*XLT(iab)

      End Do  ! jRab loop

      Return
      End
