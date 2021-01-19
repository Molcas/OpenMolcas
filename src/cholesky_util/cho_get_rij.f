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
*  CHO_get_Rij
*
*> @author F. Aquilante
*>
*> @details
*> Computes the \f$ R \f$ matrix used for the
*> maximization of Edmiston--Ruedenberg functional
*>
*> \f[ \mathit{ER} = \sum_i (ii|ii) = \mathrm{Tr}(R) \f]
*>
*> for a given set of MOs.
*>
*> \f$ R \f$ is defined from the two-electron integrals
*> computed from the MO-transformed Cholesky vectors
*>
*> \f[ R_{ij} = (ij|jj) = \sum_K L_{ij,K} L_{jj,K} \f]
*>
*> and the condition for the maximization of the ER-functional
*> is given by
*>
*> \f[ \mathrm{grad}(\mathit{ER})_{ij} = 4(R_{ij} - R_{ji}) = 0 \quad (\forall i,j) \f]
*>
*> @note
*> Requires initialization of the Cholesky information.
*>
*> @param[out]    irc     Return code
*> @param[in]     ipMO    Pointers to each symmetry block of the MO matrix, stored as \p C(k,a)
*> @param[in]     nOcc    Number of orbitals to be localized in each symmetry
*> @param[in,out] Rij     \p nOcc &times; \p nOcc symmetry blocked matrix \f$  R_{ij} = (ij|jj) \f$
*> @param[in]     timings Switch on/off timings printout
************************************************************************
      SUBROUTINE CHO_get_Rij(irc,ipMO,nOcc,Rij,timings)
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
      Logical timings,DoRead
      Integer nOcc(*),iOcc(8),iOcs(8),ipLib(8),iSkip(8),ipMO(*)
      Real*8  Rij(*)
      Real*8  tread(2),tintg(2),tmotr(2)
      Character*11  SECNAM

      parameter (SECNAM = 'CHO_get_Rij')
      parameter (zero = 0.0d0, one = 1.0d0, DoRead = .false.)

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      IREDC = -1

      JSYM=1
      If (NumCho(JSYM).lt.1) Then
         Write(6,*)SECNAM//': No total symmetric vectors present'
         irc = 77
         Return
      EndIf

      Do kS=1,nSym
         If ( nOcc(kS).gt.nBas(kS) .or. nOcc(kS).lt.0) Then
          Write(6,*)SECNAM//': Wrong nOcc in symmetry ',kS
          Write(6,*)'nOcc(',kS,')= ',nOcc(kS)
          Write(6,*)'nBas(',kS,')= ',nBas(kS)
          irc = 79
          Return
         EndIf
      End Do

      CALL CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

      do i=1,2            ! 1 --> CPU   2 --> Wall
         tread(i) = zero  !time for reading the vectors
         tmotr(i) = zero  !time for MO transformation of the vectors
         tintg(i) = zero  !time for computing the functional
      end do

C --- compute some offsets and other quantities
      iOcc(1)=0
      iOcs(1)=0
      Mocc=nOcc(1)
      Do kSym=2,nSym
         iOcc(kSym) = iOcc(kSym-1) + nOcc(kSym-1)
         iOcs(kSym) = iOcs(kSym-1) + nOcc(kSym-1)**2
         Mocc = Max(Mocc,nOcc(kSym))
      End Do

      nOcs = iOcs(nSym) + nOcc(nSym)**2

      Call FZero(Rij(1),nOcs)

      Do kS=1,nSym
         iSkip(kS) = Min(nOcc(kS),1) ! initialize skipping flags
      End Do

C --- Memory need for 1 of the half-transformed vectors : L(aj)
      Maj=0
      Do iSyma=1,nSym
         Maj = Maj + nBas(iSyma)*nOcc(iSyma)
      End Do

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
       Write(6,*)SECNAM//': cho_X_setred non-zero return code. rc= ',irc
       call abend()
      endif

      IREDC=JRED

      nRS = nDimRS(JSYM,JRED)

      Mneed = Max(nRS, Mocc**2 + 1) ! mem. for Lab and (Lij + Ljj)

      Call GetMem('MaxM','Max','Real',KDUM,LWORK)

      nVec  = Min(LWORK/(Maj+Mneed),nVrs)

      If (nVec.lt.1) Then
         WRITE(6,*) SECNAM//': Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'min. mem. need= ',Maj+Mneed
         irc = 33
         CALL Abend()
         nBatch = -9999  ! dummy assignment
      End If

      LREAD = nRS*nVec

      Call GetMem('rsL','Allo','Real',ipLab,Mneed*nVec)
      Call GetMem('Ltr','Allo','Real',ipLtr,Maj*nVec)

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

C --- First half-transformation of the vectors : Lab,J --> LiJ,b
C --------------------------------------------------------------
         iSwap = 2  ! LiK,b are returned

         kMOs = 1
         nMOs = 1

         lChoT=0
         Do kSym=1,nSym

            ipLib(kSym) = ipLtr + lChot
            lChot = lChot + nOcc(kSym)*JNUM*nBas(kSym)

         End Do


         CALL CWTIME(TCT1,TWT1)

         CALL CHO_X_getVtra(irc,Work(ipLab),LREAD,jVEC,JNUM,
     &                         JSYM,iSwap,IREDC,nMOs,kMOs,ipMO,nOcc,
     &                         ipLib,iSkip,DoRead)

               if (irc.ne.0) then
                   RETURN
               endif


         CALL CWTIME(TCT2,TWT2)
         tmotr(1) = tmotr(1) + (TCT2 - TCT1)
         tmotr(2) = tmotr(2) + (TWT2 - TWT1)


         ipLij = ipLab  ! re-use mem used for reading from red-sets

         Do kSym=1,nSym

            If (iSkip(kSym) .ne. 0) Then

               ipLjj = ipLij + nOcc(kSym)**2 * JNUM

               CALL CWTIME(TCT1,TWT1)
C ---------------------------------------------------------------------
C --- Second half-transformation  L(iK,j) = sum_b  L(iK,b) * C(j,b)
C ---------------------------------------------------------------------

              CALL DGEMM_('N','T',nOcc(kSym)*JNUM,nOcc(kSym),nBas(kSym),
     &                            One,Work(ipLib(kSym)),nOcc(kSym)*JNUM,
     &                                Work(ipMO(kSym)),nOcc(kSym),
     &                           Zero,Work(ipLij),nOcc(kSym)*JNUM)


               CALL CWTIME(TCT2,TWT2)
               tmotr(1) = tmotr(1) + (TCT2 - TCT1)
               tmotr(2) = tmotr(2) + (TWT2 - TWT1)


               CALL CWTIME(TCI1,TWI1)

               Do lj=1,nOcc(kSym)

                  Do jv=1,JNUM

                     jfrom = ipLij + nOcc(kSym)*JNUM*(lj-1)
     &                     + nOcc(kSym)*(jv-1) + lj - 1
                     jto = ipLjj + jv - 1
                     Work(jto) = Work(jfrom)

                  End Do

C --------------------------------------------------------------------
C --- Compute   R(i,j) = sum_K  L(i,K)[j] * L(K)[j]
C --------------------------------------------------------------------

                  jpR = iOcs(kSym) + nOcc(kSym)*(lj-1) + 1
                  ipLik = ipLij + nOcc(kSym)*JNUM*(lj-1)

                  CALL DGEMV_('N',nOcc(kSym),JNUM,
     &                       ONE,Work(ipLik),nOcc(kSym),
     &                           Work(ipLjj),1,ONE,Rij(jpR),1)


               End Do

               CALL CWTIME(TCI2,TWI2)
               tintg(1) = tintg(1) + (TCI2 - TCI1)
               tintg(2) = tintg(2) + (TWI2 - TWI1)

            EndIf

         End Do



      END DO  !end batch loop

C --- free memory
      Call GetMem('Ltr','Free','Real',ipLtr,Maj*nVec)
      Call GetMem('rsL','Free','Real',ipLab,Mneed*nVec)


999   Continue

      END DO   ! loop over red sets

C --- sync Rij

      Call GAdGOp(Rij,nOcs,'+')

      CALL CWTIME(TOTCPU2,TOTWALL2)
      TOTCPU = TOTCPU2 - TOTCPU1
      TOTWALL= TOTWALL2 - TOTWALL1

*
*---- Write out timing information
      if(timings)then

      Write(6,*)
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)'Timing from ', SECNAM,'           CPU      WALL  '
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
         Write(6,'(2x,A26,2f10.2)')'READ VECTORS                     '
     &                           //'         ',tread(1),tread(2)
         Write(6,'(2x,A26,2f10.2)')'TRANSFORM VECTORS                '
     &                           //'         ',tmotr(1),tmotr(2)
         Write(6,'(2x,A26,2f10.2)')'COMPUTE Rij = (ij|jj)            '
     &                           //'         ',tintg(1),tintg(2)
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif


      irc=0

      Return
      End
