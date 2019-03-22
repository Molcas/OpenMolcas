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
* Copyright (C) 2007, Francesco Aquilante                              *
************************************************************************
      SubRoutine Cho_SOSmp2_Energy(irc,EMP2,EOcc,EVir,Delete)
C
C     Francesco Aquilante, May 2007.
C
C     Purpose: compute "Scaled Opposite-Spin" MP2 energy correction from
C              MO Cholesky vectors of the matrix M(ai,bj)=(ai|bj)^2.
C
      Implicit Real*8 (a-h,o-z)
      Real*8  EMP2
      Real*8  EOcc(*), EVir(*)
      Integer irc
      Logical Delete
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"

      Character*6  ThisNm
      Character*17 SecNam
      Parameter (SecNam = 'Cho_SOSmp2_Energy', ThisNm = 'Energy')
      Parameter (zero = 0.0d0,  one = 1.0d0, two = 2.0d0)

      Integer ipWrk, lWrk, iiSoff(8), iaSoff(8)
      Integer nEnrVec(8)
*****************************************************************
      MulD2h(i,j)=iEor(i-1,j-1)+1
*****************************************************************

      Call qEnter(ThisNm)
      irc = 0

      iTyp = 2
      Call iCopy(nSym,nMP2Vec,1,nEnrVec,1)

C     Initialize SOS-MP2 energy correction.
C     -------------------------------------

      EMP2 = 0.0D0

C     Some offsets
C     ------------
      nIt=nOcc(1)
      nAt=nVir(1)
      iiSoff(1)=0
      iaSoff(1)=0
      Do iSym=2,nSym
         iiSoff(iSym)=nIt
         iaSoff(iSym)=nAt
         nIt = nIt + nOcc(iSym)
         nAt = nAt + nVir(iSym)
      End Do
      MaxNVec = nIt*nAt

      Call GetMem('Ea-Ei','Allo','Real',ip_W,MaxNVec)
      Call GetMem('iD_bj','Allo','Inte',ID_bj,MaxNVec)

C     Loop over Cholesky vector symmetries.
C     -------------------------------------

      Do jSym = 1,nSym

         Nai = nT1am(jSym)
         If (Nai.gt.0 .and. nEnrVec(jSym).gt.0) Then

            nOV=0
            Do iiSym=1,nSym
               iaSym=MulD2h(iiSym,jSym)
               Do ii=1,nOcc(iiSym)
                  iiT=iiSoff(iiSym)+ii
                  iaS=nOV+nVir(iaSym)*(ii-1)
                  Do ia=1,nVir(iaSym)
                     iaT=iaSoff(iaSym)+ia
                     iaiS=iaS+ia-1
                     Work(ip_W+iaiS)=EVir(iaT)-EOcc(iiT)
                  End Do
               End Do
               nOV=nOV+nVir(iaSym)*nOcc(iiSym) ! ... = Nai
            End Do

C           Cholesky decompsition of the Orbital Energy
C           Denominators (OED)
C           -------------------------------------------
            Call CHO_GET_ORD_bj(nOV,MaxNVec,OED_Thr,Work(ip_W),
     &                          iWork(ID_bj),NKVec,Dmax)

            If (Verbose .or. NKVec.lt.1) Then
               Write(6,'(A)')
     &         '---------------------------------------'
               Write(6,'(A,I2,A)')
     &         'Orbital energy denominators CD (sym=',jSym,')'
               Write(6,'(A)')
     &         '---------------------------------------'
               Write(6,'(1X,A,I3,A,I9,A,1P,D25.16)')
     &         'Number of vectors needed: ',NKVec,
     &         '   ( nAocc x nAvir : ',nOV,' ), max residual:',
     &         Dmax
               Call xFlush(6)
            EndIf

            If (NKVec.gt.0) Then

               Call GetMem('Yai_k','Allo','Real',ip_Y,nOV*NKVec)
!              init to one the 1st col
               call dcopy_(nOV,one,0,Work(ip_Y),1)

               Call CHO_GET_OED_cd(.true.,nOV,Work(ip_W),NKVec,
     &                             iWork(ID_bj),1,Work(ip_Y),Work(ip_Y))

               Call GetMem('GetMax','Max ','Real',ipWrk,lWrk)
               Call GetMem('GetMax','Allo','Real',ipWrk,lWrk)

C              Set up batch over Cholesky vectors.
C              -----------------------------------

               nVec = min(lWrk/(Nai+NKVec),nEnrVec(jSym))
               If (nVec .lt. 1) Then
                  Call ChoMP2_Quit(SecNam,'Insufficient memory',
     &                             'Batch setup')
               End If
               nBat = (nEnrVec(jSym)-1)/nVec + 1

               kRead = ipWrk + NKVec*nVec

C              Open Cholesky vector files.
C              ---------------------------

               Call ChoMP2_OpenF(1,iTyp,jSym)

C              Start vector batch loop.
C              ------------------------

               Do iBat = 1,nBat

                  If (iBat .eq. nBat) Then
                     NumVec = nEnrVec(jSym) - nVec*(nBat-1)
                  Else
                     NumVec = nVec
                  End If
                  jVec = nVec*(iBat-1) + 1

C                 Read vectors
C                 ------------
                  lTot = nT1am(jSym)*NumVec
                  iAdr = nT1am(jSym)*(jVec-1) + 1
                  Call ddaFile(lUnit_F(jSym,iTyp),2,
     &                         Work(kRead),lTot,iAdr)

C                 Compute   E(k,J) = sum_ai Y(ai,k) * R(ai,J)
C                 --------------------------------------------
                  Call DGEMM_('T','N',NKVec,NumVec,Nai,
     &                       one,Work(ip_Y),Nai,Work(kRead),Nai,
     &                       zero,Work(ipWrk),NKVec)

C                 Compute (unscaled) SOS-MP2 energy
C                 -----------------------------------

                  EMP2 = EMP2 + ddot_(NKVec*NumVec,Work(ipWrk),1,
     &                                            Work(ipWrk),1)

               End Do ! Cholesky vector batch

C              Close Cholesky vector files.
C              ----------------------------
               Call ChoMP2_OpenF(2,iTyp,jSym)

               Call GetMem('GetMax','Free','Real',ipWrk,lWrk)
               Call GetMem('Yai_k','Free','Real',ip_Y,nOV*NKVec)

            End If

         End If

      End Do


      Call GetMem('iD_bj','Free','Inte',ID_bj,MaxNVec)
      Call GetMem('Ea-Ei','Free','Real',ip_W,MaxNVec)

C     If requested, delete vector files.
C     ----------------------------------

      If (Delete) Then
         Do jSym = 1,nSym
            Call ChoMP2_OpenF(1,iTyp,jSym)
            Call ChoMP2_OpenF(3,iTyp,jSym)
         End Do
      End If

C     Change sign and use proper factor on energy.
C     --------------------------------------------

C-tbp, December 2012: removed factor 2
C-tbp  (wrong result for two-electron systems with factor 2)
c-tbp EMP2 = -two*EMP2
      EMP2 = -EMP2


      Call qExit(ThisNm)

      Return
      End
