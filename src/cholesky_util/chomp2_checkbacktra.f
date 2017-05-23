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
* Copyright (C) 2008, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_CheckBackTra(iTyp,COcc,CVir,lU_AO)
C
C     Thomas Bondo Pedersen, Jan. 2008.
C
C     Purpose: check backtransformation of vectors (MO->AO).
C              A summary is printed at the end of this routine.
C
C     The check is simple, comparing the quantities
C
C        X(J) = sum_alpha,beta L(J;alpha,beta)
C        Y(J) = sum_ai L(ai,J)*P(a)*P(i)
C
C     where
C
C        P(i) = sum_alpha COcc(i,alpha)
C        P(a) = sum_alpha CVir(alpha,a)
C
C     The reported abs. min., abs. max, average, and RMS errors
C     are calculated from the vector
C
C        D(J) = X(J) - Y(J)
C
      Implicit None
      Integer iTyp
      Real*8  COcc(*), CVir(*)
      Integer lU_AO(*)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Real*8   dDot_
      external ddot_

      Real*8  Err(4,8)
      Real*8  AbsMinErr, AbsMaxErr, AvgErr, RMSErr
      Integer nMP2Vec_Tot

      Integer ip_POcc, l_POcc, ip_PVir, l_PVir, ip_Q, l_Q
      Integer ip_X, l_X, ip_Y, l_Y, ip_D, l_D, ip_V, l_V
      Integer iSym, iSyma, iSymi, iSymAl, iSymBe
      Integer i, a, Al, AlBe, nAlBe, na, J
      Integer kP, kP_, kC
      Integer iOpt, lTot, iAdr
      Integer nOcc_Max

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

C     Initializations.
C     ----------------

      nMP2Vec_Tot = 0

      AvgErr = 0.0d0
      RMSErr = 0.0d0

      nOcc_Max = nOcc(1)
      Do iSym = 2,nSym
         nOcc_Max = max(nOcc_Max,nOcc(iSym))
      End Do

      l_POcc = nOccT
      l_PVir = nVirT
      l_Q = nOcc_Max
      Call GetMem('POcc','Allo','Real',ip_POcc,l_POcc)
      Call GetMem('PVir','Allo','Real',ip_PVir,l_PVir)
      Call GetMem('Q','Allo','Real',ip_Q,l_Q)

C     P(i) = sum_alpha COcc(i,alpha)
C     ------------------------------

      Call Cho_dZero(Work(ip_POcc),l_POcc)
      Do iSymi = 1,nSym
         iSymAl = iSymi
         If (nOcc(iSymi).gt.0 .and. nBas(iSymAl).gt.0) Then
            kP = ip_POcc - 1 + iOcc(iSymi)
            Do Al = 1,nBas(iSymAl)
               kC = iT1AOT(iSymi,iSymAl) + nOcc(iSymi)*(Al-1)
               Do i = 1,nOcc(iSymi)
                  Work(kP+i) = Work(kP+i) + COcc(kC+i)
               End Do
            End Do
         End If
      End Do

C     P(a) = sum_alpha CVir(alpha,a)
C     ------------------------------

      Call Cho_dZero(Work(ip_PVir),l_PVir)
      Do iSyma = 1,nSym
         iSymAl = iSyma
         If (nVir(iSyma).gt.0 .and. nBas(iSymAl).gt.0) Then
            kP_ = ip_PVir - 1 + iVir(iSyma)
            Do a = 1,nVir(iSyma)
               kP = kP_ + a
               kC = iAOVir(iSymAl,iSyma) + nBas(iSymAl)*(a-1)
               Do Al = 1,nBas(iSyma)
                  Work(kP) = Work(kP) + CVir(kC+Al)
               End Do
            End Do
         End If
      End Do

C     Check each symmetry block.
C     --------------------------

      Do iSym = 1,nSym

         If (nMP2Vec(iSym) .gt. 0) Then

C           Allocation.
C           -----------

            l_X = nMP2Vec(iSym)
            l_Y = l_X
            Call GetMem('X','Allo','Real',ip_X,l_X)
            Call GetMem('Y','Allo','Real',ip_Y,l_Y)

C           Zero result arrays.
C           -------------------

            Call Cho_dZero(Work(ip_X),l_X)
            Call Cho_dZero(Work(ip_Y),l_Y)

C           X(J) = sum_alpha,beta L(J;alpha,beta)
C           -------------------------------------

            nAlBe = 0
            Do iSymBe = 1,nSym
               iSymAl = MulD2h(iSymBe,iSym)
               nAlBe = nAlBe + nBas(iSymAl)*nBas(iSymBe)
            End Do

            l_V = nMP2Vec(iSym)
            Call GetMem('VAO','Allo','Real',ip_V,l_V)
            Do AlBe = 1,nAlBe
               iOpt = 2
               lTot = nMP2Vec(iSym)
               iAdr = nMP2Vec(iSym)*(AlBe-1) + 1
               Call ddaFile(lU_AO(iSym),iOpt,Work(ip_V),lTot,iAdr)
              Call dAXPY_(nMP2Vec(iSym),1.0d0,Work(ip_V),1,Work(ip_X),1)
            End Do
            Call GetMem('VAO','Free','Real',ip_V,l_V)

C           Y(J) = sum_ai L(ai,J)*P(a)*P(i)
C           -------------------------------

            l_V = nT1Am(iSym)
            Call GetMem('VMO','Allo','Real',ip_V,l_V)
            iOpt = 1
            Call ChoMP2_OpenF(iOpt,iTyp,iSym)
            Do J = 1,nMP2Vec(iSym)
               iOpt = 2
               lTot = nT1Am(iSym)
               iAdr = nT1Am(iSym)*(J-1) + 1
               Call ddaFile(lUnit_F(iSym,iTyp),iOpt,Work(ip_V),lTot,
     &                      iAdr)
               Do iSymi = 1,nSym
                  iSyma = MulD2h(iSymi,iSym)
                  na = max(nVir(iSyma),1)
                  Call dGeMV_('T',nVir(iSyma),nOcc(iSymi),
     &                       1.0d0,Work(ip_V+iT1Am(iSyma,iSymi)),na,
     &                             Work(ip_PVir+iVir(iSyma)),1,
     &                       0.0d0,Work(ip_Q),1)
                  Work(ip_Y-1+J) = Work(ip_Y-1+J)
     &                           + dDot_(nOcc(iSymi),Work(ip_Q),1,
     &                                  Work(ip_POcc+iOcc(iSymi)),1)
               End Do
            End Do
            iOpt = 2
            Call ChoMP2_OpenF(iOpt,iTyp,iSym)
            Call GetMem('VMO','Free','Real',ip_V,l_V)

C           Calculate errors.
C           -----------------

            l_D = nMP2Vec(iSym)
            Call GetMem('D','Allo','Real',ip_D,l_D)
            Call dCopy_(nMP2Vec(iSym),Work(ip_X),1,Work(ip_D),1)
            Call dAXPY_(nMP2Vec(iSym),-1.0d0,Work(ip_Y),1,Work(ip_D),1)
            Err(1,iSym) = abs(Work(ip_D))
            Err(2,iSym) = abs(Work(ip_D))
            Err(3,iSym) = Work(ip_D)
            Do J = 1,nMP2Vec(iSym)-1
               Err(1,iSym) = min(Err(1,iSym),abs(Work(ip_D+J)))
               Err(2,iSym) = max(Err(2,iSym),abs(Work(ip_D+J)))
               Err(3,iSym) = Err(3,iSym) + Work(ip_D+J)
            End Do
            AvgErr = AvgErr + Err(3,iSym)
            Err(3,iSym) = Err(3,iSym)/dble(nMP2Vec(iSym))
            Err(4,iSym) = dDot_(nMP2Vec(iSym),Work(ip_D),1,Work(ip_D),1)
            RMSErr = RMSErr + Err(4,iSym)
            Err(4,iSym) = Err(4,iSym)/dble(nMP2Vec(iSym))
            Err(4,iSym) = sqrt(Err(4,iSym))
            Call GetMem('D','Free','Real',ip_D,l_D)

C           Deallocation.
C           -------------

            Call GetMem('Y','Free','Real',ip_Y,l_Y)
            Call GetMem('X','Free','Real',ip_X,l_X)

         Else

            Call Cho_dZero(Err(1,iSym),4)

         End If

         If (iSym .eq. 1) Then
            AbsMinErr = Err(1,iSym)
            AbsMaxErr = Err(2,iSym)
            nMP2Vec_Tot = max(nMP2Vec(iSym),0)
         Else
            AbsMinErr = min(AbsMinErr,Err(1,iSym))
            AbsMaxErr = max(AbsMaxErr,Err(2,iSym))
            nMP2Vec_Tot = nMP2Vec_Tot + max(nMP2Vec(iSym),0)
         End If

      End Do

      AvgErr = AvgErr/dble(nMP2Vec_Tot)
      RMSErr = RMSErr/dble(nMP2Vec_Tot)
      RMSErr = sqrt(RMSErr)

C     Report.
C     -------

      Call Cho_Head('MO Vector Backtransformation Check','=',80,6)
      Write(6,'(/,2X,A,A,/,2X,A,A)')
     & 'Symmetry  Min. Abs. Error  Max. Abs. Error    Average Error',
     & '        RMS Error',
     & '-----------------------------------------------------------',
     & '-----------------'
      Do iSym = 1,nSym
         Write(6,'(4X,I2,4X,1P,4(3X,D14.6))')
     &   iSym,(Err(i,iSym),i=1,4)
      End Do
      Write(6,'(2X,A,A)')
     & '-----------------------------------------------------------',
     & '-----------------'
      Write(6,'(2X,A,2X,1P,4(3X,D14.6))')
     & 'Total:',AbsMinErr,AbsMaxErr,AvgErr,RMSErr
      Write(6,'(2X,A,A,/)')
     & '-----------------------------------------------------------',
     & '-----------------'

C     Free memory.
C     ------------

      Call GetMem('Q','Free','Real',ip_Q,l_Q)
      Call GetMem('PVir','Free','Real',ip_PVir,l_PVir)
      Call GetMem('POcc','Free','Real',ip_POcc,l_POcc)

      End
