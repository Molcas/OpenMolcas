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
      use Constants
      Implicit None
      Integer iTyp
      Real*8  COcc(*), CVir(*)
      Integer lU_AO(*)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "stdalloc.fh"

      Real*8, External:: ddot_

      Real*8  Err(4,8)
      Real*8  AbsMinErr, AbsMaxErr, AvgErr, RMSErr
      Integer nMP2Vec_Tot

      Integer iSym, iSyma, iSymi, iSymAl, iSymBe
      Integer i, a, Al, AlBe, nAlBe, na, J
      Integer kP, kP_, kC
      Integer iOpt, lTot, iAdr
      Integer nOcc_Max

      Real*8, Allocatable:: POcc(:), PVir(:), Q(:), X(:), Y(:), D(:),
     &                      V(:)

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

C     Initializations.
C     ----------------

      nMP2Vec_Tot = 0

      AvgErr = Zero
      RMSErr = Zero

      nOcc_Max = nOcc(1)
      Do iSym = 2,nSym
         nOcc_Max = max(nOcc_Max,nOcc(iSym))
      End Do

      Call mma_allocate(POcc,nOccT,Label='POcc')
      Call mma_allocate(PVir,nVirT,Label='PVir')
      Call mma_allocate(Q,nOcc_Max,Label='Q')

C     P(i) = sum_alpha COcc(i,alpha)
C     ------------------------------

      POcc(:)=Zero
      Do iSymi = 1,nSym
         iSymAl = iSymi
         If (nOcc(iSymi).gt.0 .and. nBas(iSymAl).gt.0) Then
            kP = iOcc(iSymi)
            Do Al = 1,nBas(iSymAl)
               kC = iT1AOT(iSymi,iSymAl) + nOcc(iSymi)*(Al-1)
               Do i = 1,nOcc(iSymi)
                  POcc(kP+i) = POcc(kP+i) + COcc(kC+i)
               End Do
            End Do
         End If
      End Do

C     P(a) = sum_alpha CVir(alpha,a)
C     ------------------------------

      PVir(:)=Zero
      Do iSyma = 1,nSym
         iSymAl = iSyma
         If (nVir(iSyma).gt.0 .and. nBas(iSymAl).gt.0) Then
            kP_ = iVir(iSyma)
            Do a = 1,nVir(iSyma)
               kP = kP_ + a
               kC = iAOVir(iSymAl,iSyma) + nBas(iSymAl)*(a-1)
               Do Al = 1,nBas(iSyma)
                  PVir(kP) = PVir(kP) + CVir(kC+Al)
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

            Call mma_allocate(X,nMP2Vec(iSym),Label='X')
            Call mma_allocate(Y,nMP2Vec(iSym),Label='Y')

C           Zero result arrays.
C           -------------------

            X(:)=Zero
            Y(:)=Zero

C           X(J) = sum_alpha,beta L(J;alpha,beta)
C           -------------------------------------

            nAlBe = 0
            Do iSymBe = 1,nSym
               iSymAl = MulD2h(iSymBe,iSym)
               nAlBe = nAlBe + nBas(iSymAl)*nBas(iSymBe)
            End Do

            Call mma_allocate(V,nMP2Vec(iSym),Label='V')
            Do AlBe = 1,nAlBe
               iOpt = 2
               lTot = nMP2Vec(iSym)
               iAdr = nMP2Vec(iSym)*(AlBe-1) + 1
               Call ddaFile(lU_AO(iSym),iOpt,V,lTot,iAdr)
              Call dAXPY_(nMP2Vec(iSym),1.0d0,V,1,X,1)
            End Do
            Call mma_deallocate(V)

C           Y(J) = sum_ai L(ai,J)*P(a)*P(i)
C           -------------------------------

            Call mma_allocate(V,nT1Am(iSym),Label='V')
            iOpt = 1
            Call ChoMP2_OpenF(iOpt,iTyp,iSym)
            Do J = 1,nMP2Vec(iSym)
               iOpt = 2
               lTot = nT1Am(iSym)
               iAdr = nT1Am(iSym)*(J-1) + 1
               Call ddaFile(lUnit_F(iSym,iTyp),iOpt,V,lTot,iAdr)
               Do iSymi = 1,nSym
                  iSyma = MulD2h(iSymi,iSym)
                  na = max(nVir(iSyma),1)
                  Call dGeMV_('T',nVir(iSyma),nOcc(iSymi),
     &                       1.0d0,V(1+iT1Am(iSyma,iSymi)),na,
     &                             PVir(1+iVir(iSyma)),1,
     &                       0.0d0,Q,1)
                  Y(J) = Y(J) + dDot_(nOcc(iSymi),Q,1,
     &                                POcc(1+iOcc(iSymi)),1)
               End Do
            End Do
            iOpt = 2
            Call ChoMP2_OpenF(iOpt,iTyp,iSym)
            Call mma_deallocate(V)

C           Calculate errors.
C           -----------------

            Call mma_allocate(D,nMP2Vec(iSym),Label='D')
            Call dCopy_(nMP2Vec(iSym),X,1,D,1)
            Call dAXPY_(nMP2Vec(iSym),-1.0d0,Y,1,D,1)
            Err(1,iSym) = abs(D(1))
            Err(2,iSym) = abs(D(1))
            Err(3,iSym) =     D(1)
            Do J = 1,nMP2Vec(iSym)-1
               Err(1,iSym) = min(Err(1,iSym),abs(D(1+J)))
               Err(2,iSym) = max(Err(2,iSym),abs(D(1+J)))
               Err(3,iSym) =     Err(3,iSym)   + D(1+J)
            End Do
            AvgErr = AvgErr + Err(3,iSym)
            Err(3,iSym) = Err(3,iSym)/dble(nMP2Vec(iSym))
            Err(4,iSym) = dDot_(nMP2Vec(iSym),D,1,D,1)
            RMSErr = RMSErr + Err(4,iSym)
            Err(4,iSym) = Err(4,iSym)/dble(nMP2Vec(iSym))
            Err(4,iSym) = sqrt(Err(4,iSym))
            Call mma_deallocate(D)

C           Deallocation.
C           -------------

            Call mma_deallocate(Y)
            Call mma_deallocate(X)

         Else

            Call FZero(Err(1,iSym),4)

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

      Call mma_deallocate(Q)
      Call mma_deallocate(PVir)
      Call mma_deallocate(POcc)

      End
