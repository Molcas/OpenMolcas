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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_PrtInt(iSCD,iSAB,xInt,lInt)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: Print integral shell quadruple (IfcSew=2 or 3).
C
      use ChoArr, only: iSP2F, nBstSh
      Implicit None
      Integer iSCD, iSAB
      Integer lInt
      Real*8  xInt(lInt)
#include "cholesky.fh"
#include "choptr.fh"
#include "chosew.fh"
#include "WrkSpc.fh"

      Character*10 SecNam
      Parameter (SecNam='Cho_PrtInt')

      Integer nRow(8)
      Integer iSC, iSD, iSA, iSB
      Integer nCD, nAB
      Integer AB, CD
      Integer iAB, iCD
      Integer iSym
      Integer kOffI

      Real*8 xNorm

      Integer i, j
      Integer iShP2RS, iShP2Q
      iShP2RS(i,j)=iWork(ip_iShP2RS-1+2*(j-1)+i)
      iShP2Q(i,j)=iWork(ip_iShP2Q-1+2*(j-1)+i)

      ! Set row dimension
      If (IfcSew.eq.2) Then
         Do iSym=1,nSym
            nRow(iSym)=nnBstR(iSym,2)
         End Do
      Else If (IfcSew.eq.3) Then
         Do iSym=1,nSym
            nRow(iSym)=nDim_Batch(iSym)
         End Do
      Else
         Call Cho_Quit(SecNam//': Illegal IfcSew',103)
         Do iSym=1,nSym ! avoid compiler warnings
            nRow(iSym)=0
         End Do
      End If

      ! Get full shell pair dimensions
      Call Cho_InvPck(iSP2F(iSCD),iSC,iSD,.True.)
      If (iSC.eq.iSD) Then
         nCD=nBstSh(iSC)*(nBstSh(iSC)+1)/2
      Else
         nCD=nBstSh(iSC)*nBstSh(iSD)
      End If
      Call Cho_InvPck(iSP2F(iSAB),iSA,iSB,.True.)
      If (iSA.eq.iSB) Then
         nAB=nBstSh(iSA)*(nBstSh(iSA)+1)/2
      Else
         nAB=nBstSh(iSA)*nBstSh(iSB)
      End If

      ! Loop through integral shell quadruple
      Write(LuPri,'(//,A,I4,A,I4,A,I4,A,I4,A)')
     & 'Shell Quadruple (',iSC,',',iSD,'|',iSA,',',iSB,'):'
      Do AB=1,nAB
         iAB=iShP2Q(1,AB)
         If (iAB.gt.0) Then
            iSym=iShP2Q(2,AB)
            kOffI=iOff_Col(iSym)+nRow(iSym)*(iAB-1)
            xNorm=0.0d0
            Do CD=1,nCD
               iCD=iShP2RS(1,CD)
               If (iCD.gt.0) Then
                  If (iShP2RS(2,CD).eq.iSym) Then
                     Write(Lupri,'(2X,A,I4,A,I4,A,1P,D15.6)')
     &               '(',CD,'|',AB,') =',xInt(kOffI+iCD)
                     xNorm=xNorm+xInt(kOffI+iCD)**2
                  End If
               End If
            End Do
            Write(Lupri,'(A,I4,A,1P,D15.6)')
     &      '**Norm of column',AB,':',sqrt(xNorm)
         End If
      End Do

      End
