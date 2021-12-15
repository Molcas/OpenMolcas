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
* Copyright (C) 2004,2007, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine ChoMP2_Col(Col,nDim,iCol,nCol,Buf,l_Buf)
C
C     Thomas Bondo Pedersen, Dec. 2007.
C
C     Purpose: compute specified (ai|bj) or MP2 amplitude columns.
C
      Implicit None
      Integer nDim, nCol, l_Buf
      Real*8  Col(nDim,nCol), Buf(l_Buf)
      Integer iCol(nCol)
#include "chomp2_dec.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*3  ThisNm
      Character*10 SecNam
      Parameter (SecNam = 'ChoMP2_Col', ThisNm = 'Col')

      Integer iSym

C     Check input.
C     ------------

      If (nCol.lt.1 .or. nDim.lt.1) Return
      iSym = NowSym
      If (nDim .ne. nT1am(iSym)) Then
         Write(6,*) SecNam,': inconsistent dimension. Expected: ',
     &              nT1am(iSym),'   Received: ',nDim
         Write(6,*) SecNam,': symmetry from chomp2_dec.fh: ',iSym
         Call ChoMP2_Quit(SecNam,'inconsistent dimension',' ')
      End If


C     Calculate (ai|bj) integrals.
C     ----------------------------

      Call ChoMP2_IntCol(Col,nDim,iCol,nCol,Buf,l_Buf)

C     Postprocess integrals.
C     ----------------------

      If (iOption_MP2CD .eq. 2) Then ! generate amplitudes
         Call ChoMP2_AmpFromInt(Col,nDim,iCol,nCol,
     &                          Work(ip_EOc),Work(ip_EVir))
      End If

      End
      SubRoutine ChoMP2_AmpFromInt(Col,nDim,iCol,nCol,EOcc,EVir)
C
C     Thomas Bondo Pedersen, Dec. 2007.
C
C     Purpose: scale integrals with orbital energies to get
C              (minus) MP2 amplitudes: (ai|bj)/[e(a)-e(i)+e(b)-e(j)].
C
      Implicit None
      Integer nDim, nCol
      Real*8  Col(nDim,nCol), EOcc(*), EVir(*)
      Integer iCol(nCol)
#include "chomp2_dec.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"

      Integer iSym, bj_, bj, b, iSymb, j, iSymj, iSymi, iSyma, i, ai0
      Integer a, ai
      Real*8  Ebj, DE

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

      iSym = NowSym
      Do bj_ = 1,nCol
         bj = iCol(bj_)
         Call ChoMP2_Col_Invai(bj,iSym,b,iSymb,j,iSymj)
         Ebj = EVir(iVir(iSymb)+b) - EOcc(iOcc(iSymj)+j)
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            Do i = 1,nOcc(iSymi)
               ai0 = iT1am(iSyma,iSymi) + nVir(iSyma)*(i-1)
               Do a = 1,nVir(iSyma)
                  ai = ai0 + a
                  DE = EVir(iVir(iSyma)+a) - EOcc(iOcc(iSymi)+i) + Ebj
                  Col(ai,bj_) = Col(ai,bj_)/DE
               End Do
            End Do
         End Do
      End Do

      End
      SubRoutine ChoMP2_IntCol(Col,nDim,iCol,nCol,Buf,l_Buf)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C     Renamed (from ChoMP2_Col), Thomas Bondo Pedersen, Dec. 2007.
C
C     Purpose: compute specified (ai|bj) columns.
C
#include "implicit.fh"
      Real*8  Col(nDim,nCol), Buf(l_Buf)
      Integer iCol(nCol)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_dec.fh"
#include "WrkSpc.fh"

      Character*6  ThisNm
      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_IntCol', ThisNm = 'IntCol')

      Logical DoClose

      iSym = NowSym
      If (NumCho(iSym) .lt. 1) Then
         Call Cho_dZero(Col,nDim*nCol)
         Return
      End If

      irc = 0

      If (InCore(iSym)) Then  ! old vectors available in core

         Fac = 0.0D0
         Call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,
     &                        Work(ip_OldVec),NumCho(iSym),
     &                        Buf,l_Buf,Fac,irc)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
            Call ChoMP2_Quit(SecNam,'ChoMP2_Col_Comp error','[1]')
         End If

      Else ! old vectors must be read on disk

         DoClose = .false.
         If (lUnit_F(iSym,1) .lt. 1) Then
            Call ChoMP2_OpenF(1,1,iSym)
            DoClose = .true.
         End If

         Call GetMem('MaxCol','Max ','Real',ipWrk,lWrk)

         If (l_Buf .gt. lWrk) Then ! use Buf as work space

            nVec = min(l_Buf/(nDim+1),NumCho(iSym))
            If (nVec .lt. 1) Then
               Write(6,*) SecNam,': insufficient memory for batch!'
               Call ChoMP2_Quit(SecNam,'insufficient memory','[1]')
               nBat = 0
            Else
               nBat = (NumCho(iSym) - 1)/nVec + 1
            End If

            Do iBat = 1,nBat

               If (iBat .eq. nBat) Then
                  NumV = NumCho(iSym) - nVec*(nBat - 1)
               Else
                  NumV = nVec
               End If
               iVec1 = nVec*(iBat - 1) + 1

               iOpt = 2
               lTot = nDim*NumV
               iAdr = nDim*(iVec1 - 1) + 1
               Call ddaFile(lUnit_F(iSym,1),iOpt,Buf(1),lTot,iAdr)

               If (iBat .eq. 1) Then
                  Fac = 0.0D0
               Else
                  Fac = 1.0D0
               End If

               lScr = l_Buf - lTot
               If (lWrk .gt. lScr) Then
                  lWsav = lWrk
                  Call GetMem('ColScr','Allo','Real',ipWrk,lWrk)
                  Call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,
     &                                 Buf(1),NumV,
     &                                 Work(ipWrk),lWrk,Fac,irc)
                  Call GetMem('ColScr','Free','Real',ipWrk,lWrk)
                  lWrk = lWsav
               Else
                  Call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,
     &                                 Buf(1),NumV,
     &                                 Buf(1+lTot),lScr,Fac,irc)
               End If
               If (irc .ne. 0) Then
                  Write(6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
                  Call ChoMP2_Quit(SecNam,'ChoMP2_Col_Comp error','[2]')
               End If

            End Do

         Else ! use Work as work space

            Call GetMem('ColWrk','Allo','Real',ipWrk,lWrk)

            nVec = min(lWrk/nDim,NumCho(iSym))
            If (nVec .lt. 1) Then
               Write(6,*) SecNam,': insufficient memory for batch!'
               Call ChoMP2_Quit(SecNam,'insufficient memory','[2]')
               nBat = 0
            Else
               nBat = (NumCho(iSym) - 1)/nVec + 1
            End If

            Do iBat = 1,nBat

               If (iBat .eq. nBat) Then
                  NumV = NumCho(iSym) - nVec*(nBat - 1)
               Else
                  NumV = nVec
               End If
               iVec1 = nVec*(iBat - 1) + 1

               iOpt = 2
               lTot = nDim*NumV
               iAdr = nDim*(iVec1 - 1) + 1
               Call ddaFile(lUnit_F(iSym,1),iOpt,Work(ipWrk),lTot,iAdr)

               If (iBat .eq. 1) Then
                  Fac = 0.0D0
               Else
                  Fac = 1.0D0
               End If

               lScr = lWrk - lTot
               If (l_Buf .gt. lScr) Then
                  Call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,
     &                                 Work(ipWrk),NumV,
     &                                 Buf(1),l_Buf,Fac,irc)
               Else
                  Call ChoMP2_Col_Comp(Col,nDim,iCol,nCol,
     &                                 Work(ipWrk),NumV,
     &                                 Work(ipWrk+lTot),lScr,Fac,irc)
               End If
               If (irc .ne. 0) Then
                  Write(6,*) SecNam,': ChoMP2_Col_Comp returned ',irc
                  Call ChoMP2_Quit(SecNam,'ChoMP2_Col_Comp error','[3]')
               End If

            End Do

            Call GetMem('ColWrk','Free','Real',ipWrk,lWrk)

         End If

         If (DoClose) Then
            Call ChoMP2_OpenF(2,1,iSym)
            DoClose = .false.
         End If

      End If

      End
