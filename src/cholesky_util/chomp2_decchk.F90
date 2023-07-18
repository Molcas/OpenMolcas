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
! Copyright (C) 2005,2008, Thomas Bondo Pedersen                       *
!***********************************************************************
      SubRoutine ChoMP2_DecChk(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
!
!     Thomas Bondo Pedersen, Jan. 2008.
!
!     Purpose: check Cholesky decomposition of the (ai|bj) integrals
!              or MP2 amplitudes (sym. block iSym).
!              The columns of the matrix are
!              compared nCol columns at a time. This
!              implies that the memory requirement of this routine
!              should be limited to approximately that of the
!              decomposition itself. Note, however, that since all
!              integrals are computed, this routine will consume
!              significantly more CPU time.
!              Files are assumed open.
!              On exit,
!              ErrStat(1) = min error
!              ErrStat(2) = max error
!              ErrStat(3) = rms error
!
      use ChoMP2_dec, only: iOption_MP2CD
      Implicit None
      Integer irc, iSym
      Integer nDim, nCol
      Real*8  Col(nDim,nCol)
      Integer lWrk
      Real*8  Wrk(lWrk)
      Real*8  ErrStat(3)

      Character*6  ThisNm
      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_DecChk', ThisNm = 'DecChk')


      If (iOption_MP2CD .eq. 1) Then ! (ai|bj) int
         Call ChoMP2_DecChk_1(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
      Else If (iOption_MP2CD .eq. 2) Then ! MP2 amp
         Call ChoMP2_DecChk_2(irc,iSym,Col,nDim,nCol,Wrk,lWrk,ErrStat)
      Else
         Write(6,*) SecNam,': WARNING! ',                               &
     &              'Unknown option, iOption_MP2CD = ', iOption_MP2CD
         irc = -123456
      End If


      End
      SubRoutine ChoMP2_DecChk_1(irc,iSym,Col,nDim,nCol,Wrk,lWrk,       &
     &                           ErrStat)
!
!     Thomas Bondo Pedersen, Jan. 2005.
!
!     Purpose: check MP2 decomposition of the (ai|bj) integrals
!              (sym. block iSym). The columns of the (ai|bj) matrix are
!              compared nCol columns at a time. This
!              implies that the memory requirement of this routine
!              should be limited to approximately that of the
!              decomposition itself. Note, however, that since all
!              integrals are computed, this routine will consume
!              significantly more CPU time.
!              Files are assumed open.
!              On exit,
!              ErrStat(1) = min error
!              ErrStat(2) = max error
!              ErrStat(3) = rms error
!
      use ChoMP2, only: OldVec
      use ChoMP2_dec, only: Incore
      Implicit Real*8 (a-h,o-z)
      Real*8  Col(nDim,nCol), Wrk(lWrk), ErrStat(3)
#include "cholesky.fh"
#include "chomp2.fh"

      external ddot_

      Character*8  ThisNm
      Character*15 SecNam
      Parameter (SecNam = 'ChoMP2_DecChk_1', ThisNm = 'DecChk_1')

      irc = 0

!     Check dimensions.
!     -----------------

      If (nDim.lt.1 .or. nCol.lt.1) Return
      If (nDim .ne. nT1am(iSym)) Then
         irc = -1
         Go To 1 ! exit
      End If

!     Initialize.
!     -----------

      ErrStat(1) =  9.9d15
      ErrStat(2) = -9.9d15
      ErrStat(3) =  0.0d0

!     Set up batching over columns of the (ai|bj) matrix.
!     ---------------------------------------------------

      NumCol  = min(nCol,nT1am(iSym))
      nBatCol = (nT1am(iSym) - 1)/NumCol + 1

!     Start column batch loop.
!     ------------------------

      Nai = nDim
      Do iBatCol = 1,nBatCol

!        Set batch info.
!        ---------------

         If (iBatCol .eq. nBatCol) Then
            Nbj = nT1am(iSym) - NumCol*(nBatCol - 1)
         Else
            Nbj = NumCol
         End If
         ibj1 = NumCol*(iBatCol - 1) + 1

!        Compute integrals from "new" vectors.
!        -------------------------------------

         lU     = lUnit_F(iSym,2)
         NumVec = nMP2Vec(iSym)
         Fac    = 0.0d0
         Call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,         &
     &                          Wrk,lWrk,Fac)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,' [1]'
            irc = 1
            Go To 1
         End If

!        Compute "old" and subtract "new".
!        ---------------------------------

         If (InCore(iSym)) Then
            Call DGEMM_('N','T',Nai,Nbj,NumCho(iSym),                   &
     &                 1.0d0,OldVec,Nai,OldVec(ibj1),Nai,               &
     &                 -1.0d0,Col,Nai)
         Else
            lU     = lUnit_F(iSym,1)
            NumVec = NumCho(iSym)
            Fac    = -1.0d0
            Call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,      &
     &                             Wrk,lWrk,Fac)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,   &
     &                    ' [2]'
               irc = 2
               Go To 1
            End If
         End If

!        Compute error stats.
!        --------------------

         Do kbj = 1,Nbj
            Do kai = 1,Nai
               ErrStat(1) = min(ErrStat(1),Col(kai,kbj))
               ErrStat(2) = max(ErrStat(2),Col(kai,kbj))
            End Do
         End Do
         ErrStat(3) = ErrStat(3) + dDot_(Nai*Nbj,Col,1,Col,1)

      End Do

!     Compute rms error.
!     ------------------

      xdim = dble(Nai)*dble(Nai)
      ErrStat(3) = sqrt(ErrStat(3)/xdim)

    1 Continue
      End
      SubRoutine ChoMP2_DecChk_2(irc,iSym,Col,nDim,nCol,Wrk,lWrk,       &
     &                           ErrStat)
!
!     Thomas Bondo Pedersen, Jan. 2008.
!
!     Purpose: check MP2 decomposition of the MP2 amplitudes
!              (sym. block iSym). The columns of the matrix are
!              compared nCol columns at a time. This
!              implies that the memory requirement of this routine
!              should be limited to approximately that of the
!              decomposition itself. Note, however, that since all
!              integrals are computed, this routine will consume
!              significantly more CPU time.
!              Files are assumed open.
!              On exit,
!              ErrStat(1) = min error
!              ErrStat(2) = max error
!              ErrStat(3) = rms error
!
      use ChoMP2, only: OldVec
      use ChoMP2_dec, only: EOcc, EVir, Incore
      Implicit Real*8 (a-h,o-z)
      Real*8  Col(nDim,nCol), Wrk(lWrk), ErrStat(3)
#include "cholesky.fh"
#include "chomp2.fh"

      external ddot_

      integer a, b

      Character*8  ThisNm
      Character*15 SecNam
      Parameter (SecNam = 'ChoMP2_DecChk_2', ThisNm = 'DecChk_2')

      MulD2h(k,l)=iEOr(k-1,l-1)+1

      irc = 0

!     Check dimensions.
!     -----------------

      If (nDim.lt.1 .or. nCol.lt.1) Return
      If (nDim .ne. nT1am(iSym)) Then
         irc = -1
         Go To 1 ! exit
      End If

!     Initialize.
!     -----------

      ErrStat(1) =  9.9d15
      ErrStat(2) = -9.9d15
      ErrStat(3) =  0.0d0

!     Set up batching over columns of the (ai|bj) matrix.
!     ---------------------------------------------------

      NumCol  = min(nCol,nT1am(iSym))
      nBatCol = (nT1am(iSym) - 1)/NumCol + 1

!     Start column batch loop.
!     ------------------------

      Nai = nDim
      Do iBatCol = 1,nBatCol

!        Set batch info.
!        ---------------

         If (iBatCol .eq. nBatCol) Then
            Nbj = nT1am(iSym) - NumCol*(nBatCol - 1)
         Else
            Nbj = NumCol
         End If
         ibj1 = NumCol*(iBatCol - 1) + 1

!        Compute amplitudes from "old" vectors.
!        --------------------------------------

         If (InCore(iSym)) Then
            Call DGEMM_('N','T',Nai,Nbj,NumCho(iSym),                   &
     &                 1.0d0,OldVec,Nai,OldVec(ibj1),Nai,               &
     &                 0.0d0,Col,Nai)
         Else
            lU     = lUnit_F(iSym,1)
            NumVec = NumCho(iSym)
            Fac    = 0.0d0
            Call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,      &
     &                             Wrk,lWrk,Fac)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,   &
     &                    ' [2]'
               irc = 2
               Go To 1
            End If
         End If

         ibj0 = ibj1 - 1
         Do iCol = 1,Nbj
            ibj = ibj0 + iCol
            Call ChoMP2_Col_Invai(ibj,iSym,b,iSymb,j,iSymj)
            Ebj = EVir(iVir(iSymb)+b) - Eocc(iOcc(iSymj)+j)
            Do iSymi = 1,nSym
               iSyma = MulD2h(iSymi,iSym)
               Do i = 1,nOcc(iSymi)
                  iai0 = iT1am(iSyma,iSymi) + nVir(iSyma)*(i-1)
                  Do a = 1,nVir(iSyma)
                     iai = iai0 + a
                     DE = EVir(iVir(iSyma)+a) - EOcc(iOcc(iSymi)+i)     &
     &                  + Ebj
                     Col(iai,iCol) = Col(iai,iCol)/DE
                  End Do
               End Do
            End Do
         End Do

!        Compute amplitudes from "new" vectors and subtract "old".
!        ---------------------------------------------------------

         lU     = lUnit_F(iSym,2)
         NumVec = nMP2Vec(iSym)
         Fac    = -1.0d0
         Call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,         &
     &                          Wrk,lWrk,Fac)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_DecChk_Int returned ',irc,' [1]'
            irc = 1
            Go To 1
         End If

!        Compute error stats.
!        --------------------

         Do kbj = 1,Nbj
            Do kai = 1,Nai
               ErrStat(1) = min(ErrStat(1),Col(kai,kbj))
               ErrStat(2) = max(ErrStat(2),Col(kai,kbj))
            End Do
         End Do
         ErrStat(3) = ErrStat(3) + dDot_(Nai*Nbj,Col,1,Col,1)

      End Do

!     Compute rms error.
!     ------------------

      xdim = dble(Nai)*dble(Nai)
      ErrStat(3) = sqrt(ErrStat(3)/xdim)

    1 Continue
      End
