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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************
      SubRoutine Cho_SOSmp2_DecChk(irc,iSym,Col,nDim,nCol,Wrk,lWrk,     &
     &                             ErrStat)
!
!     Francesco Aquilante, May 2007.
!
!     Purpose: check Scaled Opposite-Spin MP2 decomposition of
!              M(ai,bj)=(ai|bj)^2 for sym. block iSym.
!              The columns of the M matrix (actually, the sqrt of
!              their elements) are compared nCol columns at a time
!              with the corresponding integrals obtained from the
!              Cholesky (or RI) vectors.
!              The memory requirement of this routine
!              should be limited to approximately that of the
!              decomposition itself. Note, however, that since all
!              elements are computed, this routine will consume
!              significantly more CPU time.
!              Files are assumed open.
!              On exit,
!              ErrStat(1) = min error
!              ErrStat(2) = max error
!              ErrStat(3) = rms error
!
      use ChoMP2, only: OldVec
#include "implicit.fh"
      Real*8  Col(nDim,nCol), Wrk(lWrk), ErrStat(3)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_dec.fh"
#include "WrkSpc.fh"

      external ddot_

      Character*6  ThisNm
      Character*17 SecNam
      Parameter (SecNam = 'Cho_SOSmp2_DecChk', ThisNm = 'DecChk')

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

!        Compute  M(ai,bj)=(ai|bj)^2  from "new" vectors.
!        -----------------------------------------------

         lU     = lUnit_F(iSym,2)
         NumVec = nMP2Vec(iSym)
         Fac    = 0.0d0
         Call ChoMP2_DecChk_Int(irc,lU,Col,Nai,Nbj,ibj1,NumVec,         &
     &                          Wrk,lWrk,Fac)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': Cho_SOSmp2_DecChk_Int  rc= ',irc,' [1]'
            irc = 1
            Return
         End If

!        Obtain the (ai|bj) integrals from M(ai,bj) .
!        --------------------------------------------

         Do kbj = 1,Nbj
            Do kai = 1,Nai
               Col(kai,kbj) = sqrt(Col(kai,kbj))
            End Do
         End Do


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
               Write(6,*)SecNam,': Cho_SOSmp2_DecChk_Int returned ',irc,&
     &                   ' [2]'
               irc = 2
               Return
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

    1   Continue
      End
