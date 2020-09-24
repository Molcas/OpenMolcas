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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               Anders Ohrn                                            *
************************************************************************
      SUBROUTINE JACORD(HH,EIGVEC,NVEC,NDIM)
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      DIMENSION HH(*), EIGVEC(NDIM,NVEC)

      ThrZ=1.0D-14
      DO 100 I=1,NVEC-1
        II=(I*(I+1))/2
        EI=HH(II)
        EMIN=EI
        IMIN=I
        DO 10 J=I+1,NVEC
          JJ=(J*(J+1))/2
          EJ=HH(JJ)
          IF(EJ.GE.EMIN.or.ABS(EJ-EMIN).lt.ThrZ) GOTO 10
          EMIN=EJ
          IMIN=J
  10    CONTINUE
        IF(IMIN.EQ.I) GOTO 100
        HH(II)=EMIN
        JJ=(IMIN*(IMIN+1))/2
        HH(JJ)=EI
        DO 20 K=1,NDIM
          SWAP=EIGVEC(K,I)
          EIGVEC(K,I)=EIGVEC(K,IMIN)
          EIGVEC(K,IMIN)=SWAP
  20    CONTINUE
 100  CONTINUE
      RETURN
      END
      SUBROUTINE JACORD2(EVal,EVec,n,nB)
C      SubRoutine Sort(EVal,EVec,n,nB)
************************************************************************
*                                                                      *
*     purpose: Sort the set of eigenvalues and eigenvectors            *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       EVal    : the set of eigenvalues in random order               *
*       EVec    : the set of eigenvectors in random order              *
*       n,nB    : dimensions                                           *
*                                                                      *
*     output:                                                          *
*       EVal    : sorted set of eigenvalues                            *
*       EVec    : sorted set of eigenvectors                           *
*                                                                      *
*     called from: DCore, NewOrb, IvoGen                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
      Dimension EVal(n),EVec(nB,n)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
#endif
*
      Do 100 i = 1,n - 1
         k = i
         Do 110 j = i + 1, n
            If (EVal(j).lt.EVal(k)) k = j
110      Continue
         If (k.ne.i) Then
            Swap    = EVal(k)
            EVal(k) = EVal(i)
            EVal(i) = Swap
            Do 120 l = 1, nB
               Swap      =   EVec(l,k)
               EVec(L,K) =   EVec(l,i)
               EVec(L,I) =   Swap
120         Continue
         End If
100   Continue
*
#ifdef _DEBUG_
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
      SUBROUTINE JACORD3(EVal,EVec,n,nB)
C      SubRoutine Sort(EVal,EVec,n,nB)
************************************************************************
*                                                                      *
*     purpose: Sort the set of eigenvalues and eigenvectors, in        *
*              falling sequence.                                       *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       EVal    : the set of eigenvalues in random order               *
*       EVec    : the set of eigenvectors in random order              *
*       n,nB    : dimensions                                           *
*                                                                      *
*     output:                                                          *
*       EVal    : sorted set of eigenvalues                            *
*       EVec    : sorted set of eigenvectors                           *
*                                                                      *
*     called from: DCore, NewOrb, IvoGen                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: A. Ohrn copied JACORD2 and made minimal modification.   *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
      Dimension EVal(n),EVec(nB,n)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
#endif
*
      Do 100 i = 1,n - 1
         k = i
         Do 110 j = i + 1, n
            If (EVal(j).gt.EVal(k)) k = j
110      Continue
         If (k.ne.i) Then
            Swap    = EVal(k)
            EVal(k) = EVal(i)
            EVal(i) = Swap
            Do 120 l = 1, nB
               Swap      =   EVec(l,k)
               EVec(L,K) =   EVec(l,i)
               EVec(L,I) =   Swap
120         Continue
         End If
100   Continue
*
#ifdef _DEBUG_
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
