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
************************************************************************
      SubRoutine Sort(EVal,EVec,n,nB)
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
#define _DEBUG_
#ifdef _DEBUG_
      Call qEnter('Sort')
#endif
*
      write(6,*) 'Sort 0'
      Do 100 i = 1,n - 1
      write(6,*) 'Sort 1'
         k = i
         Do 110 j = i + 1, n
            If (EVal(j).lt.EVal(k)) k = j
            write(6,*) 'Sort 2'
110      Continue
         If (k.ne.i) Then
            write(6,*) 'Sort 3'
            Swap    = EVal(k)
            EVal(k) = EVal(i)
            EVal(i) = Swap
            Do 120 l = 1, nB
            write(6,*) 'Sort 4'
               Swap      =   EVec(l,k)
               EVec(L,K) = - EVec(l,i)
               EVec(L,I) =   Swap
120         Continue
         End If
100   Continue
*
#ifdef _DEBUG_
      Call qExit('Sort')
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
