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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Cho_P_UpdateBookmarks(iRS)
C
C     Thomas Bondo Pedersen, August 2012.
C
C     Update bookmarks for reduced set (integral pass) iRS:
C        - integral accuracy (max diagonal)
C        - number of Cholesky vectors
C
C     Note: it is assumed that array DiaMax and number of Cholesky
C     vectors are properly updated before calling this routine.
C
      Implicit None
      Integer iRS
#include "cho_para_info.fh"
#include "choglob.fh"
#include "chobkm.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"

      If (l_BkmVec.lt.1 .or. l_BkmThr.lt.1) Return

      If (Cho_Real_Par) Then
         Call Cho_UpdateBookmarks(iRS,nSym,MaxRed,NumCho_G,DiaMaxT,
     &                            iWork(ip_BkmVec),Work(ip_BkmThr))
      Else
         Call Cho_UpdateBookmarks(iRS,nSym,MaxRed,NumCho,DiaMaxT,
     &                            iWork(ip_BkmVec),Work(ip_BkmThr))
      End If
      nCol_BkmVec=nCol_BkmVec+1
      nCol_BkmThr=nCol_BkmThr+1

      End
      Subroutine Cho_UpdateBookmarks(iRS,nSym,nRS,nVec,delta,nBkmVec,
     &                               BkmThr)
      Implicit None
      Integer iRS
      Integer nRS
      Integer nSym
      Integer nVec(nSym)
      Real*8  delta(nSym)
      Integer nBkmVec(nSym,nRS)
      Real*8  BkmThr(nSym,nRS)

      Integer iSym

      Do iSym=1,nSym
         nBkmVec(iSym,iRS)=nVec(iSym)
      End Do

      Do iSym=1,nSym
         BkmThr(iSym,iRS)=delta(iSym)
      End Do

      End
