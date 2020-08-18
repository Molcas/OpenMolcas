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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE Restore_mat(n,m,lu_A0,lu_A,iD_A,Scr,lScr,Add0s)
************************************************************************
*
*     Author:  F. Aquilante
*
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer n, m, lu_A0, lu_A, iD_A(n), lScr
      Real*8  Scr(lScr)
      Logical Add0s
#include "warnings.fh"

      lmax=lScr-n
      If (lmax .lt. n) Then
         Call WarningMessage(2,'Error in Restore_mat')
         write(6,*) ' Restore_mat: too little scratch space!! '
         Call Quit(_RC_CHO_LOG_)
      Endif
*
      nMem_Col = m
      mNeed = nMem_Col*(nMem_Col+1)/2
      Do while (mNeed .gt. lmax)
         mNeed = mNeed - nMem_Col
         nMem_Col = nMem_Col - 1
      End Do
*
      kAddr=0
      ij=nMem_Col*(nMem_Col+1)/2
      Call dDaFile(lu_A0,2,Scr(1),ij,kAddr)
*
      iOff=0
      Do kCol=1,nMem_Col
         Do i=1,kCol
            iCol=iD_A(i)
            iScr=ij+iCol
            jCol=iOff+i
            Scr(iScr)=Scr(jCol)
         End Do
         Do i=kCol+1,n
            iCol=iD_A(i)
            iScr=ij+iCol
            Scr(iScr)=0.0d0
         End Do
         iAddr=n*(kCol-1)
         Call dDaFile(lu_A,1,Scr(ij+1),n,iAddr)
C        Call RecPrt('QVec',' ',Scr(ij+1),1,n)
         iOff=iOff+kCol
      End Do
*
      Do kCol=nMem_Col+1,m
         Call dDaFile(lu_A0,2,Scr(1),kCol,kAddr)
         Do i=1,kCol
            iCol=iD_A(i)
            iScr=n+iCol
            Scr(iScr)=Scr(i)
         End Do
         Do i=kCol+1,n
            iCol=iD_A(i)
            iScr=n+iCol
            Scr(iScr)=0.0d0
         End Do
         iAddr=n*(kCol-1)
         Call dDaFile(lu_A,1,Scr(n+1),n,iAddr)
C        Call RecPrt('QVec',' ',Scr(n+1),1,n)
      End Do
*
      If (Add0s) Then
         Do kCol=m+1,n   ! linearly dependent cols
            iAddr=n*(kCol-1)
            Call FZero(Scr,n)
            Call dDaFile(lu_A,1,Scr(1),n,iAddr)
         End Do
      EndIf
*
      Return
      End
