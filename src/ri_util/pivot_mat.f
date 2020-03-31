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
      SUBROUTINE Pivot_mat(n,m,lu_A0,lu_A,iD_A,Scr,lScr)
************************************************************************
*
*     Author:  F. Aquilante
*
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer n, m, lu_A0, lu_A, iD_A(n), lScr
      Real*8  Scr(lScr)
#include "warnings.fh"

      lmax=lScr-n
      If (lmax .lt. n) Then
         Call WarningMessage(2,'Error in Pivot_mat')
         write(6,*) ' Pivot_mat: too little scratch space !!'
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
      iScr=n
      Do kCol=1,nMem_Col
         jCol=iD_A(kCol)
         iAddr=n*(jCol-1)
         Call dDaFile(lu_A0,2,Scr(1),n,iAddr)
         Do i=1,kCol
            iCol=iD_A(i)
            iScr=iScr+1
            Scr(iScr)=Scr(iCol)
         End Do
      End Do
      kAddr=0
      ij=nMem_Col*(nMem_Col+1)/2
      Call dDaFile(lu_A,1,Scr(n+1),ij,kAddr)
*
      Do kCol=nMem_Col+1,m
         jCol=iD_A(kCol)
         iAddr=n*(jCol-1)
         Call dDaFile(lu_A0,2,Scr(1),n,iAddr)
         Do i=1,kCol
            iCol=iD_A(i)
            iScr=n+i
            Scr(iScr)=Scr(iCol)
         End Do
         Call dDaFile(lu_A,1,Scr(n+1),kCol,kAddr)
      End Do
*
      Return
      End
