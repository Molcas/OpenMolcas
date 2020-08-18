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
      SUBROUTINE get_pivot_idx_w(Diag,Wg,n,m,lu_A0,lu_A,iD_A,Scr,lScr,
     &                           Thr)
************************************************************************
*
*     Author:  F. Aquilante
*
*       Note:  this routine differs from Get_pivot_idx because here
*              the pivoting/convergence is decided based on weighted
*              diagonals
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer n, m, lu_A0, lu_A, iD_A(n), lScr
      Real*8  Diag(*), Wg(*), Scr(lScr)
#include "WrkSpc.fh"
#include "warnings.fh"
*
*
      Acc=Min(1.0D-12,thr*1.0D-2)
      Call GetMem('List','Allo','Inte',list,n)
      Do i=0,n-1
         iWork(list+i)=i+1
      End Do
*
      lmax=lScr-2*n
      If (lmax .lt. n) Then
         Call WarningMessage(2,'Error in Get_Pivot_idx_w')
         write(6,*) ' Get_Pivot_idx_w: too little scratch space!! '
         Call Quit(_RC_CHO_LOG_)
      Endif
*
      nMem_Col = Min(lmax/n,n)
*
      kAddr=0
      is=1+n
      ij = n*nMem_Col
      ks=is+ij
      kScr=lScr-n-ij
*
      m=0
      Do kCol = 1,n
*
         iD_Col=0
         XMax=0.0D0
         Do i=1,n
            If (Abs(Diag(i)*Wg(i)).gt.xMax+Acc) Then
               iD_Col=i
               xMax=Abs(Diag(i))
            End If
         End Do
         If (iD_Col.lt.0 .or. iD_Col.gt.n) Then
           Write(6,*) 'Get_Pivot_idx_w: Index of MaxDiag out of bounds!'
           Write(6,*) 'iD_Col = ',iD_Col
           Call Abend()
         ElseIf (iD_Col.eq.0) Then
           Go To 100
         End If
         iD_A(kCol) = iD_Col ! set the mapping
*
         js=n*(kCol-1)+is ! overlay A and Z
         If (kCol.gt.nMem_Col) js=1
*
         kAddr=n*(iD_Col-1)
         Call dDaFile(lu_A0,2,Scr(js),n,kAddr)
*
         Call CHO_FACTOR(Diag,Scr(js),iD_A,kCol,n,Scr(is),nMem_Col,lu_A,
     &                   Scr(ks),kScr,thr,lindep)
*
         If (lindep.ne.0) Goto 100
*
         iWork(list+iD_Col-1)=0
         m=m+1
*
         iAddr=n*(kCol-1)
         If (kCol.gt.nMem_Col) Call dDaFile(lu_A,1,Scr(1),n,iAddr)
*
      End Do
*
100   Continue
      iAddr=0
      Call dDaFile(lu_A,1,Scr(is),ij,iAddr)
*
      If (m.lt.n) Then
         istart=1
         Do k=m+1,n
            Do i=istart,n
               if (iWork(list+i-1).ne.0) Then
                  iD_A(k)=i
                  istart=i+1
                  goto 200
               endif
            End Do
200         Continue
         End Do
      ElseIf (m.gt.n) Then
         Write(6,*) 'Get_Pivot_idx_w: m > n is not possible!'
         Call Abend()
      EndIf
      Call GetMem('List','Free','Inte',list,n)
*
      Return
      End
