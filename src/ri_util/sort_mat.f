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
      SUBROUTINE SORT_mat(irc,Diag,nDim,nVec,iD_A,nSym,lu_A0,mode,
     &                        lScr,Scr)
************************************************************************
*
*     Author:  F. Aquilante
*
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer irc, nSym, lScr
      Integer iD_A(*), nDim(nSym), nVec(nSym), lu_A0(nSym)
      Real*8  Scr(lScr), Diag(*)
      Character*7 mode
#include "WrkSpc.fh"
      Character Name_A*6
*
C     Write (6,*) 'Mode=',Mode
      irc=0
      If (mode.eq.'GePivot') Then  ! returns iD_A
         is=1
*19112013VVP: The threshold changed from 1.d-14 to 1.d-12
        Thr=1.0D-12
*The original threshold:
*        Thr=1.0D-14
         Do iSym=1,nSym
            If (nDim(iSym).eq.0) Go To 81
            lu_A=7
            Write(Name_A,'(A4,I2.2)') 'ZMAT',iSym-1
            Call DaName_MF_WA(lu_A,Name_A)
C           Call RecPrt('Diag',' ',Diag(iS),1,nDim(iSym))
            Call get_pivot_idx(Diag(is),nDim(iSym),nVec(iSym),
     &                         lu_A0(iSym),lu_A,
     &                         iD_A(is),Scr,lScr,Thr)
            Call DaEras(lu_A) ! we do not need it
C           Call RecPrt('Diag',' ',Diag(iS),1,nDim(iSym))
            is=is+nDim(iSym)
 81         Continue
         End Do
      ElseIf (mode.eq.'DoPivot') Then ! store full-pivoted UT A-matrix
         is=1
         Do iSym=1,nSym
            If (nVec(iSym).eq.0) Go To 82
            lu_A=7
            Write(Name_A,'(A4,I2.2)') 'AMAT',iSym-1
            Call DaName_MF_WA(lu_A,Name_A)
            Call Pivot_mat(nDim(iSym),nVec(iSym),lu_A0(iSym),lu_A,
     &                     iD_A(is),Scr,lScr)
            Call DaEras(lu_A0(iSym))
            lu_A0(iSym)=lu_A
 82         Continue
            is=is+nDim(iSym)
         End Do

      ElseIf (mode.eq.'Restore') Then !store squared Q-mat (col. piv.)
         is=1
         Do iSym=1,nSym
            If (nVec(iSym).eq.0) Go To 83
            lu_A=7
            Write(Name_A,'(A4,I2.2)') 'QVEC',iSym-1
            Call DaName_MF_WA(lu_A,Name_A)
            Call Restore_mat(nDim(iSym),nVec(iSym),lu_A0(iSym),lu_A,
     &                       iD_A(is),Scr,lScr,.false.)
            Call DaEras(lu_A0(iSym))
            lu_A0(iSym)=lu_A
 83         Continue
            is=is+nDim(iSym)
         End Do

      Else
        write(6,*)' SORT_mat: invalid mode! '
        irc=66
      EndIf

      Return
      End
