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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_GetBlockedOverlapMatrix(iOpt,ip_Blocks)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: get the overlap matrix stored as a block matrix.
C              Options:
C              iOpt=0: read the full overlap matrix from disk and
C                      reorder to block matrix.
C              iOpt=1: compute each block individually.
C                      (not implemented yet.)
C
      Implicit None
      Integer iOpt
      Integer ip_Blocks

      Character*27 SecNam
      Parameter (SecNam='LDF_GetBlockedOverlapMatrix')

      If (iOpt.eq.0) Then
         Call LDF_GetBlockedOverlapMatrix_0(ip_Blocks)
      Else If (iOpt.eq.1) Then
         Write(6,'(A,A,I10,A)')
     &   SecNam,': iOpt=',iOpt,' not implemented!'
         Call LDF_NotImplemented()
      Else
         Call WarningMessage(2,SecNam//': illegal option')
         Write(6,'(A,I10)') 'iOpt=',iOpt
         Call LDF_Quit(1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_GetBlockedOverlapMatrix_0(ip_Blocks)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: get the overlap matrix stored as a block matrix.
C              Read the full overlap matrix from disk and
C              reorder to block matrix.
C
      Implicit None
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      Character*29 SecNam
      Parameter (SecNam='LDF_GetBlockedOverlapMatrix_0')

      Character*8 Label

      Logical Packed
      Integer ip_S, l_S
      Integer irc
      Integer iOpt, iComp, iSyLab

      Packed=.True.
      l_S=nBas_Valence*(nBas_Valence+1)/2+4
      Call GetMem('LDFOVLP','Allo','Real',ip_S,l_S)
      irc=-1
      iOpt=2
      Label='Mltpl  0'
      iComp=1
      iSyLab=1
      Call RdOne(irc,iOpt,Label,iComp,Work(ip_S),iSyLab)
      If (irc.ne.0) Then
         Call WarningMessage(2,
     &                      SecNam//': non-zero return code from RdOne')
         Write(6,'(A,I10)') 'irc=',irc
         Call LDF_Quit(1)
      End If
      Call LDF_Full2Blocked(Work(ip_S),Packed,ip_Blocks)
      Call GetMem('LDFOVLP','Free','Real',ip_S,l_S)

      End
