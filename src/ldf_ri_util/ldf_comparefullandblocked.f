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
      Subroutine LDF_CompareFullAndBlocked(Tol,PackedF,F,ip_FBlocks,irc)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compare two block matrices.
C
      Implicit None
      Real*8  Tol
      Logical PackedF
      Real*8  F(*)
      Integer ip_FBlocks
      Integer irc
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      Logical  BlockMatricesAreIdentical
      External BlockMatricesAreIdentical

      Integer ip_XBlocks
      Integer ip, l
      Integer iCount, i

      ! Init return code
      irc=0

C----------------------------------------------------------------
C     Check that the block matrix obtained from F is identical to
C     ip_FBlocks.
C----------------------------------------------------------------

      Call LDF_AllocateBlockMatrix('_F_',ip_XBlocks)
      Call LDF_Full2Blocked(F,PackedF,ip_XBlocks)
      If (.not.BlockMatricesAreIdentical(ip_XBlocks,ip_FBlocks,Tol))
     & Then
         irc=irc+1
      End If
      Call LDF_DeallocateBlockMatrix('_F_',ip_XBlocks)

C---------------------------------------------------------------------
C     Check that the full matrix obtained from ip_FBlocks is identical
C     to F.
C---------------------------------------------------------------------

      If (PackedF) Then
         l=nBas_Valence*(nBas_Valence+1)/2
      Else
         l=nBas_Valence**2
      End IF
      Call GetMem('F__','Allo','Real',ip,l)
      Call Cho_dZero(Work(ip),l)
      Call LDF_Blocked2Full(ip_FBlocks,PackedF,Work(ip))
      iCount=0
      Do i=1,l
         If (abs(F(i)-Work(ip-1+i)).gt.Tol) Then
            iCount=iCount+1
         End If
      End Do
      If (iCount.gt.0) Then
         irc=irc+2
      End If
      Call GetMem('F__','Free','Real',ip,l)

      End
