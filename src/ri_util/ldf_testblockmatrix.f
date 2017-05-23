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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Logical Function LDF_TestBlockMatrix(ip_Blocks,Packed,FullMatrix)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: test block matrix routines.
C
      Implicit None
      Integer ip_Blocks
      Logical Packed
      Real*8  FullMatrix(*)
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      real*8 ddot_
      external ddot_

      Real*8 Tol
      Parameter (Tol=1.0d-12)

      Real*8  z
      Integer ip, l

      ! Allocate a full matrix
      If (Packed) Then
         l=nBas_Valence*(nBas_Valence+1)/2
      Else
         l=nBas_Valence*nBas_Valence
      End If
      Call GetMem('TBMTst','Allo','Real',ip,l)

      ! Get full matrix from blocked
      Call LDF_Blocked2Full(ip_Blocks,Packed,Work(ip))

      ! Check that it is the same as FullMatrix
      Call dAXPY_(l,-1.0d0,FullMatrix,1,Work(ip),1)
      z=sqrt(dDot_(l,Work(ip),1,Work(ip),1))
      If (z.gt.Tol) Then
         LDF_TestBlockMatrix=.False.
      Else
         LDF_TestBlockMatrix=.True.
      End If

      ! Deallocate
      Call GetMem('TBMTst','Free','Real',ip,l)

      End
