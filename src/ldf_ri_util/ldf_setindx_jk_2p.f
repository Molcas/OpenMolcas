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
      Subroutine LDF_SetIndx_JK_2P(AB,CD)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Set index arrays for evaluating (J_AB | K_CD) integrals.
C     To unset: Call LDF_UnsetIndx_JK_2P().
C
C     Makes use of LDF_SetIndxG for each of the atom pairs; hence,
C     values set by LDF_SetIndxG should be properly unset before calling
C     this routine!
C
C
      Implicit None
      Integer AB, CD
#include "localdf_int.fh"

      Call LDF_SetIndxG(AB)
      Call LDF_Transfer1('AB',ip_IndxG,l_IndxG_1,l_IndxG_2,
     &                        ip_IndxG2,l_IndxG2_1,l_IndxG2_2,
     &                        ip_2CList,l_2CList_1,l_2CList_2)
      ip_IndxG=0
      l_IndxG_1=0
      l_IndxG_2=0
      ip_IndxG2=0
      l_IndxG2_1=0
      l_IndxG2_2=0
      ip_2CList=0
      l_2CList_1=0
      l_2CList_2=0
      Call LDF_SetIndxG(CD)
      Call LDF_Transfer1('CD',ip_IndxG,l_IndxG_1,l_IndxG_2,
     &                        ip_IndxG2,l_IndxG2_1,l_IndxG2_2,
     &                        ip_2CList,l_2CList_1,l_2CList_2)
      ip_IndxG=0
      l_IndxG_1=0
      l_IndxG_2=0
      ip_IndxG2=0
      l_IndxG2_1=0
      l_IndxG2_2=0
      ip_2CList=0
      l_2CList_1=0
      l_2CList_2=0
      Call LDF_SetDim_JK_2P(AB,CD)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SetDim_JK_2P(AB,CD)
      Implicit None
      Integer AB, CD
#include "localdf_int2.fh"

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      nAB=LDF_nBasAux_Pair(AB)
      nCD=LDF_nBasAux_Pair(CD)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Transfer1(Label,
     &                         ip_IndxG,l_IndxG_1,l_IndxG_2,
     &                         ip_IndxG2,l_IndxG2_1,l_IndxG2_2,
     &                         ip_2CList,l_2CList_1,l_2CList_2)
      Implicit None
      Character*2 Label
      Integer ip_IndxG,  l_IndxG_1,  l_IndxG_2
      Integer ip_IndxG2, l_IndxG2_1, l_IndxG2_2
      Integer ip_2CList, l_2CList_1, l_2CList_2
#include "localdf_int2.fh"

      If (Label.eq.'AB') Then
         ip_AB_IndxG=ip_IndxG
         l_AB_IndxG_1=l_IndxG_1
         l_AB_IndxG_2=l_IndxG_2
         ip_AB_IndxG2=ip_IndxG2
         l_AB_IndxG2_1=l_IndxG2_1
         l_AB_IndxG2_2=l_IndxG2_2
         ip_AB_2CList=ip_2CList
         l_AB_2CList_1=l_2CList_1
         l_AB_2CList_2=l_2CList_2
      Else If (Label.eq.'CD') Then
         ip_CD_IndxG=ip_IndxG
         l_CD_IndxG_1=l_IndxG_1
         l_CD_IndxG_2=l_IndxG_2
         ip_CD_IndxG2=ip_IndxG2
         l_CD_IndxG2_1=l_IndxG2_1
         l_CD_IndxG2_2=l_IndxG2_2
         ip_CD_2CList=ip_2CList
         l_CD_2CList_1=l_2CList_1
         l_CD_2CList_2=l_2CList_2
      Else
         Call WarningMessage(2,'LDF_Transfer1: unknown Label')
         Call LDF_Quit(1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Transfer2(ip_IndxG_,l_IndxG_1_,l_IndxG_2_,
     &                         ip_IndxG2_,l_IndxG2_1_,l_IndxG2_2_,
     &                         ip_2CList_,l_2CList_1_,l_2CList_2_)
      Implicit None
      Integer ip_IndxG_,  l_IndxG_1_,  l_IndxG_2_
      Integer ip_IndxG2_, l_IndxG2_1_, l_IndxG2_2_
      Integer ip_2CList_, l_2CList_1_, l_2CList_2_
#include "localdf_int.fh"

      ip_IndxG=ip_IndxG_
      l_IndxG_1=l_IndxG_1_
      l_IndxG_2=l_IndxG_2_
      ip_IndxG2=ip_IndxG2_
      l_IndxG2_1=l_IndxG2_1_
      l_IndxG2_2=l_IndxG2_2_
      ip_2CList=ip_2CList_
      l_2CList_1=l_2CList_1_
      l_2CList_2=l_2CList_2_

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetIndx_JK_2P()
      Implicit None
#include "localdf_int2.fh"

      Call LDF_Transfer2(ip_AB_IndxG,l_AB_IndxG_1,l_AB_IndxG_2,
     &                   ip_AB_IndxG2,l_AB_IndxG2_1,l_AB_IndxG2_2,
     &                   ip_AB_2CList,l_AB_2CList_1,l_AB_2CList_2)
      Call LDF_UnsetIndxG()
      Call LDF_Transfer2(ip_CD_IndxG,l_CD_IndxG_1,l_CD_IndxG_2,
     &                   ip_CD_IndxG2,l_CD_IndxG2_1,l_CD_IndxG2_2,
     &                   ip_CD_2CList,l_CD_2CList_1,l_CD_2CList_2)
      Call LDF_UnsetIndxG()

      ip_AB_IndxG=0
      l_AB_IndxG_1=0
      l_AB_IndxG_2=0
      ip_AB_IndxG2=0
      l_AB_IndxG2_1=0
      l_AB_IndxG2_2=0
      ip_AB_2CList=0
      l_AB_2CList_1=0
      l_AB_2CList_2=0

      ip_CD_IndxG=0
      l_CD_IndxG_1=0
      l_CD_IndxG_2=0
      ip_CD_IndxG2=0
      l_CD_IndxG2_1=0
      l_CD_IndxG2_2=0
      ip_CD_2CList=0
      l_CD_2CList_1=0
      l_CD_2CList_2=0

      nAB=0
      nCD=0

      End
