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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_OpenF(iOpt,iTyp,iSym)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: open (iOpt=1), close and keep (iOpt=2), or close and
C              delete (iOpt=3) Cholesky vector files for MP2 program
C              (full vectors).
C              For iOpt=0, the units are initialized (to -1).
C              iTyp=1: transformed Cholesky vectors.
C              iTyp=2: vectors from (ai|bj) decomposition.
C
#include "implicit.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"

      Character*12 SecNam
      Parameter (SecNam = 'ChoMP2_OpenF')

      Character*3 BaseNm
      Character*4 FullNm

      If (iTyp.lt.1 .or. iTyp.gt.nTypF) Then
         Call ChoMP2_Quit(SecNam,'iTyp error',' ')
      End If

C     Initialize units and return for iOpt=0.
C     ---------------------------------------

      If (iOpt .eq. 0) Then
         lUnit_F(iSym,iTyp) = -1
         Return
      End If

C     Open or close files.
C     --------------------

      If (iOpt .eq. 1) Then
         If ((nT1am(iSym).gt.0) .or.
     &       (DoDens .and. (nPQ_prod(iSym).gt.0))) Then
            If (lUnit_F(iSym,iTyp) .lt. 1) Then
               Call ChoMP2_GetBaseNm(BaseNm,iTyp)
               Write(FullNm,'(A3,I1)') BaseNm,iSym
               lUnit_F(iSym,iTyp) = 7
               Call daName_MF_WA(lUnit_F(iSym,iTyp),FullNm)
            End If
         Else
            lUnit_F(iSym,iTyp) = -1
         End If
      Else If (iOpt .eq. 2) Then
         If (lUnit_F(iSym,iTyp) .gt. 0) Then
            Call daClos(lUnit_F(iSym,iTyp))
            lUnit_F(iSym,iTyp) = -1
         End If
      Else If (iOpt .eq. 3) Then
         If (lUnit_F(iSym,iTyp) .gt. 0) Then
            Call daEras(lUnit_F(iSym,iTyp))
            lUnit_F(iSym,iTyp) = -1
         End If
      Else
         Call ChoMP2_Quit(SecNam,'iOpt out of bounds',' ')
      End If

      End
      SubRoutine ChoMP2_GetBaseNm(BaseNm,iTyp)
      Implicit None
      Character*3 BaseNm
      Integer     iTyp

      If (iTyp .eq. 1) Then
         BaseNm = '_AI'
      Else If (iTyp .eq. 2) Then
         BaseNm = '_CD'
      Else
         BaseNm = '_un'
      End If

      End
