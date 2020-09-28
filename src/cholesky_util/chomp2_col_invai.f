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
* Copyright (C) 2007, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Col_Invai(ai,iSymai,a,iSyma,i,iSymi)
C
C     Thomas Bondo Pedersen, Dec. 2007.
C
C     Purpose: calculate indices a and i (incl. symmetries)
C              from compound index ai of symmetry iSymai.
C
      Implicit None
      Integer ai, iSymai, a, iSyma, i, iSymi
#include "cholesky.fh"
#include "chomp2.fh"

      Integer iSym, i_, ai_1, ai_2

#if defined (_DEBUG_)
      Character*16 SecNam
      Parameter (SecNam = 'ChoMP2_Col_Invai')
#endif

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

C     Find iSyma and iSymi.
C     ---------------------

      iSymi = -999999
      iSyma = -999999
      iSym  = nSym + 1
      Do While (iSym .gt. 1)
         iSym  = iSym - 1
         iSymi = iSym
         iSyma = MulD2h(iSymi,iSymai)
         If (nOcc(iSymi).gt.0 .and. nVir(iSyma).gt.0 .and.
     &       ai.gt.iT1Am(iSyma,iSymi)) Then
            iSym = 0 ! Found! -- break loop
         End If
      End Do

#if defined (_DEBUG_)
      If (iSymi.lt.1 .or. iSymi.gt.nSym .or.
     &    iSyma.lt.1 .or. iSyma.gt.nSym) Then
         Call ChoMP2_Quit(SecNam,'bug detected','[1]')
      End If
      If (nOcc(iSymi).lt.1 .or. nVir(iSyma).lt.1) Then
         Call ChoMP2_Quit(SecNam,'bug detected','[2]')
      End If
#endif

C     Find a and i.
C     -------------

      i  = -999999
      a  = -999999
      i_ = 0
      Do While (i_ .lt. nOcc(iSymi))
         i_ = i_ + 1
         ai_1 = iT1Am(iSyma,iSymi) + nVir(iSyma)*(i_-1) + 1
         ai_2 = ai_1 + nVir(iSyma) - 1
         If (ai.ge.ai_1 .and. ai.le.ai_2) Then
            i = i_
            a = ai - ai_1 + 1
            i_ = nOcc(iSymi) + 1 ! Found! -- break loop
         End If
      End Do

#if defined (_DEBUG_)
      If (i.lt.1 .or. i.gt.nOcc(iSymi) .or.
     &    a.lt.1 .or. a.gt.nVir(iSyma)) Then
         Call ChoMP2_Quit(SecNam,'bug detected','[3]')
      End If
#endif

      End
