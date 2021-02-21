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
* Copyright (C) 2009, Giovanni Ghigo                                   *
************************************************************************
C!-----------------------------------------------------------------------!
C!
C!    InterSystem Crossing rate evaluation: Reduction of States
C!    Author: Giovanni Ghigo
C!            Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
C!            07 Jan-09 - XX Jan-09
C!
      Subroutine LogEVec(iPrint,nOsc,max_nOrd,minQ,
     &                        nMaxQ,nMat,lVec,nYes)
C!
C!    Generate Logical Vector of usefull States
C!
      Implicit Real*8 ( a-h,o-z )
      Integer nMat(0:max_nOrd,nOsc), lVec(0:max_nOrd)
      Integer nMaxQ(nOsc)

      If (iPrint.GE.3) then
        Write(6,*) ' Original number of States=',max_nOrd+1
      EndIf

      Do iOrd = 0,max_nOrd
        lVec(iOrd) = 1
        nSumQ = 0
        Do iOsc=1,nOsc
          If (nMat(iOrd,iOsc).GT.nMaxQ(iOsc)) lVec(iOrd) = 0
          nSumQ = nSumQ + nMat(iOrd,iOsc)
        EndDo
        If (nSumQ.LT.minQ) lVec(iOrd) = 0
      EndDo
      nYes = 0
      Do iOrd = 0,max_nOrd
        nYes = nYes + lVec(iOrd)
      EndDo
C!
      If (iPrint.GE.3) then
        Write(6,*) ' Selected number of States=',nYes
      EndIf
C!
      Return
      End


      Subroutine MkVibWind2(iPrint,nYes,iMaxYes,max_nOrd,
     &                      lVec,VibWind2)
C!
      Implicit Real*8 ( a-h,o-z )
      Integer VibWind2(nYes), lVec(0:max_nOrd)

      iMaxYes = 0
      iYes = 1
      Do iOrd = 0, max_nOrd
        If (lVec(iOrd).EQ.1) then
          VibWind2(iYes) = iOrd
          iYes = iYes + 1
          iMaxYes = Max(iMaxYes,iOrd)
        EndIf
      EndDo
C!
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iPrint)
      End
