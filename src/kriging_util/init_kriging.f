************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Init_Kriging()
      use AI
      Implicit Real*8 (a-h,o-z)
C
C     Initiate Kriging parameters.
C
      Kriging = .False.
      nspAI = 3
      anMd = .True.
      anHe = .True.
      anGr = .True.
      anHCt = .True.
      numHt = 1.0D-3
      pAI = 2
      npxAI = 1
      lb(1) = 20.0D0
      lb(2) = 20.0D0
      lb(3) = 1
      miAI = 50
      meAI = 1.0D-8
      blAI = .False.
      blaAI = .True.
      blavAI=0.04D0
*
      nInter_save=-1
      nPoints_save=-1
*
      Return
      End
