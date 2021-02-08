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
      Subroutine Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)
      Use Basis_Info
      Implicit Real*8 (a-h,o-z)
*
      Integer iSD4(0:nSD,4)
*
      Real*8  Coor(3,4)
      Integer iAOV(4), iStabs(4)
      Logical  Shijij
*
      iCnttp=iSD4(13,1)
      iCnt  =iSD4(14,1)
      jCnttp=iSD4(13,2)
      jCnt  =iSD4(14,2)
      kCnttp=iSD4(13,3)
      kCnt  =iSD4(14,3)
      lCnttp=iSD4(13,4)
      lCnt  =iSD4(14,4)
*
      If (dbsc(iCnttp)%Aux) Then
         Coor(1:3,1)=dbsc(jCnttp)%Coor(1:3,jCnt)
      Else
         Coor(1:3,1)=dbsc(iCnttp)%Coor(1:3,iCnt)
      End If
      Coor(1:3,2)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
      If (dbsc(kCnttp)%Aux) Then
         Coor(1:3,3)=dbsc(lCnttp)%Coor(1:3,lCnt)
      Else
         Coor(1:3,3)=dbsc(kCnttp)%Coor(1:3,kCnt)
      End If
      Coor(1:3,4)=dbsc(lCnttp)%Coor(1:3,lCnt)
*
      Shijij= iSD4(11,1).eq.iSD4(11,3).and.
     &        iSD4(11,2).eq.iSD4(11,4)
*
      Do iQuad = 1, 4
         iAOV(iQuad)   = iSD4( 7,iQuad)
         iStabs(iQuad) = iSD4(10,iQuad)
      End Do
*
      Return
      End
