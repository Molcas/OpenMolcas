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
      Subroutine ConvdFdRho(mGrid,dF_dRho,ndF_dRho,RhoI,RhoA,mRho)
      Implicit Real*8 (A-H,O-Z)
      Dimension dF_dRho(ndF_dRho,mGrid)
      Dimension RhoI(mRho,mGrid),RhoA(mRho,mGrid)
*
      If (mRho.eq.1) Then
        Do iGrid=1,mGrid
         dFdRho = dF_dRho(1,iGrid)
         dFdP2  = dF_dRho(2,iGrid)
*
         dF_dRho(1,iGrid)= dFdRho + RhoI(1,iGrid)*dFdP2
         dF_dRho(2,iGrid)= RhoA(1,iGrid)*dFdP2/2.0d0
        End Do
      Else If (mRho.eq.4) Then
        Do iGrid=1,mGrid
         dFdRho = dF_dRho(1,iGrid)
         dFdRhox= dF_dRho(3,iGrid)
         dFdRhoy= dF_dRho(5,iGrid)
         dFdRhoz= dF_dRho(7,iGrid)
         dFdP2  = dF_dRho(2,iGrid)
         dFdP2x = dF_dRho(4,iGrid)
         dFdP2y = dF_dRho(6,iGrid)
         dFdP2z = dF_dRho(8,iGrid)
*
         dF_dRho(1,iGrid)=dFdRho +
     &        RhoI(1,iGrid)*dFdP2  +
     &       (RhoI(2,iGrid)*dFdP2x + RhoI(3,iGrid)*dFdP2y +
     &        RhoI(4,iGrid)*dFdP2z)*2.0d0*1.0d0
         dF_dRho(3,iGrid)=dFdRhox + RhoI(1,iGrid)*dFdP2x
         dF_dRho(5,iGrid)=dFdRhoy + RhoI(1,iGrid)*dFdP2y
         dF_dRho(7,iGrid)=dFdRhoz + RhoI(1,iGrid)*dFdP2z
*
         dF_dRho(2,iGrid)=
     &       RhoA(1,iGrid)*dFdP2/2.0d0 +
     &       (RhoA(2,iGrid)*dFdP2x + RhoA(3,iGrid)*dFdP2y +
     &       RhoA(4,iGrid)*dFdP2z)*1.0d0
         dF_dRho(4,iGrid)=RhoA(1,iGrid)*dFdP2x/2.0d0
         dF_dRho(6,iGrid)=RhoA(1,iGrid)*dFdP2y/2.0d0
         dF_dRho(8,iGrid)=RhoA(1,iGrid)*dFdP2z/2.0d0
        End Do
      Else
        Call WarningMessage(2,'Somethings is wrong in ConvdFdRho')
         Call Abend()
      End If
*
      Return
      End
