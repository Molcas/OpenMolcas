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
      Function iChAtm(Coor,iChCar)
      Use Symmetry_Info, only: nIrrep, iOper
      Implicit None
      Integer iChAtm, iChCar(3)
      Real*8 Coor(3)
      Integer iCar, i, j, nOper
      If (nIrrep.eq.8) Then
         nOper=3
      Else If (nIrrep.eq.4) Then
         nOper=2
      Else If (nIrrep.eq.2) Then
         nOper=1
      Else
         nOper=0
      End If

*
*     iChAtm is an integer function which will return an integer such
*     that the three first bits will represent the characteristics of
*     the Cartesian components. If the bit is set then the Cartesian
*     component will change sign if the symmetry operator contains a
*     part which operates on that particular Cartesian direction.
*
!     Default that none of the Cartesians will change sign.
      iChAtm=0
!
!     Loop over the Cartesian components.
!
      Do iCar = 1, 3
!
*        Test if component is not zero. If zero no operator will change
!        the sign of the component.
         If (Abs(Coor(iCar)).lt.1.D-12) Cycle
*
*        Here if the Component is none zero.
*
*------- Loop over the group generators and check if there is an
*        operator that will change the sign.
*
*        The generators are stored in positions 1, 2, and 4.
*
         Do i = 1, nOper    ! skip the unit operator -- i=0
            j = 2**(i-1)
*
*           Test if symoperation will permute component
*
            If (iAnd(iOper(j),iChCar(iCar)).ne.0) Then
               iChAtm = iChAtm + 2**(iCar-1)
               Exit
            End If
         End Do
      End Do
*
      Return
      End
