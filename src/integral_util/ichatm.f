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
      Function iChAtm(Coor,iOper,nOper,iChCar)
      Implicit Real*8 (a-h,o-z)
      Integer iChAtm, iOper(0:7), iChCar(3)
      Real*8 Coor(3)
*
      iChAtm=0
      Do 2000 kxyz = 1, 3
*        Test if component is not zero
         If (Abs(Coor(kxyz)).lt.1.D-12) Go to 2000
*------- Loop over the group generators
         Do 2001 i = 1, nOper
            j = i
            If (i.eq.3) j = 4
*           Test if symoperation will permute component
            If (iAnd(iOper(j),iChCar(kxyz)).ne.0) Go To 2002
 2001    Continue
         Go To 2000
 2002    iChAtm = iChAtm + 2**(kxyz-1)
 2000 Continue
*
      Return
      End
