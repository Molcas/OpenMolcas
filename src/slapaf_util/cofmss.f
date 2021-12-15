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
      Subroutine CofMss(Coor,dMass,nsAtom,LWrite,cMass,iSym)
************************************************************************
*     Object: To calculate the molecular mass, the center of mass and  *
*             move the coordinates so origo is the center of mass.     *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 COOR(3,nsAtom), dMass(nsAtom), cMass(3)
      Integer iSym(3)
      Logical LWRITE
*
      Return
*
*     Calculate the molecular mass.
*
      TMass = Zero
      Do I = 1, nsAtom
         TMass = TMass + dMass(I) * DBLE(iDeg(Coor(1,i)))
      End Do
      iCOM=-1
      If (TMass.ge.1.D99) Then
         Do i = 1, nsAtom
            If (dMass(i).eq.1.D99) Then
               iCOM=i
               Go To 99
            End If
         End Do
      End If
 99   Continue
*
*     calculate the center of mass
*
      call dcopy_(3,[Zero],0,cMass,1)
*-----Loop over the unique centers
      Do i = 1, nsAtom
         Do j = 1, 3
*-----------Add contribution
            If (iSym(j).eq.0) cMass(j) = cMass(j) +
     &         dMass(i) *  Coor(j,i) *
     &         DBLE(iDeg(Coor(1,i)))
         End Do
      End Do
*
      Do i = 1, 3
         cMass(i) = cMass(i) / TMass
      End Do
      If (iCOM.ge.1.and.iCOM.le.nsAtom)
     &   call dcopy_(3,COOR(1,iCOM),1,cMass,1)
*
      If (LWrite) Write(6,100) (cMass(i),i=1,3), TMass
 100  FORMAT(//,' Center of Mass (Bohr) ',3F10.5,/,
     &          ' Molecular Mass   (au) ',1F15.5)
*
*     translate the center of mass to origo
*
      Do i = 1, nsAtom
         Do j = 1, 3
            Coor(j,i) = Coor(j,i) - cMass(j)
         End Do
      End Do
*
      Return
      End
