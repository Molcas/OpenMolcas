!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Integer Function iChxyz(Coord,iGen,nGen)
      use Constants
      Implicit None
      Integer nGen
      Real*8 Coord(3)
      Integer iGen(nGen), iChCar(3)

      Integer iCar
!
      Call ChCar(iChCar,iGen,nGen)
      iChxyz=0
      Do iCar = 1, 3
         If (Coord(iCar).ne.Zero) iChxyz = iChxyz + iChCar(iCar)
      End Do
!
      Return
      End Function iChxyz
