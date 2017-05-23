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
      Module vRys_RW
      Integer MaxRys
      Parameter(MaxRys=9)
      Real*8, Dimension(:), Allocatable :: TMax, ddx, x0
      Real*8, Dimension(:), Allocatable :: HerR2, HerW2, Cff
      Integer, Dimension(:), Allocatable :: Map
      Integer, Dimension(:), Allocatable :: iHerR2, iHerW2
      Integer   iMap(MaxRys), nMap(MaxRys), ix0(MaxRys), nx0(MaxRys),
     &          iCffR(0:6,MaxRys), iCffW(0:6,MaxRys), nMxRys
      End Module vRys_RW
