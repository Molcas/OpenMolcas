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
      Subroutine DeDe_Mck(FD,nFD,ipOffD,nOffD,DDen,lDDen,mDeDe,mIndij)
      use k2_arrays, only: MaxDe
      Real*8 FD(nFD), DDen(lDDen)
      Integer ipOffD(nOffD)
      Logical Special_NoSym, DFT_Storage
*
      Special_NoSym=.False.
      DFT_Storage=.False.
      nr_of_Densities=1
*
      ipDeDe=1
      ipD00=1
!     ipDijS is controlled in the calling routine
      Call mk_DeDe(FD,nFD,nr_of_Densities,ipOffD,nOffD,ipDeDe,ipD00,
     &             MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,
     &             DDen,lDDen)
*
      Return
      End
