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
      use k2_arrays, only: ipD00, MaxDe
#include "WrkSpc.fh"
      Real*8 FD(nFD), DDen(lDDen)
      Integer ipOffD(nOffD)
      Logical Special_NoSym, DFT_Storage
      Call QEnter('DeDe_Mck')
*
      Special_NoSym=.False.
      DFT_Storage=.False.
      nr_of_Densities=1
      jpDeDe=ip_of_Work(DDen(1))
      ipD00=jpDeDe
      Call DeDe(FD,nFD,nr_of_Densities,ipOffD,nOffD,jpDeDe,ipD00,MaxDe,
     &          mDeDe,mIndij,Special_NoSym,DFT_Storage,Work,1)
*
      Call QExit('DeDe_Mck ')
      Return
      End
