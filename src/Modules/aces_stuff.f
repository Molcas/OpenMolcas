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
      Module Aces_Stuff
*
*---- Stuff for Aces 2 read of Gamma file
*
      Logical                Gamma_On
      Integer                LuGamma, lBin, nGamma
      Real*8, Allocatable:: G_Toc(:), Bin(:,:)
      Integer, Allocatable:: SO2cI(:,:)
*
      End Module Aces_Stuff
