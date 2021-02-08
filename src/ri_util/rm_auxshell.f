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
      Subroutine rm_AuxShell(iCnttp)
************************************************************************
*                                                                      *
*     Remove an auxiliary basis set by making it empty.                *
*                                                                      *
************************************************************************
      Use Basis_Info, only: dbsc, Shells
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
#include "real.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Do k = 0, dbsc(iCnttp)%nShells-1
         iShll = dbsc(iCnttp)%iVal + k
*
         Shells(iShll)%nExp=0
         Shells(iShll)%nBasis  =0
         Shells(iShll)%nBasis_c=0
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
