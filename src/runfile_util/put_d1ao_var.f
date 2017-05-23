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
      Subroutine Put_D1ao_Var(D1ao,nD1ao)
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
      Real*8       D1ao(nD1ao)
      Character*24 Label
*
*     Write the variational 1st order density matrix
*     density matrix in AO/SO basis
*
      Label='D1aoVar'
      Call Put_dArray(Label,D1ao,nD1ao)

      Return
      End
