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
      Subroutine Get_D1ao_ab(ipD1ao,nDens)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "SysDef.fh"

      Character*24 Label
      Logical      Found

      Call Get_iScalar('System BitSwitch',iOption)
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*                                                                      *
************************************************************************
*                                                                      *
      Label='D1ao_ab'
      Call qpg_dArray(Label,Found,nDens)
      If (.not.Found .or.nDens.eq.0) Then
         Call SysAbendMsg('get_d1ao_ab','Could not locate:',Label)
      End If
      Call GetMem('Dens_ab','Allo','Real',ipD1ao,nDens)
      Call get_dArray(Label,Work(ipD1ao),nDens)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
