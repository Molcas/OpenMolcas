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
      Subroutine Get_Density_Matrix_mpprop(ip_D,nDens,nBas,nSym)

      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Integer nBas(nSym)

      If (nSym.eq.1) Then
       Call Get_D1ao(ip_D,nDens)
#ifdef _DEBUGPRINT_
       Call RecPrt('D',' ',Work(ip_D),1,nDens)
#endif
      Else
       Write(6,*) 'MpProp cannot handle symmetry'
       Call Abend()
      EndIf
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nBas)
      End
*
