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
      Subroutine Free_DeDe_Funi(Dens,nDens,ipDq)
      use k2_arrays
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "setup.fh"
      Real*8 Dens(nDens)
*
*     If (nIrrep.eq.1) Then
*.... Fix the diagonal elements of D
*        Call DScal_(nDens,Two,Dens,1)
*        ij=0
*        Do i=1,nBas(0)
*          ij=ij+i
*          Dens(ij)  =Half*Dens(ij)
*        End Do
*        mDens=nBas(0)**2
*        Call GetMem('DENQ','Free','Real',ipDq,mDens)
*     End If
*
      Call mma_deallocate(ipOffD)
      Call GetMem('DeDe2','Free','Real',ipDeDe,nDeDe_DFT)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Dens)
         Call Unused_integer(ipDq)
      End If
      End
