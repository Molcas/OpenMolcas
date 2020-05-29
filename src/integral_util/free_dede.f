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
      Subroutine Free_DeDe(Dens,TwoHam,nDens,ipDq,ipFq)
      use k2_arrays
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "setup.fh"
*
      Real*8 Dens(nDens), TwoHam(nDens)
*
      If (nIrrep.eq.1) Then
* symmetrize fock matrix
*.... Fix the diagonal elements of D and F
         Call DScal_(nDens,Two,Dens,1)
         nc=nbas(0)
         mDens=nBas(0)**2
         ijq=ipFq-1
         jiq=ipFq-nc
         ij=0
         do i=1,nc
           do j=1,i
             ij=ij+1
             TwoHam(ij) = Half*(Work(ijq+j) + Work(jiq+j*nc))
           end do
           Dens(ij)  =Half*Dens(ij)
           jiq = jiq + 1
           ijq = ijq + nc
         end do
         Call GetMem('FMAQ','Free','Real',ipFq,mDens)
         Call GetMem('DENQ','Free','Real',ipDq,mDens)
      End If
*
      Call mma_deallocate(ipOffD)
      Call GetMem('DeDe2','Free','Real',ipDeDe,nDeDe)
*
      Return
      End
