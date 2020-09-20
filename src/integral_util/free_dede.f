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
      Subroutine Free_DeDe(Dens,TwoHam,nDens)
      use k2_arrays
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "setup.fh"
*
      Real*8, Target:: Dens(nDens), TwoHam(nDens)
*
      Nullify(pDq)
      Nullify(pFq)
*
      If (nIrrep.eq.1) Then
* symmetrize fock matrix
*.... Fix the diagonal elements of D and F
         Call DScal_(nDens,Two,Dens,1)
         nc=nbas(0)
         mDens=nBas(0)**2
         ijq=0
         jiq=1-nc
         ij=0
         do i=1,nc
           do j=1,i
             ij=ij+1
             TwoHam(ij) = Half*(Fq(ijq+j) + Fq(jiq+j*nc))
           end do
           Dens(ij)  =Half*Dens(ij)
           jiq = jiq + 1
           ijq = ijq + nc
         end do
         Call mma_deallocate(Dq)
         Call mma_deallocate(Fq)
      End If
*
      Call mma_deallocate(ipOffD)
      Call mma_deallocate(DeDe)
*
      Return
      End
