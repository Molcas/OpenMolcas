************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      SubRoutine AddGrad2(rMat,idsym,fact)
*
*     Purpose:
*             Adds the contribution from the gradient to
*              [2]
*             E   Kappa. This is done to insure us about
*             a beautiful convergence of the PCG,
*             which is just the case if E is symmetric.
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 rMat(*)
      Real*8, Allocatable:: Temp(:)

      Do iS=1,nSym
        If (nOrb(is)*nOrb(is)==0) Cycle
        Call mma_allocate(Temp,nOrb(is)**2,Label='Temp')
*
*    T=Brillouin matrix
*

        Call DGeSub(Work(ipF0SQMO+ipCM(is)-1),nOrb(is),'N',
     &              Work(ipF0SQMO+ipCM(is)-1),nOrb(is),'T',
     &              Temp,nOrb(is),
     &              nOrb(is),nOrb(is))
*
*               t           t
*   +1/2 { Kappa T - T kappa  }
*
*
        Call DaXpY_(nOrb(is)**2,-2.0d0*Fact,Temp,1,
     &              rMat(ipMat(is,is)),1)
        Call mma_deallocate(Temp)
      End Do

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(idsym)
      End
