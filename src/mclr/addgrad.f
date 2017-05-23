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
      SubRoutine AddGrad(rKappa,rMat,idsym,fact)
*
*     Purpose:
*             Adds the contribution from the gradient to
*              [2]
*             E   Kappa. This is done to insure us about
*             a beautifull convergence of the PCG,
*             which is just the case if E is symmetric.
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
      Real*8 rkappa(*),rMat(*)
      Do iS=1,nSym
        js=iEor(is-1,idsym-1)+1
        If (nOrb(is)*nOrb(js).gt.0) Then
        Call Getmem('OJTEMP','ALLO','REAL',ipTempi,nOrb(is)**2)
        Call Getmem('OJTEMP','ALLO','REAL',ipTempj,nOrb(js)**2)
*
*    T=Brillouin matrix
*

        Call DGeSub(Work(ipF0SQMO+ipCM(is)-1),nOrb(is),'N',
     &              Work(ipF0SQMO+ipCM(is)-1),nOrb(is),'T',
     &              Work(ipTempi),nOrb(is),
     &              nOrb(is),nOrb(is))
        Call DGeSub(Work(ipF0SQMO+ipCM(js)-1),nOrb(js),'N',
     &              Work(ipF0SQMO+ipCM(js)-1),nOrb(js),'T',
     &              Work(ipTempj),nOrb(js),
     &              nOrb(js),nOrb(js))
*
*               t           t
*   +1/2 { Kappa T - T kappa  }
*
*
        Call DGEMM_('T','N',nOrb(is),nOrb(js),nOrb(js),
     &             0.5d0*fact,rkappa(ipMat(js,is)),nOrb(js),
     &             Work(ipTempj),nOrb(js),1.0d0,
     &             rMat(ipMat(is,js)),nOrb(is))
        Call DGEMM_('N','T',nOrb(is),nOrb(js),nOrb(is),
     &             -0.5d0*fact,Work(ipTempi),nOrb(is),
     &             rKappa(ipmat(js,is)),nOrb(js),1.0d0,
     &             rMat(ipMat(is,js)),nOrb(is))
        Call Getmem('OJTEMP','FREE','REAL',ipTempi,nOrb(is)**2)
        Call Getmem('OJTEMP','FREE','REAL',ipTempj,nOrb(is)**2)
      End If
      End Do
      Return
      End
