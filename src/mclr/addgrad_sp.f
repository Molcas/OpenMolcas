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
      SubRoutine AddGrad_sp(rKappa,rMat,F,idsym,r1,r2)
*
*     Purpose:
*             Adds the contribution from the gradient to
*              [2]
*             E   Kappa. This is done to insure us about
*             a beautifull convergence of the PCG,
*             which is just the case if E is symmetric.
*
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Pointers.fh"

#include "Input.fh"
#include "stdalloc.fh"
      Real*8 rkappa(*),rMat(*),F(*)
      Real*8, Allocatable:: Tempi(:), Tempj(:)
      Real*8, Allocatable:: K(:), M(:)

      Call mma_allocate(K,ndens2,Label='K')
      Call mma_allocate(M,ndens2,Label='M')
      M(:)=Zero
      Call Unc(rkappa,K,idsym,r1)

      Do iS=1,nSym
        js=iEor(is-1,idsym-1)+1
        If (nOrb(is)*nOrb(js)==0) Cycle
        Call mma_allocate(Tempi,nOrb(is)**2,Label='Tempi')
        Call mma_allocate(Tempj,nOrb(js)**2,Label='Tempj')
*
*    T=Brillouin matrix
*

        Call DGeSub(F(ipCM(is)),nOrb(is),'N',
     &              F(ipCM(is)),nOrb(is),'T',
     &              Tempi,nOrb(is),
     &              nOrb(is),nOrb(is))
        Call DGeSub(F(ipCM(js)),nOrb(js),'N',
     &              F(ipCM(js)),nOrb(js),'T',
     &              Tempj,nOrb(js),
     &              nOrb(js),nOrb(js))
*
*               t             t
*        { Kappa T  r2 T kappa  }
*
*
        Call DGEMM_('T','N',nOrb(is),nOrb(js),nOrb(js),
     &             One,K(ipMat(js,is)),nOrb(js),
     &                 Tempj,nOrb(js),
     &             Zero,M(ipMat(is,js)),nOrb(is))
        Call DGEMM_('N','T',nOrb(is),nOrb(js),nOrb(is),
     &             r2,Tempi,nOrb(is),
     &                K(ipmat(js,is)),nOrb(js),
     &             One,M(ipMat(is,js)),nOrb(is))
        Call mma_deallocate(Tempi)
        Call mma_deallocate(Tempj)
      End Do
      Call COMPRESS(M,rmat,idsym)
      Call mma_deallocate(K)
      Call mma_deallocate(M)
      Return
      End
