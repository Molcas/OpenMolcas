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
      Subroutine LToCore(F,nalpha,ishll,la,iAng,nvecac)
*******************************************************************************
*
* Transformation kernel to atomic orbials in normailized spherical harmonics
*
*******************************************************************************
*    @parameter F  The cartesian components of <A|core>
*    @parameter nAlpha Number of exponents
*    @parameter ishll Shell number for ECP
*    @parameter la angular momenta LS
*    @parameter iAng angular momenta core
*    @parameter Number of derivatives
*
      use Real_Spherical
      use Basis_Info
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 F(*)
      Real*8, Allocatable:: Tmp1(:), Tmp2(:)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nac=nelem(la)*nelem(iang)
      nExpi=Shells(iShll)%nExp
      nBasisi=Shells(iShll)%nBasis
      Call mma_allocate(Tmp1,nExpi*nac*nVecAC*nalpha,Label='Tmp1')
      Call mma_allocate(Tmp2,nExpi*nac*nVecAC*nalpha,Label='Tmp2')
*--------------From the lefthandside overlap, form iKaC from ikac by
*              1) i,kac -> k,aci
*
      n = nExpi*nac*nVecAC
      Call DgeTMo(F,nAlpha,
     &            nAlpha, n,
     &            Tmp1,n)
*
*--------------2) aciK =  k,aci * k,K (Contract over core orbital)
*
      Call DGEMM_('T','N',
     &            nac*nVecAC*nAlpha,nBasisi,nExpi,
     &            One,Tmp1,nExpi,
     &                Shells(iShll)%pCff,nExpi,
     &            Zero,Tmp2,nac*nVecAC*nAlpha)
*
*--------------3) Mult by shiftoperators aci,K -> Bk(K) * aci,K
*
      Do iBk = 1, nBasisi
         Call DYaX(nac*nVecAC*nAlpha,Shells(iShll)%Bk(iBk),
     &              Tmp2((iBk-1)*nac*nVecAC*nAlpha+1),1,
     &              Tmp1((iBk-1)*nac*nVecAC*nAlpha+1),1)
      End Do
*
*--------------4) a,ciK -> ciKa
*
      Call DgeTMo(Tmp1,nElem(la),nElem(la),
     &            nElem(iAng)*nVecAC*nAlpha*nBasisi,
     &            Tmp2,
     &            nElem(iAng)*nVecAC*nAlpha*nBasisi)
*
*--------------5) iKa,C = c,iKa * c,C
*
      Call DGEMM_('T','N',
     &            nVecAC*nAlpha*nBasisi*nElem(la),
     &            (2*iAng+1),nElem(iAng),
     &            One,Tmp2,nElem(iAng),
     &                 RSph(ipSph(iAng)),nElem(iAng),
     &            Zero,Tmp1,nVecAC*nAlpha*nBasisi*nElem(la))
*
      Call DgeTMo(Tmp1,nVecAC,nVecAC,
     &            nAlpha*nBasisi*nElem(la)*(2*iAng+1),
     &            F,
     &            nAlpha*nBasisi*nElem(la)*(2*iAng+1))

      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp1)

      Return
      End
