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
      Subroutine RToCore(F,nBeta,ishll,lb,iAng,nveccb)
*******************************************************************************
*
* Transformation kernel to atomic orbials in normailized spherical harmonics
*
*******************************************************************************
*    @parameter F  The cartesian components of <A|core>
*    @parameter nBeta Number of exponents
*    @parameter ishll Shell number for ECP
*    @parameter lb angular momenta Ket
*    @parameter iAng angular momenta core
*    @parameter nveccb Number of derivatives
*******************************************************************************
*
      use Real_Spherical
      use Basis_Info
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "stdalloc.fh"
      Real*8 F(*)
      Real*8, Allocatable:: Tmp1(:), Tmp2(:)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

*******************************************************************************
      ncb=nelem(lb)*nelem(iang)
      nExpi=Shells(iShll)%nExp
      nBasisi=Shells(iShll)%nBasis
      Call mma_allocate(TMP1,nExpi*ncb*nVecCB*nBeta,Label='Tmp1')
      Call mma_allocate(TMP2,nExpi*ncb*nVecCB*nBeta,Label='Tmp2')
*
*--------------And (almost) the same thing for the righthand side, form
*              KjCb from kjcb
*              1) jcb,K = k,jcb * k,K
*
      Call DGEMM_('T','N',
     &            nBeta*ncb*nVecCB,nBasisi,nExpi,
     &            One,F,nExpi,
     &                Shells(iShll)%pCff,nExpi,
     &            Zero,Tmp1,nBeta*ncb*nVecCB)
*
*--------------2)  j,cbK -> cbK,j
*
      Call DgeTMo(Tmp1,nBeta,nBeta,
     &            ncb*nVecCB*nBasisi,Tmp2,
     &            ncb*nVecCB*nBasisi)
*
*--------------3) bKj,C = c,bKj * c,C
*
      Call DGEMM_('T','N',
     &            nElem(lb)*nVecCB*nBasisi*nBeta,
     &            (2*iAng+1),nElem(iAng),
     &            One,Tmp2,nElem(iAng),
     &                RSph(ipSph(iAng)),nElem(iAng),
     &            Zero,Tmp1,
     &            nElem(lb)*nVecCB*nBasisi*nBeta)
*
*--------------4) b,KjC -> KjC,b
*
      Call DgeTMo(Tmp1,nElem(lb)*nVecCB,
     &            nElem(lb)*nVecCB,
     &            nBasisi*nBeta*(2*iAng+1),F,
     &            nBasisi*nBeta*(2*iAng+1))
*
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp1)
      Return
      End
