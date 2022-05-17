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
      Subroutine RToSph(F,nBeta,ishll,lb,iAng,nveccb)
*******************************************************************************
*
*   Transform <core|B> from Cartesian components to Sperical harmonics
*
*******************************************************************************
*    @parameter F  The cartesian components of <core|B>(in)
*                  The spherical components of <core|B>(out)
*    @parameter nBeta Number of exponents
*    @parameter ishll Shell number for ECP
*    @parameter lb angular momenta Ket
*    @parameter iAng angular momenta core
*    @parameter Number of derivatives
*******************************************************************************
      use Real_Spherical
      use Basis_Info, only: Shells
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "stdalloc.fh"
      Dimension F(*)
      Real*8, Allocatable:: Tmp1(:), Tmp2(:)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

*******************************************************************************

      ncb=nelem(lb)*nelem(iang)
      nExpi=Shells(iShll)%nExp
      Call mma_allocate(TMP1,nExpi*ncb*nVecCB*nBeta,Label='Tmp1')
      Call mma_allocate(TMP2,nExpi*ncb*nVecCB*nBeta,Label='Tmp2')


*-------------1) kj,cbx -> cbx,kj
*
      Call DgeTMo(F,
     &            nBeta*nExpi,nBeta*nExpi,
     &            ncb*nVecCB,Tmp1,ncb*nVecCB)
*
*--------------2) bxkj,C = c,bxkj * c,C
*
      Call DGEMM_('T','N',
     &            nElem(lb)*nVecCB*nExpi*nBeta,
     &            (2*iAng+1),nElem(iAng),
     &            One,Tmp1,nElem(iAng),
     &                RSph(ipSph(iAng)),nElem(iAng),
     &            Zero,Tmp2,nElem(lb)*nVecCB*nExpi*nBeta)
*
*--------------3) bx,kjC -> kjC,bx
*
      Call DgeTMo(Tmp2,nElem(lb)*nVecCB,
     &            nElem(lb)*nVecCB,
     &            nExpi*nBeta*(2*iAng+1),Tmp1,
     &            nExpi*nBeta*(2*iAng+1))

      call dcopy_(nExpi*
     &           nBeta*(2*iAng+1)*nElem(lb)*nVecCB,
     &           Tmp1,1,F,1)

      Call mma_deallocate(Tmp1)
      Call mma_deallocate(Tmp2)
      Return
      End
