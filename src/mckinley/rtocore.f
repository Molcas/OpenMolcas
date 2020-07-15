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

#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
      Real*8 F(*)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

*******************************************************************************
      ncb=nelem(lb)*nelem(iang)
      nExpi=Shells(iShll)%nExp
      Call Getmem('TMP1','ALLO','REAL',iptmp,
     &             nExpi*ncb*nVecCB*nBeta)
      Call Getmem('TMP2','ALLO','REAL',ipF,
     &             nExpi*ncb*nVecCB*nBeta)
*
*--------------And (almost) the same thing for the righthand side, form
*              KjCb from kjcb
*              1) jcb,K = k,jcb * k,K
*
      Call DGEMM_('T','N',
     &            nBeta*ncb*nVecCB,nBasis(iShll),nExpi,
     &            1.0d0,F,nExpi,
     &            Shells(iShll)%pCff,nExpi,
     &            0.0d0,Work(ipTmp),nBeta*ncb*nVecCB)
*
*--------------2)  j,cbK -> cbK,j
*
      Call DgeTMo(Work(ipTmp),nBeta,nBeta,
     &            ncb*nVecCB*nBasis(iShll),Work(ipF),
     &            ncb*nVecCB*nBasis(iShll))
*
*--------------3) bKj,C = c,bKj * c,C
*
      Call DGEMM_('T','N',
     &            nElem(lb)*nVecCB*nBasis(iShll)*nBeta,
     &            (2*iAng+1),nElem(iAng),
     &            1.0d0,Work(ipF),nElem(iAng),
     &            RSph(ipSph(iAng)),nElem(iAng),
     &            0.0d0,Work(ipTmp),
     &            nElem(lb)*nVecCB*nBasis(iShll)*nBeta)
*
*--------------4) b,KjC -> KjC,b
*
      Call DgeTMo(Work(ipTmp),nElem(lb)*nVecCB,
     &            nElem(lb)*nVecCB,
     &            nBasis(iShll)*nBeta*(2*iAng+1),F,
     &            nBasis(iShll)*nBeta*(2*iAng+1))
*
      Call Getmem('TMP1','FREE','REAL',iptmp,
     &            nExpi*ncb*nVecCB*nBeta)
      Call Getmem('TMP2','FREE','REAL',ipF,
     &            nExpi*ncb*nVecCB*nBeta)
       Return
       End
