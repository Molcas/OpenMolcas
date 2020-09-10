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
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
      Real*8 F(*)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nac=nelem(la)*nelem(iang)
      nExpi=Shells(iShll)%nExp
      nBasisi=Shells(iShll)%nBasis
      Call Getmem('TMP1','ALLO','REAL',iptmp,
     &             nExpi*nac*nVecAC*nalpha)
      Call Getmem('TMP2','ALLO','REAL',ipF,
     &             nExpi*nac*nVecAC*nalpha)
*--------------From the lefthandside overlap, form iKaC from ikac by
*              1) i,kac -> k,aci
*
      n = nExpi*nac*nVecAC
      Call DgeTMo(F,nAlpha,
     &            nAlpha, n,
     &            Work(ipTmp),n)
*
*--------------2) aciK =  k,aci * k,K (Contract over core orbital)
*
      Call DGEMM_('T','N',
     &            nac*nVecAC*nAlpha,nBasisi,nExpi,
     &            1.0d0,Work(ipTmp),nExpi,
     &            Shells(iShll)%pCff,nExpi,
     &            0.0d0,Work(ipF),nac*nVecAC*nAlpha)
*
*--------------3) Mult by shiftoperators aci,K -> Bk(K) * aci,K
*
      Do iBk = 1, nBasisi
         Call DYaX(nac*nVecAC*nAlpha,Shells(iShll)%Bk(iBk),
     &              Work((iBk-1)*nac*nVecAC*nAlpha+ipF),1,
     &              Work((iBk-1)*nac*nVecAC*nAlpha+ipTmp),1)
      End Do
*
*--------------4) a,ciK -> ciKa
*
      Call DgeTMo(Work(ipTmp),nElem(la),nElem(la),
     &            nElem(iAng)*nVecAC*nAlpha*nBasisi,
     &            Work(ipF),
     &            nElem(iAng)*nVecAC*nAlpha*nBasisi)
*
*--------------5) iKa,C = c,iKa * c,C
*
      Call DGEMM_('T','N',
     &            nVecAC*nAlpha*nBasisi*nElem(la),
     &            (2*iAng+1),nElem(iAng),
     &            1.0d0,Work(ipF),nElem(iAng),
     &            RSph(ipSph(iAng)),nElem(iAng),
     &            0.0d0,Work(ipTmp),
     &            nVecAC*nAlpha*nBasisi*nElem(la))
*
      Call DgeTMo(Work(ipTmp),nVecAC,nVecAC,
     &            nAlpha*nBasisi*nElem(la)*(2*iAng+1),
     &            F,
     &            nAlpha*nBasisi*nElem(la)*(2*iAng+1))

      Call Getmem('TMP1','FREE','REAL',iptmp,
     &            nExpi*nac*nVecAC*nalpha)

      Call Getmem('TMP2','FREE','REAL',ipF,
     &            nExpi*nac*nVecAC*nalpha)
      Return
      End
