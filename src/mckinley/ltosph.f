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
      Subroutine LToSph(F,nalpha,ishll,la,iAng,nvecac)
*******************************************************************************
*
*   Transform <A|core> from Cartesian components to Sperical harmonics
*
*              Observe that as opposed to the projection operator that this
*              contraction is done in the primitive basis.
*
*******************************************************************************
*    @parameter F  The cartesian components of <A|core>(in)
*                  The spherical components of <A|core>(out)

*    @parameter nAlpha Number of exponents
*    @parameter ishll Shell number for ECP
*    @parameter la angular momenta LS
*    @parameter iAng angular momenta core
*    @parameter Number of derivatives
*******************************************************************************
*
      use Basis_Info
      use Real_Spherical
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"
      Dimension F(*)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*******************************************************************************
*
      nExpi=Shells(iShll)%nExp
      nac=nelem(la)*nelem(iang)
      Call Getmem('TMP1','ALLO','REAL',iptmp,
     &             nExpi*nac*nVecAC*nalpha)
      Call Getmem('TMP2','ALLO','REAL',ipF,
     &             nExpi*nac*nVecAC*nalpha)


      Call DgeTMo(F,nAlpha*nExpi*nElem(la),
     &            nAlpha*nExpi*nElem(la),
     &            nElem(iAng)*nVecAC,Work(ipTmp),
     &            nElem(iAng)*nVecAC)
*
*--------------2) xika,C = c,xika * c,C
*
      Call DGEMM_('T','N',
     &            nVecAC*nAlpha*nExpi*nElem(la),
     &            (2*iAng+1),nElem(iAng),
     &            1.0d0,Work(ipTmp),nElem(iAng),
     &            RSph(ipSph(iAng)),nElem(iAng),
     &            0.0d0,Work(ipF),
     &            nVecAC*nAlpha*nExpi*nElem(la))
*
*--------------3) x,ikaC -> ikaC,x
*
      Call DGetMo(Work(ipF),nVecAC,nVecAC,
     &            nAlpha*nExpi*nElem(la)*(2*iAng+1),
     &            Work(ipTmp),
     &            nAlpha*nExpi*nElem(la)*(2*iAng+1))
      call dcopy_(nVecAC*
     &           nAlpha*nExpi*nElem(la)*(2*iAng+1),
     &           Work(ipTmp),1,F,1)
*
      Call Getmem('TMP1','FREE','REAL',iptmp,
     &             nExpi*nac*nVecAC*nalpha)
      Call Getmem('TMP2','FREE','REAL',ipF,
     &             nExpi*nac*nVecAC*nalpha)
      Return
      End
