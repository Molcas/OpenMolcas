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
      Real*8 Function E2_td(FockI,rMo,loper,idisp)
      use Arrays, only: G1t, G2sq
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "disp_mclr.fh"

*
      Real*8 FockI(nCMO),rMO(*)
      Logical Go
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      E22=0.0d0
      If (loper.eq.0) Then
         Go = idisp.lt.0
         If (.Not.Go) Go = iAnd(ntpert(idisp),2**2).eq.4
         If (Go) Then
            Do i=1,nna
               Do j=1,nna
                  ij2=i+(j-1)*nna
                  ij=itri(i,j)
                  Do k=1,nna
                     Do l=1,nna
                        ijkl=itri(ij,itri(k,l))
                        kl2=k+(l-1)*nna
                        ijkl2=ij2+(kl2-1)*nna**2
                        E22=E22+0.5d0*G2sq(ijkl2)*rmo(ijkl)
                     End Do
                  End Do
               End Do
            End Do
         End If
         Do is=1,nSym
            Do iA=1,Nash(is)
               iAA=nA(is)+ia
               iAB=ia+nish(is)
               js=is
               Do jA=1,nAsh(js)
                  jAA=ja+nA(js)
                  jAB=jA+nIsh(js)
                  ipF=(iab-1)*nbas(is)+jab+ipCM(is)-1
                  E22=E22+Focki(ipf)*G1t(itri(iaa,jaa))
               End Do
            End Do
         End Do
      End If
*
      e2_td=e22
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
