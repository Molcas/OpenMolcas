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
      Real*8 Function E2(FockI,rMo,loper,idisp)
      use Arrays, only: G1t
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "WrkSpc.fh"

      Real*8 FockI(nCMO),rMO(*)
      Logical Go
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
      E22=0.0d0
      If (loper.eq.0) Then
         Go = iDisp.lt.0
         If (.Not.Go) Go = iAnd(ntpert(idisp),2**2).eq.4
         If (Go) Then
            Do i=1,nna
               Do j=1,nna
                  ij=itri(i,j)
                  Do k=1,nna
                     Do l=1,nna
                        ijkl=itri(ij,itri(k,l))
                        E22=E22+0.5d0*work(ipg2+ijkl-1)*rmo(ijkl)
                     End Do
                  End Do
               End Do
            End Do
         End If
         Do is=1,nSym
            Do iA=1,nAsh(is)
               iAA=nA(iS)+ia
               iAB=ia+nIsh(iS)
               js=is
               Do jA=1,nAsh(js)
                  jAA=ja+nA(js)
                  jAB=jA+nIsh(js)
                  ipF=(iab-1)*norb(is)+jab+ipCM(is)-1
                  E22=E22+Focki(ipf)*G1t(itri(iaa,jaa))
               End Do
            End Do
         End Do
      End If
*
      e2=e22
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
