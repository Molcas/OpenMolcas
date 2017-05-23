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
      Subroutine Pickmo_td(rmo,rmoaa,idsym)
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
      real*8 rmo(*),rmoaa(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      If (.not.timedep) Then
      Do iS=1,nSym
       Do jS=1,iS
        Do kS=1,is
         ls=ieor(iEor(is-1,js-1),iEor(ks-1,idsym-1))+1
         If (ls.le.ks) Then
          Do iA=1,nAsh(is)
           iAA=iA+nA(is)
           Do jA=1,nAsh(js)
            jAA=jA+nA(js)
            ijAA=itri(iAA,jAA)
            Do kA=1,nAsh(ks)
             kAA=kA+nA(ks)
             Do lA=1,nAsh(ls)
              lAA=lA+nA(ls)
              klAA=iTri(kAA,lAA)
              If (ijAA.ge.klAA) Then
               ijkl=iTri(ijAA,klAA)
               ipi=ipMO(js,ks,ls)+nIsh(is)+iA-1+
     &          nBas(is)*(jA-1)+nBas(is)*nAsh(js)*
     6          (kA-1)+nBas(is)*nAsh(js)*nAsh(ks)*(lA-1)
*              ipi=ipMO(js,ks,ls)+
*    &          (nBas(is)*(jA-1+nAsh(js)*
*    6          (kA-1+nAsh(ks)*(lA-1))))+nish(is)+iA-1
               rmoaa(ijkl)=rmo(ipi)
              End If
             End Do
            End Do
           End Do
          End Do
         End If
        End Do
       End Do
      End Do
      Else
      Do iS=1,nSym
       Do jS=1,nsym
        Do kS=1,nsym
         ls=ieor(iEor(is-1,js-1),iEor(ks-1,idsym-1))+1
          Do iA=1,nAsh(is)
           iAA=iA+nA(is)
           Do jA=1,nAsh(js)
            jAA=jA+nA(js)
            ijAA=iAA+(jaa-1)*ntash
            Do kA=1,nAsh(ks)
             kAA=kA+nA(ks)
             Do lA=1,nAsh(ls)
              lAA=lA+nA(ls)
              klAA=kAA+(laa-1)*ntash
              If (ijAA.ge.klAA) Then
               ijkl=iTri(ijAA,klAA)
               ipi=ipMO(js,ks,ls)+nIsh(is)+iA-1+
     &          nBas(is)*(jA-1)+nBas(is)*nAsh(js)*
     6          (kA-1)+nBas(is)*nAsh(js)*nAsh(ks)*(lA-1)
*              ipi=ipMO(js,ks,ls)+
*    &          (nBas(is)*(jA-1+nAsh(js)*
*    6          (kA-1+nAsh(ks)*(lA-1))))+nish(is)+iA-1
               rmoaa(ijkl)=rmo(ipi)
              End If
             End Do
            End Do
           End Do
          End Do
        End Do
       End Do
      End Do
      End IF
      Return
      End
