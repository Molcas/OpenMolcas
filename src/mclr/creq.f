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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
       SubRoutine creq(q,rint,G2,idsym)
*
*      Constructs the Q matrix
*
       Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
       Real*8 Q(nDens2),rint(*),G2(*)
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*
*      Q = (pj|kl)d
*       pi         ijkl
*
       call dcopy_(ndens2,0.0d0,0,Q,1)
       Do iS=1,nSym
        ipS=iEOr(is-1,idsym-1)+1
         if (norb(ips).ne.0) Then

        Do jS=1,nsym
         ijS=iEOR(is-1,js-1)+1
         Do kS=1,nSym
          ls=iEOr(ijs-1,ks-1)+1
          Do iAsh=1,nAsh(is)
           Do jAsh=1,nAsh(js)
            iij=itri(iAsh+nA(is),jAsh+nA(jS))
            Do kAsh=1,nAsh(ks)
             Do lAsh=1,nAsh(ls)
              ikl=itri(lAsh+nA(lS),kAsh+nA(kS))
              ipQ=ipMatba(ips,is)+norb(ips)*(iAsh-1)
              ipG=itri(iij,ikl)
              ipi=ipMO(js,ks,ls)+
     &         (norb(ips)*(jAsh-1+nAsh(js)*
     6          (kAsh-1+nAsh(ks)*(lAsh-1))))
             call daxpy_(norb(ips),G2(ipG),rint(ipI),1,
     &                  Q(ipQ),1)
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do
        end if
       End Do
       Return
       end
