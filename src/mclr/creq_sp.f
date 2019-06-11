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
       SubRoutine creq_sp(q,rint,G2,idsym)
*
*      Constructs the Q matrix
*
       Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
       Real*8 Q(nDens2),rint(*),G2(nna,nna,nna,nna)
*
*      Q = (pj|kl)d
*       pi         ijkl
*
*      call dcopy_(ndens2,[0.0d0],0,Q,1)
       Do iS=1,nSym
        ipS=iEOr(is-1,idsym-1)+1
         if (norb(ips).ne.0) Then

        Do jS=1,nsym
         ijS=iEOR(is-1,js-1)+1
         Do kS=1,nSym
          ls=iEOr(ijs-1,ks-1)+1
          Do iAsh=1,nAsh(is)
           Do jAsh=1,nAsh(js)
            Do kAsh=1,nAsh(ks)
             Do lAsh=1,nAsh(ls)
              ipQ=ipMat(ips,is)+norb(ips)*(nish(is)+iAsh-1)
              ipi=ipMO(js,ks,ls)+
     &         (norb(ips)*(jAsh-1+nAsh(js)*
     6          (kAsh-1+nAsh(ks)*(lAsh-1))))
               rd=G2(na(is)+iash,na(js)+jash,na(ks)+kash,na(ls)+lash)
             call daxpy_(norb(ips),rd,rint(ipI),1,
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
