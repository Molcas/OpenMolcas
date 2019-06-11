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
       SubRoutine creq2(q,G2,idSym,Temp,Scr,n2)
*
*      Constructs the Q matrix
*
       Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
       Real*8 Q(nDens2),G2(*),Temp(n2),Scr(n2)
*                                                                      *
************************************************************************
*                                                                      *
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*      Q = (pj|kl)d
*       pi         ijkl
*
       call dcopy_(ndens2,[0.0d0],0,Q,1)
*
       Do iS=1,nSym
          ipS=iEOr(is-1,idsym-1)+1
          If (norb(ips).ne.0) Then
             Do jS=1,nsym
                ijS=iEOR(is-1,js-1)+1
                Do kS=1,nSym
                   ls=iEOr(ijs-1,ks-1)+1
*
                   Do kAsh=1,nAsh(ks)
                      Do lAsh=1,nAsh(ls)
                         ikl=itri(lAsh+nA(lS),kAsh+nA(kS))
*
                         Call Coul(ipS,jS,kS,lS,
     &                             nIsh(kS)+kAsh,nIsh(lS)+lAsh,Temp,Scr)
*
                         Do iAsh=1,nAsh(is)
                            ipQ=ipMatba(ipS,iS)+nOrb(ipS)*(iAsh-1)
                            Do jAsh=1,nAsh(jS)
                               iij=iTri(iAsh+nA(iS),jAsh+nA(jS))
                               ipG=iTri(iij,ikl)
                               ipI = (nIsh(jS)+jAsh-1)*nOrb(ipS) + 1
*
                              call daxpy_(nOrb(ipS),G2(ipG),Temp(ipI),1,
     &                                   Q(ipQ),1)
*
                            End Do
                         End Do
                      End Do
                   End Do
*
                End Do
             End Do
          End If
       End Do
*                                                                      *
************************************************************************
*                                                                      *
       Return
       End
