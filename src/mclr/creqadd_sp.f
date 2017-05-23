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
       SubRoutine creqadd_sp(q,G2,idsym,Temp,Scr,n2)
*
*      Constructs the Q matrix
*
       Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
       Real*8 Q(nDens2),G2(nna,nna,nna,nna), Temp(n2),Scr(n2)
*                                                                      *
************************************************************************
*                                                                      *
*      Q = (pj|kl)d
*       pi         ijkl
*
       Do iS=1,nSym
          ipS=iEOr(is-1,idsym-1)+1
          If (norb(ips).ne.0) Then
             Do jS=1,nsym
             ijS=iEOR(is-1,js-1)+1
             Do kS=1,nSym
                ls=iEOr(ijs-1,iEor(ks-1,idsym-1))+1
*
                Do kAsh=1,nAsh(ks)
                   kAA=kAsh+nIsh(ks)
                   Do lAsh=1,nAsh(ls)
                      lAA=lAsh+nIsh(ls)
*
                      Call Coul(ipS,jS,kS,lS,kAA,lAA,Temp,Scr)
*
                      Do iAsh=1,nAsh(is)
                         iAA=iAsh+nIsh(is)
                         ipQ=ipMat(ips,is)+(iAA-1)*nOrb(ips)
                         Do jAsh=1,nAsh(js)
                            jAA=jAsh+nIsh(js)
                            ipM=(jAA-1)*nOrb(ipS)+1
*
                            rd=G2(iAsh+na(is),jAsh+na(js),
     &                            kAsh+na(ks),lash+na(ls))
                            call daxpy_(nOrb(ipS),rd,Temp(ipM),1,
     &                                 Q(ipQ),1)
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
