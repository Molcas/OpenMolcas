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
      Subroutine negp(ipdia,ipsigma,rout)
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "negpre.fh"
      integer opout
      real*8 rout(*)
*
      idisk=0
      irc=opout(ipdia)
      Call Getmem('Tmp','ALLO','REAL',ipTmp,nconf)
      Call Getmem('Tmp2','ALLO','REAL',ipTmp2,2*lroots)
      Call Getmem('Tmp3','ALLO','REAL',ipTmp3,2*lroots)
      Do i=1,lroots
       Call dDAFILE(luciv,2,Work(ipTmp),nconf,idisk)
       Work(ipTmp2+2*i-2)=DDOT_(nconf,rout,1,Work(ipTmp),1)
       Work(ipTmp2+2*i-1)=DDOT_(nconf,Work(ipin(ipSIgma)),1,
     &                         Work(ipTmp),1)
      End Do
      irc=ipout(ipsigma)
      Call dGeMV_('N',2*lroots,2*lroots,1.0d0,
     &                     Work(ipSS),2*lroots,Work(ipTmp2),1,
     &                     0.0d0,Work(ipTmp3),1)

      idisk=0
      Do i=1,lroots
       Call dDAFILE(luciv,2,Work(ipTmp),nconf,idisk)
       Call Exphinvv(Work(ipin(ipdia)),Work(ipTmp),rout,
     &               1.0d0,Work(iptmp3+2*i-2))
       call daxpy_(nConf,Work(iptmp3+2*i-1),Work(ipTmp),1,rout,1)
      End Do
      Call Getmem('Tmp', 'FREE','REAL',ipTmp,nconf)
      Call Getmem('Tmp2','FREE','REAL',ipTmp2,2*lroots)
      Call Getmem('Tmp3','FREE','REAL',ipTmp3,2*lroots)
*
      Return
      End
