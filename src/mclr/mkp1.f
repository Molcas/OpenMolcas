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
      Subroutine mkp1(nEX,lst,rMat,rdiag)
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "negpre.fh"
#include "WrkSpc.fh"
      Real*8 rMat(*),rdiag(*)
      Integer lst(nex)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      idisk=0
      Call Getmem('TMP1','ALLO','REAL',iptmp1,nconf)
      Call Getmem('TMP2','ALLO','REAL',iptmp2,nconf)
      Do i=1,lroots
       Call dDaFile(LuCIV,2,Work(ipTmp1),nconf,iDisk)
       jdisk=0
       Do j=1,i
         Call dDafile(luciv,2,Work(ipTmp2),nconf,jDisk)
         rTmp=0.0d0
         Do k=1,nex
          do l=1,nex
           kk=lst(k)-1
           ll=lst(l)-1
           rtmp=rtmp+Work(ipTmp1+kk)*Work(ipTmp2+ll)*rmat(itri(k,l))
          End Do
         End Do
         Do k=1,nconf
          rtmp=rtmp+Work(ipTmp1+k-1)*Work(ipTmp2+k-1)*
     &               rdiag(k)
         End Do
         If (i.eq.j) rtmp=rtmp-ERASSCF(1)
         Do  k=1,nEx
          kk=lst(k)-1
          rtmp=rtmp-Work(ipTmp1+kk)*Work(ipTmp2+kk)*
     &              (rdiag(kk+1)-ERASSCF(1))
         End Do
         P1(itri(i,j))=rtmp
        End Do
       End Do
       Call Getmem('TMP1','FREE','REAL',iptmp1,nconf)
       Call Getmem('TMP2','FREE','REAL',iptmp2,nconf)
       Return
       end
