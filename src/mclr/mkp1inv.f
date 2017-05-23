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
      Subroutine mkp1inv(rdia)
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "negpre.fh"
#include "WrkSpc.fh"
#include "incdia.fh"
      Real*8 rdia(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      idisk=0
      Call Getmem('TMP1','ALLO','REAL',iptmp1,nconf)
      Call Getmem('TMP2','ALLO','REAL',iptmp2,nconf)
      Do i=1,lroots
       jdisk=idisk
       Call dDaFile(LuCIV,2,Work(ipTmp1),nconf,iDisk)
       Call ExpHinvv(rdia,Work(ipTmp1),Work(ipTmp1),0.0d0,1.0d0)
       Do j=i,lroots
         Call dDafile(luciv,2,Work(ipTmp2),nconf,jDisk)
         p1INV(itri(i,j))=DDOT_(nconf,Work(ipTmp2),1,
     &                          Work(ipTmp1),1)
       End Do
      End Do
      Call Getmem('TMP1','FREE','REAL',iptmp1,nconf)
      Call Getmem('TMP2','FREE','REAL',iptmp2,nconf)
      Return
      end
