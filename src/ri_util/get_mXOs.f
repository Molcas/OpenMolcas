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
      Subroutine Get_mXOs(kOrb,XO,locc,nSkal,nIrrep,nOcc)
      Implicit Real*8 (a-h,o-z)
      Integer kOrb, nOcc(nIrrep), nSkal
      Real*8 XO(locc,nSkal,nIrrep)
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "exterm.fh"
#include "WrkSpc.fh"
************************************************************************
*                                                                      *
*     Statement function
*
      NBASSH(I,J)=IWORK(ip_NBASSH-1+NSYM*(J-1)+I)
*                                                                      *
************************************************************************
*                                                                      *
      Call FZero(XO,locc*nSkal*nIrrep)
*
*
*     Loop over irreps
*
      Do ir=1,nIrrep
*
*        Pointer to the next block of X_i,mu
*
         jpCMO = ip_CMOi(kOrb) - 1 + iOff_CMOi(ir,kOrb)
*        Call RecPrt('X_i,mu',' ',Work(jpCMO+1),nOcc(ir),nBas(ir))
*
*        Loop over all valence shells
*
         iOff=0
         Do isk=1,nSkal
*
*           Loop over all basis functions of this shell in this
*           irrep.
*
*           Write (*,*) 'isk,nBasSh(ir,isk)=',isk,nBasSh(ir,isk)
*
            Do ib=1,nBasSh(ir,isk)
               kb=iOff+ib ! relative SO index in this irrepp
               js=nOcc(ir)*(kb-1) ! pointer to block
*              Write (*,*) 'js=',js
*
*              Loop over all the occupied MOs and pick up the largest
*              coefficient for shell isk
*
               Do iok=1,nOcc(ir)
                  jok=jpCMO+js+iok
*                 Write (*,*) 'jok=',js+iok
*                 Write (*,*) 'Work(jok)=',Work(jok)
                  XO(iok,isk,ir)=Max(XO(iok,isk,ir),abs(Work(jok)))
               End Do
            End Do
            iOff=iOff+nBasSh(ir,isk)
         End Do
*        Call RecPrt('XO(*,*,ir)',' ',XO(1,1,ir),locc,nskal)
      End Do
*
      Return
      End
