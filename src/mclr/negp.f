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
#include "stdalloc.fh"
#include "real.fh"
#include "negpre.fh"
      integer opout
      Real*8 rout(*)
      Real*8, Allocatable:: Tmp(:), Tmp2(:,:), Tmp3(:,:)
*
      idisk=0
      irc=opout(ipdia)

      Call mma_allocate(Tmp,nConf,Label='Tmp')
      Call mma_allocate(Tmp2,2,lRoots,Label='Tmp2')
      Call mma_allocate(Tmp3,2,lRoots,Label='Tmp3')

      Do i=1,lroots
         Call dDAFILE(luciv,2,Tmp,nconf,idisk)
         Tmp2(1,i)=DDOT_(nconf,rout,1,Tmp,1)
         Tmp2(2,i)=DDOT_(nconf,Work(ipin(ipSigma)),1,Tmp,1)
      End Do
      irc=ipout(ipsigma)
      Call dGeMV_('N',2*lroots,2*lroots,One,
     &                     Work(ipSS),2*lroots,Tmp2,1,
     &                     Zero,Tmp3,1)

      idisk=0
      Do i=1,lroots
         Call dDAFILE(luciv,2,Tmp,nconf,idisk)
         Call Exphinvv(Work(ipin(ipdia)),Tmp,rout,One,Tmp3(1,i))
         call daxpy_(nConf,Tmp3(2,i),Tmp,1,rout,1)
      End Do

      Call mma_deallocate(Tmp)
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp3)
*
      Return
      End
