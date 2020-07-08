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
      Subroutine Kap_CI(iph1,iph2,ips1)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
      Real*8, Allocatable :: R(:,:)
      Call CISigma_sa(0,state_sym,state_sym,iph1,iph2,
     &                    idum,ipCI,ips1,'N')
      Call DSCAL_(nroots*ncsf(STATE_SYM),2.0d0,
     &           Work(ipin(ips1)),1)
      Call mma_allocate(R,[0,nroots-1],[0,nroots-1],label='R')
      Do i=0,nroots-1
       Do j=0,nroots-1
        R(i,j)=ddot_(nconf1,Work(ipin(ips1)+nconf1*i),1,
     &                     Work(ipin(ipci)+nconf1*j),1)
       End Do
      End Do
      Do i=0,nroots-1
       Do j=0,nroots-1
       call daxpy_(nconf1,-R(i,j),
     &                   Work(ipin(ipci)+i*nconf1),1,
     *                   Work(ipin(ipS1)+j*nconf1),1)
       End Do
      End Do
      Call mma_deallocate(R)
      Return
      End
