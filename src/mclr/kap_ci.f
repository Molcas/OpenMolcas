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
      Subroutine Kap_CI(h1,nh1,h2,nh2,ipS1)
      use ipPage, only: W
      Implicit Real*8(a-h,o-z)

#include "real.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
      Real*8, Allocatable :: R(:,:)
      Real*8 h1(nh1), h2(nh2)
      Real*8 rDum(1)
*                                                                      *
************************************************************************
*                                                                      *
       Interface
       SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,
     &                       Int2a,nInt2a,ipCI1,ipCI2, Have_2_el)
       Integer iispin, iCsym, iSSym
       Integer nInt1, nInt2s, nInt2a
       Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
       Integer ipCI1, ipCI2
       Logical Have_2_el
       End SubRoutine CISigma_sa
       End Interface
*                                                                      *
************************************************************************
*                                                                      *
      Call CISigma_sa(0,state_sym,state_sym,h1,nh1,h2,nh2,rdum,1,ipCI,
     &                ipS1,.True.)

      irc=ipin(ipS1)
      irc=ipin(ipCI)

      Call DSCAL_(nroots*ncsf(STATE_SYM),Two,W(ipS1)%Vec,1)
      Call mma_allocate(R,[0,nroots-1],[0,nroots-1],label='R')

      Do i=0,nroots-1
       Do j=0,nroots-1
        R(i,j)=ddot_(nconf1,W(ipS1)%Vec(1+nconf1*i),1,
     &                     W(ipCI)%Vec(1+nconf1*j),1)
       End Do
      End Do

      Do i=0,nroots-1
       Do j=0,nroots-1
       call daxpy_(nconf1,-R(i,j),
     &                   W(ipCI)%Vec(1+i*nconf1),1,
     *                   W(ipS1)%Vec(1+j*nconf1),1)
       End Do
      End Do

      Call mma_deallocate(R)

      Return
      End
