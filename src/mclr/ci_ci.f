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
      Subroutine Ci_Ci(ipcid,ips2)
      use ipPage, only: W
      use Arrays, only: FIMO, INT2
      Implicit Real*8(A-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
      Real*8 rDum(1)

      Call CISigma_sa(0,state_sym,state_sym,FIMO,SIZE(FIMO),
     &                Int2,SIZE(Int2),rDum,1,ipCId,ips2,.True.)
      irc=ipin(ipCId)
      irc=ipin(ipS2)
      Do i=0,nroots-1
         EC=(rin_ene+potnuc-ERASSCF(i+1))*Weight(i+1)
         Call Daxpy_(ncsf(State_Sym),EC,
     &               W(ipCId)%Vec(1+i*ncsf(state_sym)),1,
     &               W(ipS2)%Vec(1+i*ncsf(state_sym)),1)
      End Do
      Call DSCAL_(nroots*ncsf(state_SYM),2.0d0,W(ipS2)%Vec,1)
      Return
      End
