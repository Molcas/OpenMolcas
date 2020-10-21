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
      Implicit Real*8(A-h,o-z)
#include "WrkSpc.fh"

#include "Input.fh"
#include "Pointers.fh"
      Call CISigma_sa(0,state_sym,state_sym,Work(ipFimo),Work(k2int),
     &                    idum,ipCId,ips2,'N')
      Do i=0,nroots-1
         EC=(rin_ene+potnuc-ERASSCF(i+1))*Weight(i+1)
       Call Daxpy_(ncsf(State_Sym),EC,
     &            Work(ipin(ipCId)+i*ncsf(state_sym)),1,
     &            Work(ipin(ipS2)+i*ncsf(state_sym)),1)
      End Do
      Call DSCAL_(nroots*ncsf(state_SYM),2.0d0,Work(ipin(ipS2)),1)
      Return
      End
