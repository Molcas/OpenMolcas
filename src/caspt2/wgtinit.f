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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      subroutine wgtinit(H)
      implicit real(8) (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      real(8) H(Nstate,Nstate)

      CALL QENTER('wgtinit')

      if (IPRGLB.GE.DEBUG) then
        write(6,*)' Entered wgtinit.'
      end if

* Initialize array of weights with all zeros
      call dcopy_(Nstate**2,[0.0D0],0,WORK(LDWGT),1)

* Main loop over all states to compute the weights
      do I=1,Nstate

* If it is an XDW-CASPT2 calculation, the weights are computed
        if (IFDW.and.zeta.ge.0.0d0) then
          Ebeta = H(I,I)
* Compute normalization factor
          do J=1,Nstate
            Ealpha = H(J,J)
            factor = 0.0D0
            do K=1,Nstate
              Egamma = H(K,K)
              factor = factor + exp(-zeta*(Ealpha - Egamma)**2)
            end do
            IJ = (I-1) + Nstate*(J-1)
* Compute weight according to XDW prescription
            WORK(LDWGT+IJ) = exp(-zeta*(Ealpha - Ebeta)**2)/factor
          end do

* If it is an XMS-CASPT2 calculation, all the weights are equal,
* i.e. they all are 1/Nstate
        else if (IFXMS.and.(.not.IFDW)) then
          call dcopy_(Nstate**2,[1.0D0/Nstate],0,WORK(LDWGT),1)

* If it is a normal MS-CASPT2 or a (X)DW-CASPT2 with zeta->infinity
* the weight vectors are the standard unit vectors e_1, e_2, ...
        else
          WORK(LDWGT + (Nstate*(I-1)) + (I-1)) = 1.0d0
        end if

* End of loop over states
      end do

* In case it is a XDW calculation, print out the weights
      if (IFDW.and.(IPRGLB.ge.VERBOSE)) then
        if (IFEFOCK) then
          write(6,*)' Weights calculated with <I|H0|I>:'
        else
          write(6,*)' Weights calculated with <I|H|I>:'
        end if
        call prettyprint(WORK(LDWGT),Nstate,Nstate)
      end if

      CALL QEXIT('wgtinit')

      return
      end
