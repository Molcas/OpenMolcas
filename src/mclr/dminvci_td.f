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
      SubRoutine DMinvCI_td(rin,rout,rome,idsym)
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
#include "negpre.fh"
#include "incdia.fh"

      Real*8 rout(*),rin(*)
*
*                                    -1           -1
*                               (H -E) |0><0|(H -E)|Sigma>
*                  -1             0            0
*     |rNew>=(H - E) |Sigma> - ----------------------------
*              0                               -1
*                                      <0|(H -E) |0>
*                                           0
*
      if (nconf1.gt.1) then

        Do i=1,nconf1
          rout(i)=rin(i)/(Work(ipin(ipdia)+i-1)+rome)
        End Do

*
* To asure orthogonal response if response is in same symmetry as wavefunction
*
        If (idsym.eq.1) Then
          r1=ddot_(nconf1,Work(ipin(ipci)),1,rout,1)

          r2=0.0d0
          Do i=0,nconf1-1
            r2=r2+Work(ipIn(ipci)+i)**2/(Work(ipin(ipDia)+i)+rome)
          End Do

          Do i=1,nconf1
            rout(i)=rout(i)-r1/r2*Work(ipin(ipci)+i-1)/
     &               (Work(ipin(ipDia)+i-1)+rome)
          end do
        end if

      else
        call dcopy_(nconf1,rin,1,rout,1)
      end if

      Call DSCAL_(nconf1,0.5d0,rout,1)
      return
      end
