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
      use ipPage, only: W
      use negpre
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "real.fh"
#include "Pointers.fh"
#include "incdia.fh"

      Real*8 rout(*),rin(*)
*                                    -1           -1
*                               (H -E) |0><0|(H -E)|Sigma>
*                  -1             0            0
*     |rNew>=(H - E) |Sigma> - ----------------------------
*              0                               -1
*                                      <0|(H -E) |0>
*                                           0
*
      if (nconf1.gt.1) then
         irc=ipin(ipdia)
         Do i=1,nconf1
            rout(i)=rin(i)/(W(ipdia)%Vec(i)+rome)
         End Do
*
*        To asure orthogonal response if response is in same symmetry as
*        wavefunction
*
         If (idsym.eq.1) Then
            irc=ipin(ipCI)
            r1=ddot_(nconf1,W(ipCI)%Vec,1,rout,1)

            r2=0.0d0
            irc=ipin(ipDia)
            Do i=1,nconf1
               r2=r2+W(ipCI)%Vec(i)**2/(W(ipDia)%Vec(i)+rome)
            End Do

            Do i=1,nconf1
               rout(i)=rout(i)-r1/r2*W(ipCI)%Vec(i)/
     &               (W(ipDia)%Vec(i)+rome)
            end do
         end if

      else

        rout(1:nConf1) = rin(1:nConf1)

      end if

      rout(1:nConf1) = Half * rout(1:nConf1)

      return
      end
