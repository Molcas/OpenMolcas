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
      SubRoutine DMinvCI(ipSigma,rout,rC_HE_C,idsym)
      use Exp, only: NewPre
      use ipPage, only: W
      use negpre
      Implicit Real*8(a-h,o-z)

#include "real.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "incdia.fh"

      Real*8 rout(*)
      integer opout
*
*                                    -1           -1
*                               (H -E) |0><0|(H -E)|Sigma>
*                  -1             0            0
*     |rNew>=(H - E) |Sigma> - ----------------------------
*              0                               -1
*                                      <0|(H -E) |0>
*                                           0
*
      If (nconf1.gt.1) Then

         irc=ipin(ipdia)
         irc=ipin(ipSigma)
         Call exphinvv(W(ipdia)%Vec,W(ipsigma)%Vec,rout,Zero,One)
         irc=ipout(ipsigma)
         irc=opout(ipdia)
*
*        OBS <0|(H-E)|Sigma>=0 if idsym=/=1

         If (NewPre.and.idsym.eq.1) Then
*                    -1
*           rcoeff=<0|(H -E) |Sigma>
*                       0
*                 -------------------
*                            -1
*                    <0|(H -E) |0>
*                         0
*
            If (.not.ngp) Then
               irc=ipin(ipCI)
               rcoeff=ddot_(nconf1,rout,1,W(ipCI)%Vec,1)/rC_HE_C
*
*                                     -1
*              rout=rout-rocoeff*(H -E) |0>
*                                  0
               irc=ipin(ipdia)
               Call exphinvv(W(ipdia)%Vec,W(ipci)%Vec,rOUT,One,-rcoeff)
               irc=opout(ipCI)
            Else
               Call NEGP(ipdia,ipSigma,rout)
            End If

         End If

         Call DSCAL_(nconf1,Half,rout,1)

      Else

         irc=ipin(ipsigma)
         rout(1:nConf1)=W(ipSigma)%Vec(:)

      End if

      return
      end
