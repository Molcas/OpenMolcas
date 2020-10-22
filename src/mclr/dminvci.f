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
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
#include "negpre.fh"
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
      Call exphinvv(Work(ipin(ipdia)),Work(ipin(ipsigma)),
     &               rout,0.0d0,1.0d0)
      irc=ipout(ipsigma)
      irc=opout(ipdia)
*
*      OBS <0|(H-E)|Sigma>=0 if idsym=/=1


      If (NewPre.and.idsym.eq.1) Then
*
*                    -1
*    rcoeff=<0|(H -E) |Sigma>
*                0
*          -------------------
*                     -1
*             <0|(H -E) |0>
*                  0
*
       If (.not.ngp) Then
        rcoeff=ddot_(nconf1,rout,1,Work(ipin(ipCI)),1)
     &         /rC_HE_C
*
*                            -1
*      rout=rout-rocoeff*(H -E) |0>
*                        0
*
        Call exphinvv(Work(ipin(ipdia)),Work(ipin(ipci)),
     &              rOUT,1.0d0,-rcoeff)
        irc=opout(ipci)
       Else
        Call NEGP(ipdia,ipsigma,rout)
       End If
      End If
      Call DSCAL_(nconf1,0.5d0,rout,1)
      else
       call dcopy_(nconf1,Work(ipin(ipsigma)),1,rout,1)
      end if
      return
      end
