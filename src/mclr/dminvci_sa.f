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
      SubRoutine DMinvCI_sa(ipSigma,rout,rC_HE_C,idsym,S)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
#include "negpre.fh"
#include "incdia.fh"

      Real*8 rout(*),rcoeff(mxroot),alpha(mxRoot),
     &       S(nroots,nroots,nroots)
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
       ipD=ipin(ipdia)
       ip=ipin(ipsigma)
       k=1
       Do i=1,nroots
          E=ERASSCF(i)
          Do j=0,ncsf(state_SYM)-1
           rout(k)=Work(ip+k-1)/(Work(ipD+j)-E)
           k=k+1
          End Do
       End Do
       Do iR=1,nroots

        W=weight(iR)
        E=ERASSCF(iR)
        Do jR=1,nroots
         rcoeff(jR)=ddot_(nconf1,rout(1+(iR-1)*ncsf(State_Sym)),1,
     &             Work(ipin(ipCI)+(jR-1)*ncsf(State_Sym)),1)
        End Do

        Do i=1,nroots
         alpha(i)=0.0d0
         Do j=1,nroots
          alpha(i)=alpha(i)+S(i,j,iR)*rcoeff(j)
         End Do
        End Do

        Do i=1,nroots
         Do j=1,ncsf(State_Sym)
          rout(j+(iR-1)*ncsf(State_Sym))=
     &    rout(j+(iR-1)*ncsf(State_Sym))-
     &    Work(ipin(ipCI)+j-1+(i-1)*ncsf(State_Sym))*alpha(i)
     &    /(Work(ipD+j-1)-E)
         End Do
        End Do

       End Do
      Else
       call dcopy_(nconf1*nroots,[0.0d0],0,rout,1)
      End If
      return
c Avoid unused argument warnings
      If (.False.) Then
       Call Unused_real(rC_HE_C)
       Call Unused_integer(idsym)
      End If
      end
