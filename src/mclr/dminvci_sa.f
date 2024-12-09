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
      SubRoutine DMinvCI_sa(ipSigma,rout,S)
      use ipPage, only: W
      use negpre
      use MCLR_Data, only: nConf1, ipCI
      use MCLR_Data, only: ipDia
      use input_mclr, only: nRoots,ERASSCF,nCSF,State_Sym
      Implicit None
      Integer ipSigma
      Real*8 rout(*), S(nroots,nroots,nroots)
#include "rasdim.fh"
      Real*8 rcoeff(mxroot),alpha(mxRoot)
      Integer, External:: ipIn
      Real*8, External:: DDot_
      Integer irc, i, j, k, iR, jR
      Real*8 E
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
       irc=ipin(ipsigma)
       k=0
       Do i=1,nroots
          E=ERASSCF(i)
          Do j=1,ncsf(state_SYM)
           k=k+1
           rout(k)=W(ipSigma)%Vec(k)/(W(ipdia)%Vec(j)-E)
          End Do
       End Do
       Do iR=1,nroots

!       We=weight(iR)
        E=ERASSCF(iR)
        irc=ipin(ipCI)
        Do jR=1,nroots
         rcoeff(jR)=ddot_(nconf1,rout(1+(iR-1)*ncsf(State_Sym)),1,
     &             W(ipCI)%Vec(1+(jR-1)*ncsf(State_Sym)),1)
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
     &    W(ipCI)%Vec(j+(i-1)*ncsf(State_Sym))*alpha(i)
     &    /(W(ipdia)%Vec(j)-E)
         End Do
        End Do

       End Do
      Else
       call dcopy_(nconf1*nroots,[0.0d0],0,rout,1)
      End If
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      end SubRoutine DMinvCI_sa
