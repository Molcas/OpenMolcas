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
      Subroutine inter1_faiemp(Label,iBas_Lab,Coor,ZNUC,N_Cent,ipInf)
************************************************************************
*                                                                      *
* Object: A routine similar to Inter1 but the list of atoms it builds  *
*         includes the expanded fragment's atoms.                      *
*         The original Inter1 routine is not modified as it is called  *
*         from numerous places and this routine is only used for the   *
*         generation of the MOLDEN file during SCF.                    *
*                                                                      *
************************************************************************
      Implicit None
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
      Real*8 A(3),Coor(3,*),ZNUC(*)
      integer Ibas_Lab(*),N_Cent,ipInf
      Character*(LENIN) Lbl
      Character*(LENIN) Label(*)
      Logical DSCF
      Integer nDiff,mdc,ndc,iCnttp,ixyz,iCnt,kop,iCo
      Real*8  A1,A2,A3
      Integer NrOpr,iPrmt
      External NrOpr,iPrmt

*
      DSCF=.False.
      nDiff=0
      Call IniSew(ipInf,DSCF,nDiff)
*
      mdc=0
      ndc=0
      Do iCnttp=1,nCnttp
         If (pChrg(iCnttp).or.AuxCnttp(iCnttp).or.
     &       FragCnttp(iCnttp)) Then
           mdc = mdc + nCntr(iCnttp)
           Go To 99
         End If
         ixyz = ipCntr(iCnttp)
         Do iCnt=1,nCntr(iCnttp)
            mdc=mdc+1
            Lbl=LblCnt(mdc)(1:LENIN)
            call dcopy_(3,Work(ixyz),1,A,1)
            Do iCo=0,nIrrep/nStab(mdc)-1
               ndc=ndc+1
               kop=iCoSet(iCo,0,mdc)
               A1=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),1))*A(1)
               A2=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),2))*A(2)
               A3=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),4))*A(3)
               Label(ndc)=Lbl(1:LENIN)
               iBas_Lab(ndc)=icnttp
               Coor(1,ndc)=A1
               Coor(2,ndc)=A2
               Coor(3,ndc)=A3
               ZNUC(ndc)=DBLE(iAtmNr(iCnttp))
            End Do
            ixyz=ixyz+3
         End Do
 99      Continue
      End Do
      n_cent=ndc
*
      Return
      End
