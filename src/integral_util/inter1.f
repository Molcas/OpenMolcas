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
      Subroutine inter1(Label,iBas_Lab,Coor,ZNUC,N_Cent)
      Use Basis_Info
      Use Center_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8(a-h,o-z)
#include "Molcas.fh"
      Real*8 A(3),Coor(3,*),ZNUC(*)
      integer Ibas_Lab(*)
      Character*(LENIN) Lbl
      Character*(LENIN) Label(*)
      Logical DSCF
*
      DSCF=.False.
      nDiff=0
      Call IniSew(DSCF,nDiff)
*
      mdc=0
      ndc=0
      Do iCnttp=1,nCnttp
         If(dbsc(iCnttp)%Aux.or.
     &      dbsc(iCnttp)%Frag.or.
     &      dbsc(iCnttp)%pChrg) Then
           mdc = mdc + dbsc(iCnttp)%nCntr
           Go To 99
         End If
         Do iCnt=1,dbsc(iCnttp)%nCntr
            mdc=mdc+1
            Lbl=dc(mdc)%LblCnt(1:LENIN)
            A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
            Do iCo=0,nIrrep/dc(mdc)%nStab-1
               ndc=ndc+1
               kop=dc(mdc)%iCoSet(iCo,0)
               Call OA(kOp,A,Coor(1:3,ndc))
               Label(ndc)=Lbl(1:LENIN)
               iBas_Lab(ndc)=iCnttp
               ZNUC(ndc)=DBLE(dbsc(iCnttp)%AtmNr)
            End Do
         End Do
 99      Continue
      End Do
      n_cent=ndc
*
      Return
      End
