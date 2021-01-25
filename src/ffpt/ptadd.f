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
      Subroutine PtAdd(H0,Ovlp,RR,nSize,Temp,nTemp)
*
************************************************************************
*                                                                      *
*     Objective: Construct the modified Hamiltonian                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*

#include "input.fh"
*
      Real*8 H0(nSize), Ovlp(nSize), RR(nSize), Temp(nTemp)
      Character*8 Label
      Logical Debug
      Data Debug /.False./
      Dimension idum(1)
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Read nuclear attraction and kinteic energy integrals.            *
*     Combine them to generate the one-electron Hamiltonian.           *
*     Finally read the overlap matrix.                                 *
*                                                                      *
*----------------------------------------------------------------------*
*
*
      iOpt1=1
      iOpt2=2
      iComp=1
      iSyLbl=nSym
      If(LCumulate) Then
         Label='OneHam  '
         Write(6,*)
         Write(6,*)'Adding perturbation cumulatively'
         Write(6,*)
      Else
         Label='OneHam 0'
      EndIf
      iRc=-1
      Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
      nInts=idum(1)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'PtAdd: Error reading ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      If (nInts+4.ne.nSize) Then
         Write (6,*) 'PtAdd: nInts+4.ne.nSize',nInts+4,nSize
         Call Abend
      End If
      iRc=-1
      Call RdOne(iRc,iOpt2,Label,iComp,H0,iSyLbl)
      If ( Debug ) Then
         Call PrDiOp('One Hamiltonian intgrl',nSym,nBas,H0)
         Write (6,*) 'PotNuc=',H0(nInts+4)
      End If
*
*----------------------------------------------------------------------*
*     Loop over all possible commands and branch to "special purpose"  *
*     subroutines to add perturbations.                                *
*----------------------------------------------------------------------*
*
      Call PtRela(H0,Ovlp,RR,nSize,Temp,nTemp)
      Call PtDipo(H0,Ovlp,RR,nSize,Temp,nTemp)
      Call PtQuad(H0,Ovlp,RR,nSize,Temp,nTemp)
      Call PtOkt0(H0,Ovlp,RR,nSize,Temp,nTemp)
      Call PtEfld(H0,Ovlp,RR,nSize,Temp,nTemp)
      Call PtEfgr(H0,Ovlp,RR,nSize,Temp,nTemp)
      Call PtGLbl(H0,Ovlp,RR,nSize,Temp,nTemp)
*
*----------------------------------------------------------------------*
*     If the user have requested a local (a la LoProp) perturbation    *
*     then make some modifications to the perturbation matrix.         *
*----------------------------------------------------------------------*
*                                                                      *
      If(ComStk(4,0,0,0))Call SelectLoc(H0,nSize)
*
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*
      If ( Debug ) Then
         Call PrDiOp('Core Hamiltonian',nSym,nBas,H0)
         Write (6,*) 'PotNuc=',H0(nInts+4)
      End If
      iRc=-1
      iOpt=0
      iComp=1
      Label='OneHam  '
      Call WrOne(iRc,iOpt,Label,iComp,H0,iSyLbl)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'PtAdd: Error writing to ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
*     Call Put_PotNuc(H0(nInts+4))
      Call Put_dScalar('PotNuc',H0(nInts+4))
*
*----------------------------------------------------------------------*
*     Normal Exit                                                      *
*----------------------------------------------------------------------*
*
      Return
      End
