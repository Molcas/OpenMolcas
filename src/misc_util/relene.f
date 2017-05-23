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
      Subroutine RelEne(ErelMV,ErelDC,nSym,nBas,CMO,OCC,D,OP)
************************************************************************
*                                                                      *
*     Purpose: Compute relativistic correction to energy for a given   *
*              set of natural orbitals and occupation numbers          *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
*
      Dimension nBas(*),CMO(*),OCC(*),D(*),OP(*)
*----------------------------------------------------------------------*
*     Compute 1-particle density matrix in AO basis                    *
*----------------------------------------------------------------------*
      ij=0
      iOff1=0
      iOff2=0
      Do 10 iSym=1,nSym
        nBs=nBas(iSym)
        If ( nBs.eq.0 ) Goto 10
        Do 20 iBas=1,nBs
          Do 30 jBas=1,iBas
            ij=ij+1
            D(ij)=0.0D0
            Do 40 iOrb=1,nBs
              iOff3=iOff1+(iOrb-1)*nBs
              D(ij)=D(ij)
     *             +OCC(iOff2+iOrb)*CMO(iBas+iOff3)*CMO(jBas+iOff3)
40          Continue
            If ( iBas.ne.jBas ) D(ij)=2.0d0*D(ij)
30        Continue
20      Continue
        iOff1=iOff1+nBs*nBs
        iOff2=iOff2+nBs
10    Continue
*----------------------------------------------------------------------*
*     Compute energy contributions                                     *
*----------------------------------------------------------------------*
      lOp=0
      Do iSym=1,nSym
         nBs=nBas(iSym)
         lOp=lOp+(nBs**2+nBs)/2
      End Do
      ErelMV=0.0
      iRc=-1
      iOpt=1
      iComp=1
      Call RdOne(iRc,iOpt,'MassVel ',iComp,OP,iSyLbl)
      If ( iRc.eq.0 ) then
        iRc=-1
        iOpt=6
        iComp=1
        Call RdOne(iRc,iOpt,'MassVel ',iComp,OP,iSyLbl)
        ErelMV = DDOT_(lOP,D,1,OP,1)
      End If
      ErelDC=0.0
      iRc=-1
      iOpt=1
      iComp=1
      Call RdOne(iRc,iOpt,'Darwin  ',iComp,OP,iSyLbl)
      If (iRc.eq.0) then
        iRc=-1
        iOpt=6
        iComp=1
        Call RdOne(iRc,iOpt,'Darwin  ',iComp,OP,iSyLbl)
        ErelDC = DDOT_(lOP,D,1,OP,1)
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
