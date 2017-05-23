************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Rd1Int_fciqmc(ipOvlp,ipHOne,ipKine)
************************************************************************
*                                                                      *
*     Objective: Read the header of the one-electron integral file     *
*                Extract also symmetry and basis set information.      *
*                In addition read the overlap, the nuclear attraction  *
*                and kinetic integrals.                                *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)
c      Include 'files_motra.inc'
#include "fciqmc_global.fh"
#include "WrkSpc.fh"
*
      Character*8 OneLbl
      Logical Found
*
      Call qEnter('Rd1Int')
*----------------------------------------------------------------------*
*     Read one-electron integral file header etc.                      *
*----------------------------------------------------------------------*
      Call Get_cArray('Seward Title',Header,144)
*----------------------------------------------------------------------*
*     Read no.of symm. species                                         *
*----------------------------------------------------------------------*
      Call Get_iScalar('nSym',nSym)
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
      Call Get_iArray('Symmetry operations',iOper,nSym)
*----------------------------------------------------------------------*
*     Read no. of basis functions per symm. species                    *
*----------------------------------------------------------------------*
      Call Get_iArray('nBas',nBas,nSym)
*----------------------------------------------------------------------*
*     Read no. of basis functions per symm. species                    *
*----------------------------------------------------------------------*
      nDim=0
      Do iSym=1,nSym
         nDim=nDim+nBas(iSym)
      End Do
      Call Get_cArray('Unique Basis Names',BsLbl,(LENIN4)*nDim)
*----------------------------------------------------------------------*
*     Read no. of unique atoms in the system                           *
*----------------------------------------------------------------------*
      Call Get_iScalar('Unique atoms',nAtoms)
*----------------------------------------------------------------------*
*     Read atom labels                                                 *
*----------------------------------------------------------------------*
c      Call Get_cArray('Unique Atom Names',AtLbl,(LENIN)*nAtoms)
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      Call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
*----------------------------------------------------------------------*
*     Read nuclear repulsion energy                                    *
*----------------------------------------------------------------------*
*     Call Get_PotNuc(PotNuc)
      Call Get_dScalar('PotNuc',PotNuc)
*----------------------------------------------------------------------*
*     Allocate memory for one-electron integrals                       *
*----------------------------------------------------------------------*
      nTot=0
      nTot1=0
      nTot2=0
      n2max=0
      Do iSym=1,nSym
        iBas=nBas(iSym)
        nTot=nTot+iBas
        nTot1=nTot1+iBAs*(iBas+1)/2
        nTot2=nTot2+iBas*iBas
        n2max=Max(n2max,iBas*iBas)
      End Do
*
      Call GetMem('Ovlp','Allo','Real',ipOvlp,nTot1+4)
      Call GetMem('Kine','Allo','Real',ipKine,nTot1+4)
      Call GetMem('HOne','Allo','Real',ipHOne,nTot1+4)
*----------------------------------------------------------------------*
*     Read overlap integrals                                           *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      OneLbl='Mltpl  0'
      Call RdOne(iRc,iOpt,OneLbl,iComp,Work(ipOvlp),iSyLbl)

      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read core Hamiltonian                                            *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      OneLbl='OneHam  '
      Call RdOne(iRc,iOpt,OneLbl,iComp,Work(ipHone),iSyLbl)

      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     Read kinetic energy integrals                                    *
*----------------------------------------------------------------------*
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      OneLbl='Kinetic '
      Call RdOne(iRc,iOpt,OneLbl,iComp,Work(ipKine),iSyLbl)

      If ( iRc.ne.0 ) Goto 991
*----------------------------------------------------------------------*
*     If this is a perturbative reaction field calculation then        *
*     modifiy the one-electron Hamiltonian by the reaction field and   *
*     the nuclear attraction by the cavity self-energy                 *
*----------------------------------------------------------------------*
      If ( iRFpert.ne.0 ) then
         nTemp=0
         Do iSym=1,nSym
            nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
         End Do
         Call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call Get_dScalar('RF Self Energy',ERFself)
         Call Get_dArray('Reaction field',Work(lTemp),nTemp)
         If (Found) Call NameRun('RUNFILE')
         PotNuc=PotNuc+ERFself
         Call Daxpy_(nTemp,1.0D0,Work(lTemp),1,Work(ipHone),1)
         Call GetMem('RFFLD','Free','Real',lTemp,nTemp)
      End If
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Call qExit('Rd1Int')
      Return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
991   Write (6,*) 'Rd1Int: Error reading from ONEINT'
      Write (6,*) 'OneLbl=',OneLbl
      Call QTrace()
      Call Abend()
      End
