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
      Subroutine Get_Molecule(ipCM,ipCoor,ipGrd,AtomLbl,nsAtom,mxdc)
      Implicit Real*8 (a-h,o-z)
#include "sbs.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "weighting.fh"
#include "Molcas.fh"
      Character*(LENIN) AtomLbl(mxdc)
      Logical TransVar, RotVar, Found
      Integer Columbus
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('Get_Molecule')
*                                                                      *
************************************************************************
*                                                                      *
*     Read initial data
*
      Call Get_iScalar('Unique atoms',nsAtom)
      If (nsAtom.gt.mxdc) Then
         Call WarningMessage(2,'Init: nsAtom.gt.mxdc')
         Write (6,*) 'nsAtom,mxdc=',nsAtom,mxdc
         Call Abend()
      End If
*
      Call Allocate_Work(ipCoor,3*nsAtom)
      Call Get_dArray('Unique Coordinates',Work(ipCoor),3*nsAtom)
*
      Call Allocate_Work(ipCM  ,nsAtom)
      Call Get_dArray('Nuclear charge',Work(ipCM),nsAtom)
*
      Call Get_iScalar('Grad ready',iGO)
      iJustGrad = iAnd(iGO, 2**0)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate gradient (it will be read later)
*     (This should eventually be removed, as ipGrd is unused...)
*
      Call Get_iScalar('Columbus',columbus)
      If ((iJustGrad.eq.1).and.(columbus.eq.1)) Then
*
*        C&M mode
*
         Call Get_iScalar('ColGradMode',iMode)
         If (iMode.eq.0) Then
            Call Get_Grad(ipGrd,Length)
         Else If (iMode.le.3) Then
            Call qpg_dArray('Grad State1',Found,Length)
            If (.not.Found .or. Length.eq.0) Then
               Call SysAbendmsg('Get_Molecule','Did not find:',
     &                          'Grad State1')
            End If
            Call GetMem('Grad','Allo','Real',ipGrd,Length)
            Call Get_dArray('Grad State1',Work(ipGrd),Length)
*
         End If
         If ( length.ne.3*nsAtom ) Then
            Call WarningMessage(2,'Init: length.ne.3*nsAtom')
            Write (6,*) 'Grad'
            Write (6,*) 'length,nsAtom=',length,nsAtom
            Call Abend()
         End If
         iJustGrad = 0
         iGO = iOr(iGO,iJustGrad)
         Call Put_iScalar('Grad ready',iGO)
      Else
*
*        M mode
*
         Call GetMem('GRAD','Allo','Real',ipGrd,3*nsAtom)
         Call FZero(Work(ipGrd),3*nsAtom)
      End If

      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*     Check if method is translational or rotational variant.
*
      TransVar=iAnd(iSBS,2**7).eq.2**7
      RotVar  =iAnd(iSBS,2**8).eq.2**8
*
      iPL=iPrintLevel(-1)
      If ((TransVar.or.RotVar).and.(iPL.gt.0)) Then
         Write (6,*)
         If (TransVar)
     &      Write (6,*) '    Gradient is translational variant!'
         If (RotVar)
     &      Write (6,*) '    Gradient is rotational variant!'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Read weights
*
      Call Qpg_dArray('Weights',Found,nData)
      If (Found.And.(nData.ge.nsAtom)) Then
*        The weights array length is actually the total number of atoms,
*        not just symmetry-unique, but the symmetry-unique ones are first
         Call GetMem('Weights','Allo','Real',ipWeights,nData)
         Call Get_dArray('Weights',Work(ipWeights),nData)
      Else
         Call SysAbendMsg('Get_Molecule',
     &        'No or wrong weights were found in the RUNFILE.','')
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Get_Molecule')
      Return
      End
