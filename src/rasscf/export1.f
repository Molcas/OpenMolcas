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
* Copyright (C) 1994, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Export1(iFinal,CMO,DA,PA,DAO,Focc)
************************************************************************
*                                                                      *
*     purpose: Save all information relevant to geometry               *
*              optimizations                                           *
*                                                                      *
*     calling arguments:                                               *
*     iFinal  : Switch including routing information                   *
*     CMO     : MO coefficients in last CI                             *
*     DA      : 1-density of active orbitals                           *
*     PA      : 2-density of active orbitals                           *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1994                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_cfg
#endif
*
      Implicit Real*8 (a-h,o-z)
*...  Define global variables .........................................*
#include "rasdim.fh"
#include "rasscf.fh"
#include "gas.fh"
#include "general.fh"
#include "wadr.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
      Dimension CMO(*),DA(*),PA(*),DAO(*),Focc(*)
*...  Define local variables ..........................................*
      Character*8 RlxLbl,Method
      Logical SCF, Found
      Integer nTemp(8)
      Character(Len=16) mstate
*
#ifndef _DMRG_
      logical :: doDMRG = .false.
#endif
*----------------------------------------------------------------------*
*     Prologue                                                         *
*----------------------------------------------------------------------*
      Call qEnter('Export1')
*----------------------------------------------------------------------*
*     Save information pertinent to the gradient calculation           *
*----------------------------------------------------------------------*
*...  Add elementary information ......................................*
      SCF=.false.
      If (nac.eq.0.or.2*nac.eq.nactel) SCF=.true.
*
      If (SCF) then
         Do iS=1,nSym
            nTemp(is)=nAsh(is)+nIsh(is)
         End Do
         Call Put_iArray('nIsh',nTemp,nSym)
         Do iS=1,nSym
            nTemp(is)=0
         End Do
         Call Put_iArray('nAsh',nTemp,nSym)
      Else
         Call Put_iArray('nIsh',nIsh,nSym)
         Call Put_iArray('nAsh',nAsh,nSym)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Find the correct method label
*
      Method='CASSCF  '
      If (KSDFT.ne.'SCF') Method='CASDFT  '
*
*     Set flags for State-Average or Single root
*
      If (nRoots.ne.1) Then
*
*        iSA=-1 non-equivalent multi state SA-CASSCF
*        iSA=0  equivalent multi state SA-CASSCF
*        iSA=2  single root SA-CASSCF
*
         Method='CASSCFSA'
*
*        Check if equal weight SA-CASSCF
*
         iSA=0
         Do iR = 2, nRoots
            If (Weight(1).ne.Weight(iR)) iSA=-1
         End Do
         If (iSA.ne.0) Then
*
*           Check if SA-CASSCF is optimized for just one root.
*
            nW=0
            Do iR = 1, nRoots
               If (Weight(iR).ne.0.0D0) nW = nW + 1
            End Do
            If (nW.eq.1) iSA=2
         End If
         Call Put_iScalar('SA ready',iSA)
         If ((iSA.eq.0).or.(iSA.eq.-1)) Then
            mstate = '****************'
            Call Put_cArray('MCLR Root',mstate,16)
         End If
      End If
*
*     Check if it is a RASSCF function and not a CASSCF
      If (nHole1.ne.0 .or. nElec3.ne.0) Method(1:1)='R'
*     Check if it is a GASSCF function
      If (iDoGAS) Method(1:1)='G'
*     Check if it is a DMRGSCF function
      if(doDMRG)then
                        Method='DMRGSCF '
        if(nroots.ne.1) Method='DMRGSCFS'
      endif
*
      Call Put_cArray('Relax Method',Method,8)
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_iScalar('nSym',i)
      Call Put_iArray('nFro',nFro,i)
      Call Put_iArray('nDel',nDel,i)
*...  Add MO-coefficients .............................................*
      Call Put_CMO(CMO,NTOT2)
*...  Add one body density matrix in AO/SO basis ......................*
      Call Put_D1AO(DAO,NTOT1)
*...  Add one body density matrix in MO, active orbitals only .........*
      Call Put_D1MO(DA,NACPAR)
*...  Add two body density matrix in MO basis, active orbitals only ...*
      If ( .not.SCF ) Call Put_P2MO(PA,NACPR2)
*...  Next version of MOLCAS add the state to relax file ..............*
      Call Qpg_iScalar('Relax Original ro',Found)
      If (Found) Then
         Call Get_iScalar('Relax Original ro',irlxroot1)
         Call Get_iScalar('Relax CASSCF root',irlxroot2)
         If (irlxroot1.eq.irlxroot2) Then
            Call Put_iScalar('Relax Original ro',irlxroot)
         End If
      Else
         Call Put_iScalar('Relax Original ro',irlxroot)
      End If
      Call Put_iScalar('Relax CASSCF root',irlxroot)
*...  Remove overlaps (computed by rassi) .............................*
      Call Put_darray('State Overlaps',Work(ip_Dummy),0)
      Call Put_lscalar('Track Done',.False.)
*...  Add generalized Fock matrix .....................................*
      If ( ifinal.ge.1 ) then
         Call Put_Fock_Occ(Focc,ntot1)
         RlxLbl='Thrs    '
         tmp=Max(thrte,thrsx)
         Call put_dscalar(RlxLbl,tmp)
      End If
*----------------------------------------------------------------------*
*     Epilogue                                                         *
*----------------------------------------------------------------------*
      Call qExit('Export1')
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
