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
* Copyright (C) 1991, Roland Lindh                                     *
*               1996, Per Ake Malmqvist                                *
************************************************************************
      SubRoutine Drv1El_DMET(DMET_f,DMET_h,nBfn)
************************************************************************
*                                                                      *
* Object: driver for computation of one-electron matrices.             *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              OneEl                                                   *
*              GetDens                                                 *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1991                                             *
************************************************************************
      Use GeoList
      Use MpmC
      Use PrpPnt
      Implicit Real*8 (A-H,O-Z)
      External KnEInt
      External KnEMem
*     ipList: list of pointers to the integrals of each component
*             of the operator
*     OperI: list which irreps a particular component of the operator
*            belongs to
*     OperC: list the character of each component of the operator
*     CoorO: list of origins of the operator, one for each component
      Integer, Dimension(:), Allocatable :: ipList, OperI, OperC
      Real*8, Dimension(:), Allocatable :: CoorO, Nuc
      logical lECPnp,lPAM2np
      Real*8 DMET_f(nBfn,nBfn), DMET_h(nBfn,nBfn)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "wldata.fh"
#include "property_label.fh"
#include "oneswi.fh"
#include "warnings.fh"
      Character*8 Label
      Character*512 FName
      Integer nComp
*
      iRout = 131
      iPrint = nPrint(iRout)
*     Call qEnter('Drv1El')
*
      Call StatusLine(' Seward:',' Computing 1-electron integrals')
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
#ifdef _FDE_
      if (embPot) call EmbPotRdRun
#endif

*
*     set center selector in OneSwi to all centers (default)
*
      NDDO = .FALSE.
      If (Prprt.and.DKroll) Then
         Call WarningMessage(2,
     &               'Prprt and DKroll options can not be combined!')
         Call Quit_OnUserError()
      End If
*
*     We will always compute the following one-electron integrals per
*     default.
*     1) Multipole moments up to quadrupole moments
*     2) Kinetic energy
*     3) Nuclear Attraction
*     4) ECP contributions
*     5) One-Electron Hamiltonian
*     6) Mass-Velocity
*     7) Darwin 1-electron contact term
*
      lECPnp = lECP
      lPAM2np = lPAM2
      if (DKroll.and.Primitive_Pass) then
         lECPnp = .False.
      endif
      If (Prprt) Then
         FName=SW_FileOrb
         IF (mylen(FName).eq.0) FName='INPORB'
         Call GetDens(FName(:mylen(FName)),short,iPrint)
         Call CollapseOutput(1,'   Molecular properties:')
         Write (6,'(3X,A)')    '   ---------------------'
         Write (6,*)
      End If
************************************************************************
************************************************************************
*2)                                                                    *
*                                                                      *
*     Kinetic energy, nuclear attraction and ECP/PP integrals          *
*                                                                      *
*     Mass-velocity and One-electron Darwin contact term integrals.    *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=One
      nComp=1
*
      If (.Not.Prprt) Then
         Call Allocate_Auxiliary()
         Call dcopy_(3,Zero,0,CoorO,1)
         OperI(1) = 1
         OperC(1) = iChBas(1)
*
         Label='Kinetic '
         nOrdOp = 2
         Call OneEl_DMET(KnEInt,KnEMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Zero,rHrmt,OperC,
     &              DMET_h,nBfn)
*
         write(6,*) 'One_el Integrals Done'
         Call Deallocate_Auxiliary()
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
* 21) Atomic Fock matrix
* *
*                                                                      *
************************************************************************
************************************************************************
       Call Gen_RelPointers(-(Info-1))
       PLabel=' '
       rHrmt=One
       nComp=1
       nOrdOp = 0
       write(6,*) 'before fock'
       If (.Not.Prprt.and..Not.Primitive_Pass.and.Do_FckInt) Then
           write(6,*) 'after if fock'
       Call Allocate_Auxiliary()
       Call dcopy_(3,Zero,0,CoorO,1)
           write(6,*) 'copy coord fock'
       OperI(1) = 1
       OperC(1) = iChBas(1)
*
       Label='FckInt  '
           write(6,*) 'before drv_fck'
       Call Drv_Fck_DMET(Label,ipList,OperI,nComp,
     &                   CoorO,nOrdOp,Zero,rHrmt,OperC,
     &                   DMET_f,nBfn)
*
           write(6,*) 'after drv_fck'
       Call Deallocate_Auxiliary()
       End If
       Call Gen_RelPointers(Info-1)
************************************************************************
*                                                                      *
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*     Call qExit('Drv1El')
      Return
*
      Contains
      Subroutine Allocate_Auxiliary()
      Implicit None
*
      Call mma_Allocate(ipList,nComp,label='ipList')
      Call mma_Allocate(OperI,nComp,label='OperI')
      Call mma_Allocate(OperC,nComp,label='OperC')
      Call mma_Allocate(CoorO,3*nComp,label='CoorO')
      Call mma_Allocate(Nuc,nComp,label='Nuc')
*
      Return
      End Subroutine Allocate_Auxiliary
      Subroutine Deallocate_Auxiliary()
      Implicit None
*
      Call mma_Deallocate(OperC)
      Call mma_Deallocate(OperI)
      Call mma_Deallocate(ipList)
      Call mma_Deallocate(CoorO)
      Call mma_Deallocate(Nuc)
*
      Return
      End Subroutine Deallocate_Auxiliary
*
      End Subroutine Drv1el_DMET
