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
      SubRoutine Drv1El_DMET(NAele,DMET_s,DMET_f,DMET_h,nBfn)
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
      Real*8 DMET_s(nBfn,nBfn),DMET_f(nBfn,nBfn), DMET_h(nBfn,nBfn)
      Real*8 Ccoor(3)
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
      Integer nComp, NAele
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
       If (.Not.Prprt.and..Not.Primitive_Pass.and.Do_FckInt) Then
       Call Allocate_Auxiliary()
       Call dcopy_(3,Zero,0,CoorO,1)
       OperI(1) = 1
       OperC(1) = iChBas(1)
*
       Label='FckInt  '
       Call Drv_Fck_DMET(Label,ipList,OperI,nComp,
     &                   CoorO,nOrdOp,Zero,rHrmt,OperC,
     &                   DMET_f,nBfn)
*
       Call Deallocate_Auxiliary()
       End If
       Call Gen_RelPointers(Info-1)

************************************************************************
************************************************************************
*
*   Multipole Moments starting with the overlap. If SEWARD is run in
*         the property mode we will skip the overlap integrals.
*
************************************************************************
************************************************************************
       PLabel=' '
       rHrmt=One
       iLow = 0
       mMltpl=nMltpl
       If (Prprt) iLow = 1
*     If not Douglas-Kroll and primitive pass do no property integrals
       If (Primitive_Pass) Then
           iLow=0
           If (.Not. DKroll) mMltpl=-1
       End If

       Do 10 iMltpl = iLow, nMltpl
          Write (Label,'(A,I2)') 'Mltpl ', iMltpl
          nComp = (iMltpl+1)*(iMltpl+2)/2
          Call DCopy_(3,Coor_MPM(1,iMltpl+1),1,Ccoor,1)
          Call Allocate_Auxiliary()
          iComp=0
          Do 11 ix = iMltpl, 0, -1
              If (Mod(ix,2).eq.0) Then
                  iSymX=1
              Else
                  ixyz=1
                  iSymX=2**IrrFnc(ixyz)
                  If (Ccoor(1).ne.Zero) iSymX = iOr(iSymX,1)
              End If
              Do 12 iy = iMltpl-ix, 0, -1
                 If (Mod(iy,2).eq.0) Then
                     iSymY=1
                 Else
                     ixyz=2
                     iSymY=2**IrrFnc(ixyz)
                     If (Ccoor(2).ne.Zero) iSymY = iOr(iSymY,1)
                 End If
                 iz = iMltpl-ix-iy
                 If (Mod(iz,2).eq.0) Then
                     iSymZ=1
                 Else
                     ixyz=4
                     iSymZ=2**IrrFnc(ixyz)
                     If (Ccoor(3).ne.Zero) iSymZ = iOr(iSymZ,1)
                 End If
                 iChO = Mod(ix,2)*iChBas(2)
     &              + Mod(iy,2)*iChBas(3)
     &              + Mod(iz,2)*iChBas(4)
                 OperI(1+iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ,
     &                            nIrrep),nIrrep)
                 OperC(1+iComp) = iChO
                 Call DCopy_(3,Coor_MPM(1,iMltpl+1),1,
     &                    CoorO(1+iComp*3),1)
                 iComp = iComp + 1
 12         Continue
 11      Continue

*         If(iMltpl.eq.0) Then
             Nuc = NAele
             Call Put_dScalar('Total Nuclear Charge',Nuc)
*         End If
*        write(6,*) "MltInt", MltInt
*        write(6,*) "MltMem", Mltmem
*        write(6,*) "label", label
*        write(6,*) "iplist", iplist
*        write(6,*) "OperI", OperI
*        write(6,*) "nComp",nComp
*        write(6,*) "CoorO",CoorO
*        write(6,*) "nOrdOp",nOrdOp
*        write(6,*) "Nuc",  Nuc
*        write(6,*) "OperC",OperC


         nOrdOp=iMltpl
         Call OneEl_DMET_s(MltInt,MltMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Zero,rHrmt,OperC,
     &              DMET_s,nBfn)

         Call Deallocate_Auxiliary()

 10   Continue

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
