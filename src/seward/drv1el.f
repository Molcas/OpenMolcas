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
      SubRoutine Drv1El(DInf,nDInf,Info)
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
      External MltInt, KnEInt, MVeInt,  VeInt,  D1Int,  NAInt,  EFInt,
     &         OAMInt, OMQInt, DMSInt, WelInt, XFdInt,  PrjInt,
     &          M1Int,  M2Int, SROInt, AMPInt,  PXPInt,  PXInt,
     &          VPInt, PPInt, CntInt, EMFInt,
     &         MltInt_GIAO,
     &         KneInt_GIAO,
     &          NAInt_GIAO,
     &          dTdmu_Int
      External MltMem, KnEMem, MVeMem,  VeMem,  D1Mem,  NAMem,  EFMem,
     &         OAMMem, OMQMem, DMSMem, WelMem, XFdMem,  PrjMem,
     &          M1Mem, M2Mem, SROMem, AMPMem, PXPmem,  PXMem,
     &          VPMem, PPMem, CntMem, EMFMem,
     &         MltMem_GIAO,
     &         KneMem_GIAO,
     &          NAMem_GIAO,
     &          dTdmu_Mem
      Real*8 DInf(nDInf)
#ifdef _FDE_
      ! Embedding
      External embPotMem, embPotKernel
      Real*8, Dimension(:), Allocatable :: Emb_Int
#endif
*     ipList: list of pointers to the integrals of each component
*             of the operator
*     OperI: list which irreps a particular component of the operator
*            belongs to
*     OperC: list the character of each component of the operator
*     CoorO: list of origins of the operator, one for each component
      Integer, Dimension(:), Allocatable :: ipList, OperI, OperC
      Real*8, Dimension(:), Allocatable :: CoorO, Nuc, KnE_Int, NA_Int,
     &                                     FragP, OneHam, PtEl,
     &                                     PtNuc, SumEl, SumNuc
      Real*8, Dimension(:,:), Allocatable :: PAMexp
      External PAM2Int, FragPint, PAM2Mem, FragPMem
      External P_Int, EPEInt,
     &         P_Mem, EPEMem
      logical lECPnp,lPAM2np
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
#ifdef _FDE_
#include "embpotdata.fh"
#endif
      Integer iSymR(0:3)
      Character*8 Label
      Character*512 FName
      Real*8 Ccoor(3)
      Integer iAtmNr2(mxdbsc), nComp
      Real*8 Charge2(mxdbsc)
      Dimension dum(1),idum(1)
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
         Call GetDens(FName(:mylen(FName)),short,iPrint)
         Call CollapseOutput(1,'   Molecular properties:')
         Write (6,'(3X,A)')    '   ---------------------'
         Write (6,*)
      End If
************************************************************************
************************************************************************
*1)                                                                    *
*                                                                      *
*     Multipole Moments starting with the overlap. If SEWARD is run in *
*     the property mode we will skip the overlap integrals.            *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=One
      iLow = 0
      mMltpl=nMltpl
      If (Prprt) iLow = 1
*
*     If not Douglas-Kroll and primitive pass do no property integrals
*
      If (Primitive_Pass) Then
         iLow=0
         If (.Not. DKroll) mMltpl=-1
      End If
*
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
*
               OperI(1+iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ,
     &                            nIrrep),nIrrep)
               OperC(1+iComp) = iChO
               Call DCopy_(3,Coor_MPM(1,iMltpl+1),1,
     &                    CoorO(1+iComp*3),1)
               iComp = iComp + 1
 12         Continue
 11      Continue
*
         Call MltNuc(CoorO,Chrg,Centr,kCentr,
     &               Nuc,iMltpl,nComp)
*--- pow hack
         If(iMltpl.eq.0) Then
            Call Put_dScalar('Total Nuclear Charge',Nuc(1))
         End If
*--- pow hack
         nOrdOp=iMltpl
         Call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*                                                                      *
************************************************************************
*                                                                      *
*        Write FMM multipole moments to disk
*
         If (.Not.Prprt.and.DoFMM) Then
            Write (Label,'(A,I2)') 'FMMInt', iMltpl
            Call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
*
*           FMM overlap distribution centres:
*           Pretend they are 1-e integrals with three (x,y,z)
*           components and write to disk in canonical order
*
            If(iMltpl.eq.0) Then
               Write (Label,'(A)') 'FMMCnX'
               Call OneEl(MltInt,MltMem,Label,ipList,OperI,
     &                 nComp,CoorO,nOrdOp+1,Nuc,rHrmt,
     &                 OperC,dum,1,dum,idum,0,0,
     &                 dum,1,0)
               Write (Label,'(A)') 'FMMCnY'
               Call OneEl(MltInt,MltMem,Label,ipList,OperI,
     &                 nComp,CoorO,nOrdOp+1,Nuc,rHrmt,
     &                 OperC,dum,1,dum,idum,0,0,
     &                 dum,1,0)
               Write (Label,'(A)') 'FMMCnZ'
               Call OneEl(MltInt,MltMem,Label,ipList,OperI,
     &                 nComp,CoorO,nOrdOp+1,Nuc,rHrmt,
     &                 OperC,dum,1,dum,idum,0,0,
     &                 dum,1,0)
            End If
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        For picture-change corrected integrals.
*
         If (DKroll.and.Primitive_Pass) Then
            Write (Label,'(A,I2)') 'pMp   ', iMltpl
            PLabel='MltInt'
            Call FZero(Nuc,nComp)
            Call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp+2,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If
*                                                                      *
************************************************************************
*                                                                      *
         If (iMltpl.eq.0) Then
*           these are overlap integrals...
*           set center selector in OneSwi to single center...
            NDDO = .TRUE.
            Write (Label,'(A,I2)') 'MltplS', iMltpl
            Call OneEl(MltInt,MltMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
*           reset center selector in OneSwi to all centers...
            NDDO = .FALSE.
         End If
*
         Call Deallocate_Auxiliary()
*
 10   Continue
************************************************************************
************************************************************************
*1a)                                                                   *
*                                                                      *
*     PAM integrals                                                    *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
************************************************************************
      If (lPAM2np.and..Not.Primitive_Pass) Then
         Call molcas_open(28,'R_vect')
         PLabel=' '
         rHrmt = One
         iPAMcount=1
        Do 348 kCnttpPAM = 1, nCnttp

           iAddr=ipPAM2xp(kCnttpPAM)
           nPAMltpl=nPAM2(kCnttpPAM)

           If(nPAMltpl.lt.0) Go To 348
           Do 347 iPAMltpl=0,nPAMltpl
              nOrdOp= iPAMltpl
              nComp =(iPAMltpl+1)*(iPAMltpl+2)/2
              iPAMPrim=Int(Work(iAddr))
              iPAMBas =Int(Work(iAddr+1))

              if(iPAMBas.eq.0.or.iPAMPrim.eq.0) go to 3471
              Call dcopy_(3,[Zero],0,Ccoor,1)
              Call Allocate_Auxiliary()
              Do iComp=0,nComp-1
                 Call dcopy_(3,Work(ipCntr(kCnttpPAM)),
     &                      1,CoorO(1+3*iComp),1)
              End Do
*
*****    Define symmetry properties of the operator:
*
         iComp=0
         Do 111 ix = iPAMltpl, 0, -1
            If (Mod(ix,2).eq.0) Then
               iSymX=1
            Else
               ixyz=1
               iSymX=2**IrrFnc(ixyz)
               If (Ccoor(1).ne.Zero) iSymX = iOr(iSymX,1)
            End If
            Do 121 iy = iPAMltpl-ix, 0, -1
               If (Mod(iy,2).eq.0) Then
                  iSymY=1
               Else
                  ixyz=2
                  iSymY=2**IrrFnc(ixyz)
                  If (Ccoor(2).ne.Zero) iSymY = iOr(iSymY,1)
               End If
               iz = iPAMltpl-ix-iy
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
*
               OperC(1+iComp) = iChO
               OperI(1+iComp) = 1
*
               iComp = iComp + 1
 121        Continue
 111     Continue
*
*****    Loop over basis finction
*
         Call mma_allocate(PAMexp,iPAMPrim,2,label='PAMexp')
         Call dcopy_(iPAMPrim,Work(iAddr+2),1,PAMexp(1,1),1)
         Do iPAMf=1,iPAMBas
            Call dcopy_(iPAMPrim,Work(iAddr+2+iPAMPrim*iPAMf),1,
     &                  PAMexp(1,2),1)
            Write (Label,'(A,I2.2,I1.1,I2.2)')
     &             'PAM', kCnttpPAM,iPAMltpl,iPAMf
*

            Call dcopy_(nComp,[Zero],0,Nuc,1)

            Call OneEl(PAM2Int,PAM2Mem, Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)

c           iPAMcount=iPAMcount+1
*
         End Do
         Call mma_deallocate(PAMexp)
         Call Deallocate_Auxiliary()
 3471    iAddr=iAddr+2+iPAMPrim*(iPAMBas+1)
 347     Continue
 348    Continue
      close(28)
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
         Call dcopy_(3,[Zero],0,CoorO,1)
         OperI(1) = 1
         OperC(1) = iChBas(1)
*
         Label='Kinetic '
         nOrdOp = 2
         Call OneEl(KnEInt,KnEMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         nOrdOp = 0
*
         Label='Attract '
         Call OneEl(NAInt,NAMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,[PotNuc],rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
#ifdef _FDE_
         ! Embedding
         if (embPot) then
          Label='Embpot '
          Call OneEl(EmbPotKernel,EmbPotMem,Label,ipList,OperI,nComp,
     &               CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &               dum,1,dum,idum,0,0,
     &               dum,1,0)
         end if
#endif

*        set center selector in OneSwi to two center NA Int...
         NDDO = .TRUE.
         Label='AttractS'
         Call OneEl(NAInt,NAMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,[PotNuc],rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*        reset center selector in OneSwi to all centers...
         NDDO = .FALSE.
         If (lECPnp.and..Not.Primitive_Pass) Then
*
            Label='PrjInt  '
            Call OneEl(PrjInt,PrjMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
            Label='M1Int   '
            Call OneEl(M1Int, M1Mem, Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
            Label='M2Int   '
            Call OneEl(M2Int, M2Mem, Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
            Label='SROInt  '
            Call OneEl(SROInt,SROMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If     ! lECPnp
         If (lPP.and..Not.Primitive_Pass) Then
            Label='PPInt   '
            Call OneEl(PPInt,PPMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If
         If (lXF.and..Not.Primitive_Pass) Then
            mOrdOp=nOrd_XF
            Label='XFdInt  '
            Call OneEl(XFdInt,XFdMem,Label,ipList,OperI,nComp,
     &                 CoorO,mOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If     ! lXF
         If (lRel.and..Not.Primitive_Pass) Then
            Label='MassVel '
            nOrdOp=4
            Call OneEl(MVeInt,MVeMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
            Label='Darwin  '
            nOrdOp=0
            Call OneEl(D1Int,D1Mem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If     ! lRel
*
         Call Deallocate_Auxiliary()
      End If
************************************************************************
************************************************************************
*8a)                                                                   *
*                                                                      *
*     Velocity integrals.                                              *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=-One
      If (Vlct.and..Not.Primitive_Pass) Then
         nOrdOp = 1
         Label='Velocity'
         nComp = 3
         Call Allocate_Auxiliary()
         Call dcopy_(3*nComp,[Zero],0,CoorO,1)
         ixyz=1
         OperI(1  ) = 2**IrrFnc(ixyz)
         OperC(1  ) = iChBas(2)
         ixyz=2
         OperI(1+1) = 2**IrrFnc(ixyz)
         OperC(1+1) = iChBas(3)
         ixyz=4
         OperI(1+2) = 2**IrrFnc(ixyz)
         OperC(1+2) = iChBas(4)
*
         Call dcopy_(3,[Zero],0,Nuc,1)
         Call OneEl(VeInt,VeMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
      End If    ! Vlct
************************************************************************
************************************************************************
*8b)                                                                   *
*                                                                      *
*     Electromagnetic field radiation integrals.                       *
*                                                                      *
*     Note that the integral is not symmetric or antisymmetric!        *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=-One ! Note used
      If (EMFR.and..Not.Primitive_Pass) Then
*
*        The second term in eq. 14 of Bernadotte et al
*
         nOrdOp = 0
         Label='EMFR0'
         nComp = 2
         Call Allocate_Auxiliary()
*        Here we put in the k-vector
         Call FZero(CoorO,3*nComp)
         Call dcopy_(3,KVector,1,CoorO,1)
*
*        The electromagnetic field operator contributes to all
*        irreducible irreps, hence OperI=255. Since the operator
*        itself is not symmetry adapted OperC is set to a dummy value.
*
         OperI(1   ) = 255
         OperI(1+1 ) = 255
         OperC(1   ) = 0 ! Dummy
         OperC(1+1 ) = 0 ! Dummy
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
*
*        The first (dominating) term in eq. 14 of Bernadotte et al
*
         nOrdOp = 1
         Label='EMFR'
         nComp = 12
         Call Allocate_Auxiliary()
*        Here we put in the k-vector
         Call FZero(CoorO,3*nComp)
         Call dcopy_(3,KVector,1,CoorO,1)
*
*        The electromagnetic field operator contributes to all
*        irreducible irreps, hence OperI=255. Since the operator
*        itself is not symmetry adapted OperC is set to a dummy value.
*
         OperI(1   ) = 255
         OperI(1+1 ) = 255
         OperI(1+2 ) = 255
         OperI(1+3 ) = 255
         OperI(1+4 ) = 255
         OperI(1+5 ) = 255
         OperI(1+6 ) = 255
         OperI(1+7 ) = 255
         OperI(1+8 ) = 255
         OperI(1+9 ) = 255
         OperI(1+10) = 255
         OperI(1+11) = 255
         OperC(1   ) = 0 ! Dummy
         OperC(1+1 ) = 0 ! Dummy
         OperC(1+2 ) = 0 ! Dummy
         OperC(1+3 ) = 0 ! Dummy
         OperC(1+4 ) = 0 ! Dummy
         OperC(1+5 ) = 0 ! Dummy
         OperC(1+6 ) = 0 ! Dummy
         OperC(1+7 ) = 0 ! Dummy
         OperC(1+8 ) = 0 ! Dummy
         OperC(1+9 ) = 0 ! Dummy
         OperC(1+10) = 0 ! Dummy
         OperC(1+11) = 0 ! Dummy
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
      End If    ! EMFR
************************************************************************
************************************************************************
*9)                                                                    *
*10)                                                                   *
*                                                                      *
*     Electric field integrals.                                        *
*     Electric field gradient integrals.                               *
*                                                                      *
************************************************************************
************************************************************************
      ixyz=1
      iSymX = 2**IrrFnc(ixyz)
      ixyz=2
      iSymY = 2**IrrFnc(ixyz)
      ixyz=4
      iSymZ = 2**IrrFnc(ixyz)
      ixyz=3
      iSymXY = 2**IrrFnc(ixyz)
      ixyz=5
      iSymXZ = 2**IrrFnc(ixyz)
      ixyz=6
      iSymYZ = 2**IrrFnc(ixyz)
      ixyz=7
      iSyXYZ = 2**IrrFnc(ixyz)
*
      PLabel=' '
      rHrmt=One
      Do nOrdOp = 0, nOrdEF
*
         nComp = (nOrdOp+1)*(nOrdOp+2)/2
*
         Call Allocate_Auxiliary()
*
      Do iEF = 1, nEF
*
*        Note that this parsing is a bit different here!
*
         Write (Label,'(A,I1,I5)') 'EF',nOrdOp,iEF
         Call dcopy_(3,Work(ipEF+(iEF-1)*3),1,Ccoor,1)
*
         iSymR(0) = 1
         If (Ccoor(1).ne.Zero) iSymR(0) = iOr(iSymR(0),iSymX)
         If (Ccoor(2).ne.Zero) iSymR(0) = iOr(iSymR(0),iSymY)
         If (Ccoor(3).ne.Zero) iSymR(0) = iOr(iSymR(0),iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymR(0) = iOr(iSymR(0),iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymR(0) = iOr(iSymR(0),iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymR(0) = iOr(iSymR(0),iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymR(0) = iOr(iSymR(0),iSyXYZ)
*
         ixyz=1
         iSym=2**IrrFnc(ixyz)
         If (Ccoor(1).ne.Zero ) iSym = iOr(iSym,1)
         iSymR(1)=iSym
*
         ixyz=2
         iSym=2**IrrFnc(ixyz)
         If (Ccoor(2).ne.Zero ) iSym = iOr(iSym,1)
         iSymR(2)=iSym
*
         ixyz=4
         iSym=2**IrrFnc(ixyz)
         If (Ccoor(3).ne.Zero ) iSym = iOr(iSym,1)
         iSymR(3)=iSym
*
         iComp=0
         Do ix = nOrdOp, 0, -1
            Do iy = nOrdOp-ix, 0, -1
               iz=nOrdOp-ix-iy
               iComp = iComp + 1
*
               iSymX=1
               If (Mod(ix,2).ne.0) iSymX=iSymR(1)
               iSymCX=MltLbl(iSymR(0),iSymX,nIrrep)
               iSymY=1
               If (Mod(iy,2).ne.0) iSymY=iSymR(2)
               iSymCXY=MltLbl(iSymCX,iSymY,nIrrep)
               iSymZ=1
               If (Mod(iz,2).ne.0) iSymZ=iSymR(3)
*
               OperI(1+(iComp-1)) = MltLbl(iSymCXY,iSymZ,nIrrep)
               OperC(1+(iComp-1)) = Mod(ix,2)*iChBas(2)
     &                              + Mod(iy,2)*iChBas(3)
     &                              + Mod(iz,2)*iChBas(4)
*
               Call dcopy_(3,Ccoor,1,CoorO(1+(iComp-1)*3),1)
            End Do
         End Do
*
         Call EFNuc(CoorO,Chrg,Centr,kCentr,
     &               Nuc,nOrdOp)
         Call OneEl(EFInt,EFMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*                                                                      *
************************************************************************
*                                                                      *
*        For picture-change corrected property integrals.
*
         If (DKroll.and.Primitive_Pass) Then
            Write (Label,'(A,I1,I5)') 'PP',nOrdOp,iEF
            PLabel='EFInt '
            Call FZero(Nuc,nComp)
            Call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp+2,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End Do
         Call Deallocate_Auxiliary()
*
* For a properties calculation, read the EF values saved in a temp file
* and write the sum through Add_Info
*
         If (PrPrt.and.nEF.gt.0) Then
           Call mma_allocate(PtEl,nComp,label='PtEl')
           Call mma_allocate(PtNuc,nComp,label='PtNuc')
           Call mma_allocate(SumEl,nComp,label='SumEl')
           Call mma_allocate(SumNuc,nComp,label='SumNuc')
           Call FZero(SumEl,nComp)
           Call FZero(SumNuc,nComp)
*          Read and sum the values
           LuTmp=10
           Call DaName(LuTmp,'TMPPRP')
           iDisk=0
           Do iEf=1,nEF
             Call dDaFile(LuTmp,2,PtEl,nComp,iDisk)
             Call dDaFile(LuTmp,2,PtNuc,nComp,iDisk)
             Call DaXpY_(nComp,One,PtEl,1,SumEl,1)
             Call DaXpY_(nComp,One,PtNuc,1,SumNuc,1)
           End Do
           Call DaClos(LuTmp)
*          set the tolerance according to the total number of centers
*          (assuming error scales with sqrt(nEF))
           iTol=5
           iTol=iTol-NInt(Half*Log10(Dble(nEF)))
           Write (label,'(a,i1,a)') 'EF',nOrdOp,'   el'
           Call Add_Info(label,SumEl,nComp,iTol)
           Write (label,'(a,i1,a)') 'EF',nOrdOp,'  nuc'
           Call Add_Info(label,SumNuc,nComp,iTol)
           Call mma_deallocate(PtEl)
           Call mma_deallocate(PtNuc)
           Call mma_deallocate(SumEl)
           Call mma_deallocate(SumNuc)
         End If
*
      End Do
************************************************************************
************************************************************************
*12)                                                                   *
*                                                                      *
*     Orbital angular momentum integrals.                              *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=-One
      If (lOAM.and..Not.Primitive_Pass) Then
         Label='AngMom  '
         nComp = 3
         nOrdOp = 2
         Call Allocate_Auxiliary()
         Call dcopy_(3,Work(ipOAM),1,CoorO(1  ),1)
         Call dcopy_(3,Work(ipOAM),1,CoorO(1+3),1)
         Call dcopy_(3,Work(ipOAM),1,CoorO(1+6),1)
         Call dcopy_(3,Work(ipOAM),1,Ccoor,1)
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
*
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iChOx = iChBas(3) + iChBas(4)
         OperI(1  ) = iSymLx
         OperC(1  ) = iChOx
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iChOy = iChBas(4) + iChBas(2)
         OperI(1+1) = iSymLy
         OperC(1+1) = iChOy
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
         iChOz = iChBas(2) + iChBas(3)
         OperI(1+2) = iSymLz
         OperC(1+2) = iChOz
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(OAMInt,OAMMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
      End If   ! lOAM
************************************************************************
************************************************************************
*12b)                                                                  *
*                                                                      *
*     Orbital Magnetic Quadrupole integrals.                           *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=-One
      If (lOMQ.and..Not.Primitive_Pass) Then
         Label='OMQ     '
         nComp = 9
         nOrdOp = 3
         Call Allocate_Auxiliary()
*
         Call dcopy_(nComp,[Work(ipOMQ  )],0,CoorO(1  ),3)
         Call dcopy_(nComp,[Work(ipOMQ+1)],0,CoorO(1+1),3)
         Call dcopy_(nComp,[Work(ipOMQ+2)],0,CoorO(1+2),3)
         Call dcopy_(3,Work(ipOMQ),1,Ccoor,1)
*
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
*
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
*
* Calculates M_ij = r_j*L_i + L_i*r_j = 2*r_j*L_i + i hbar E_ijk r_k
* Since the i hbar is included outside we could do
* M_ij = r_j*L_i + L_i*r_j = 2*r_j*L_i + E_ijk r_k
* We use all nine components even if we only need 6 of these later.
*
* Mxx
         iSymxLx = MltLbl(iSymCx,iSymLx,nIrrep)
         iChOxx = iChBas(15)
         OperI(1  ) = iSymxLx
         OperC(1  ) = iChOxx
* Mxy
         iSymxLy = iSymCz
         iChOxy = iChBas(4)
         OperI(1+1) = iSymxLy
         OperC(1+1) = iChOxy
* Mxz
         iSymxLz = iSymCy
         iChOxz = iChBas(3)
         OperI(1+2) = iSymxLz
         OperC(1+2) = iChOxz
* Myx
         iSymyLx = iSymCz
         iChOyx = iChBas(4)
         OperI(1+3) = iSymyLx
         OperC(1+3) = iChOyx
* Myy
         iSymyLy = MltLbl(iSymCy,iSymLy,nIrrep)
         iChOyy = iChBas(15)
         OperI(1+4) = iSymyLy
         OperC(1+4) = iChOyy
* Myz
         iSymyLz = iSymCx
         iChOyz = iChBas(4)
         OperI(1+5) = iSymyLz
         OperC(1+5) = iChOyz
* Mzx
         iSymzLx = iSymCy
         iChOzx = iChBas(3)
         OperI(1+6) = iSymzLx
         OperC(1+6) = iChOzx
* Mzy
         iSymzLy = iSymCx
         iChOzy = iChBas(4)
         OperI(1+7) = iSymzLy
         OperC(1+7) = iChOzy
* Mzz
         iSymzLz = MltLbl(iSymCz,iSymLz,nIrrep)
         iChOzz = iChBas(15)
         OperI(1+8) = iSymzLz
         OperC(1+8) = iChOzz
*
         Call DCopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(OMQInt,OMQMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
      End If   ! lOMQ
************************************************************************
************************************************************************
*13)                                                                   *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
************************************************************************
      If (DKroll.and.Primitive_Pass) then
         rHrmt=One
         Label='pVp     '
         PLabel='NAInt '
         nOrdOp=2
         nComp = 1
         Call Allocate_Auxiliary()
         Call dcopy_(3,[Zero],0,CoorO,1)
         OperI(1  ) = 1
         OperC(1  ) = iChBas(1)

         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)

         Call Deallocate_Auxiliary()
*
         If (BSS) Then
            rHrmt=-One
            nOrdOp = 1
            nComp = 3
            Call Allocate_Auxiliary()
*
            Call dcopy_(3*nComp,[Zero],0,CoorO,1)
            Call dcopy_(3,[Zero],0,Nuc,1)
*
            ixyz=1
            OperI(1  ) = 2**IrrFnc(ixyz)
            OperC(1  ) = iChBas(2)
            ixyz=2
            OperI(1+1) = 2**IrrFnc(ixyz)
            OperC(1+1) = iChBas(3)
            ixyz=4
            OperI(1+2) = 2**IrrFnc(ixyz)
            OperC(1+2) = iChBas(4)
*
            Label='pV      '
            PLabel='NAInt '
            Call OneEl(PXInt,PXMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
*
            Label='Vp      '
            Call OneEl(VPInt,VPMem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
*
            Call Deallocate_Auxiliary()
         End If    ! BSSInt
      End If
************************************************************************
************************************************************************
*14)                                                                   *
*                                                                      *
*     Diamagnetic shielding integrals.                                 *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=One
*
      mDMS=nDMS
      If (Primitive_Pass) mDMS=0
*
      Do 1500 iDMS = 1, mDMS
         Write (Label,'(A,I2,I2)') 'DMS ',1,iDMS
         nComp = 9
         nOrdOp=2
         Call dcopy_(3,Work(ipDMS+(iDMS-1)*3),1,Ccoor,1)
         Call Allocate_Auxiliary()
         iSymC = 1
         If (Ccoor(1).ne.Zero) iSymC = iOr(iSymC,iSymX)
         If (Ccoor(2).ne.Zero) iSymC = iOr(iSymC,iSymY)
         If (Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSyXYZ)
*
         iComp = 0
         iC = 0
         Do 1510 ix = 1, 0, -1
         Do 1510 iy = 1-ix, 0, -1
            iz=1-ix-iy
            iC = iC + 1
            iChO1 = iChBas(iC+1)
            ixyz=0
            If (Mod(ix,2).ne.0) ixyz=iOr(ixyz,1)
            If (Mod(iy,2).ne.0) ixyz=iOr(ixyz,2)
            If (Mod(iz,2).ne.0) ixyz=iOr(ixyz,4)
            iSym = 2**IrrFnc(ixyz)
            If (Ccoor(iC).ne.Zero) iSym = iOr(iSym,1)
            iD = 0
            Do 1511 jx = 1, 0, -1
            Do 1511 jy = 1-jx, 0, -1
               jz=1-jx-jy
               iD = iD + 1
               iChO2 = iChBas(iD+1)
               jxyz=0
               If (Mod(jx,2).ne.0) jxyz=iOr(jxyz,1)
               If (Mod(jy,2).ne.0) jxyz=iOr(jxyz,2)
               If (Mod(jz,2).ne.0) jxyz=iOr(jxyz,4)
               iSymD = 2**IrrFnc(jxyz)
               If (Dxyz(iD).ne.Zero) iSymD = iOr(iSymD,1)
               If (iC.eq.iD) Then
                  i2 = iD + 1
                  If (i2.gt.3) i2 = i2 - 3
                  i3 = iD + 2
                  If (i3.gt.3) i3 = i3 - 3
                  iChO = iAnd(iChBas(i2+1),iChBas(i3+1))
               Else
                  iChO = iOr(iChO1,iChO2)
               End If
               OperI(1+iComp) = MltLbl(iSymD,MltLbl(iSym,iSymC,
     &                            nIrrep),nIrrep)
               OperC(1+iComp) = iChO
               Call dcopy_(3,Ccoor,1,CoorO(1+iComp*3),1)
               iComp = iComp + 1
 1511       Continue
 1510    Continue
         Call dcopy_(3,Dxyz,1,CoorO(1+3),1)
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(DMSInt,DMSMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
1500  Continue
************************************************************************
************************************************************************
*15)                                                                   *
*                                                                      *
*     Nuclear attraction integrals for finite centers.                 *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=One
      Label='Center  '
************************************************************************
************************************************************************
*16)                                                                   *
*                                                                      *
*     Spherical well integrals.                                        *
*                                                                      *
************************************************************************
************************************************************************
      PLabel=' '
      rHrmt=One
      If (.Not.Prprt.and..Not.Primitive_Pass) Then
         nComp=1
         iWel = 0
         Call Allocate_Auxiliary()
         Call dcopy_(3,[Zero],0,CoorO,1)
         OperI(1) = 1
         OperC(1) = iChBas(1)
         Do 1600 iWel = 1, nWel
            r0   = Work(ipWel+(iWel-1)*3  )
            ExpB = Work(ipWel+(iWel-1)*3+1)
            Write (Label,'(A,I4)') 'Well',iWel
            Call OneEl(WelInt,WelMem,Label,ipList,OperI,nComp,
     &                 CoorO,iWel,[Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
 1600    Continue
         Call Deallocate_Auxiliary()
      End If  ! .Not.Prprt
************************************************************************
************************************************************************
*5)                                                                    *
*                                                                      *
*     One-electron Hamiltonian integrals.                              *
*                                                                      *
************************************************************************
************************************************************************
      If (.Not.Prprt.and..Not.Primitive_Pass) Then
         Call mma_allocate(KnE_Int,n2Tri(1)+4,label='KnE_Int')
         Call mma_allocate(NA_Int,n2Tri(1)+4,label='NA_Int')
#ifdef _FDE_
         ! Embedding
         if (embpot) Call mma_allocate(Emb_Int,n2Tri(1)+4,
     &                                 label='Emb_Int')
#endif
         iOpt = 0
         iRC = -1
         Label='Kinetic '
         Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
         If (iRC.ne.0) then
            Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
            Call Quit(_RC_IO_ERROR_READ_)
         End If
         Label='Attract '
         iRC = -1
         Call RdOne(iRC,iOpt,Label,1,NA_Int,lOper)
         If (iRC.ne.0) then
            Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
            Call Quit(_RC_IO_ERROR_READ_)
         End If
         Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
#ifdef _FDE_
         ! Embedding
         if (embpot) then
          if (embPotInBasis) then
!           write(*,*) "ENTER"
           ! If the potential is given in basis set representation it
           ! has not been calculated with a OneEl call and is just read
           ! from file here.
           iunit = isFreeUnit(1)
           call molcas_open(iunit, embPotPath)
           do iEmb=1, n2Tri(1)
            read(iunit,*) Emb_Int(iEmb)
!            write(*,*) iEmb-1, ": ", Emb_Int(iEmb)
           end do
           close(iunit)
          else
           Label='embpot  '
           iRC=-1
           Call RdOne(iRC,iOpt,Label,1,Emb_Int,lOper)
           If (iRC.ne.0) then
              Call WarningMessage(2,
     &                    'Drv1El: Error reading ONEINT;'
     &                  //'Label='//Label)
              Call Quit(_RC_IO_ERROR_READ_)
           End If
          end if
          Call DaXpY_(n2Tri(1)+4,One,Emb_Int,1,NA_Int,1)
         end if
#endif
*
*--------Add contribution from ECP
*
         If (lECPnp) Then
            Label='PrjInt  '
            iRC = -1
            Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
            If (iRC.ne.0) then
               Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
               Call Quit(_RC_IO_ERROR_READ_)
            End If
            Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
            Label='M1Int   '
            iRC = -1
            Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
            If (iRC.ne.0) then
               Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
               Call Quit(_RC_IO_ERROR_READ_)
            End If
            Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
            Label='M2Int   '
            iRC = -1
            Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
            If (iRC.ne.0) then
               Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
               Call Quit(_RC_IO_ERROR_READ_)
            End If
            Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
            Label='SROInt  '
            iRC = -1
            Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
            If (iRC.ne.0) then
               Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
               Call Quit(_RC_IO_ERROR_READ_)
            End If
            Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
         End If   ! lECPnp
*
*--------Add contributions from the Pseudo Potential
*
         If (lPP) Then
            Label='PPInt   '
            iRC = -1
            Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
            If (iRC.ne.0) then
               Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
               Call Quit(_RC_IO_ERROR_READ_)
            End If
            Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
         End If
*
*--------Add contributions from the external field
*
         If (lXF) Then
            Label='XFdInt  '
            iRC = -1
            Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
            If (iRC.ne.0) then
               Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
               Call Quit(_RC_IO_ERROR_READ_)
            End If
            Call DaXpY_(n2Tri(1)+4,One,KnE_Int,1,NA_Int,1)
         End If ! lXF
*
*--------Add contributions from Spherical wells
*
         If (nWel.ne.0) Then
            Do iWel = 1, nWel
               Fact=Work(ipWel-1+(iWel-1)*3+3)
               Write (Label,'(A,I4)') 'Well',iWel
               iRC = -1
               Call RdOne(iRC,iOpt,Label,1,KnE_Int,lOper)
               If (iRC.ne.0) then
                  Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
                  Call Quit(_RC_IO_ERROR_READ_)
               End If
               Call DaXpY_(n2Tri(1)+4,Fact,KnE_Int,1,NA_Int,1)
            End Do
         End If  ! nWel.ne.0
*
         Label='OneHam  '
         If (iPrint.ge.10) Call PrMtrx(Label,[lOper],1,[1],NA_Int)
         iRC = -1
         Call WrOne(iRC,iOpt,Label,1,NA_Int,lOper)
         If (iRC.ne.0) then
            Call WarningMessage(2,
     &                  'Drv1El: Error writing ONEINT;'
     &                //'Label='//Label)
            Call Quit(_RC_IO_ERROR_WRITE_)
         End If

#ifdef _FDE_
         ! Embedding
         if (embpot) then
          Label='embpot  '
          iRC = -1
          Call WrOne(iRC,iOpt,Label,1,Work(ipEmb),lOper)
          If (iRC.ne.0) then
             Call WarningMessage(2,
     &                   'Drv1El: Error writing ONEINT;'
     &                 //'Label='//Label)
             Call Quit(_RC_IO_ERROR_WRITE_)
          End If
         end if
#endif
*
         Label='OneHam 0'
         iRC = -1
         Call WrOne(iRC,iOpt,Label,1,NA_Int,lOper)
         If (iRC.ne.0) then
            Call WarningMessage(2,
     &                  'Drv1El: Error writing ONEINT;'
     &                //'Label='//Label)
            Call Quit(_RC_IO_ERROR_WRITE_)
         End If
         Call mma_deallocate(NA_Int)
         Call mma_deallocate(KnE_Int)

#ifdef _FDE_
         ! Embedding
         if (embPot) Call mma_deallocate(Emb_Int)
#endif
      End If
************************************************************************
************************************************************************
*17)                                                                   *
*                                                                      *
*     Angular momentum products (PAM)                                  *
*                                                                      *
************************************************************************
************************************************************************
* Hermitized products of angular momentum integrals
* Component(1) is Lx*Lx
* Component(2) is (Lx*Ly+Ly*Lx)/2, etc.
* Coded P-A Malmqvist, Garching, Nov 1996
      PLabel=' '
      rHrmt=-One
      If (lAMP.and..Not.Primitive_Pass) Then
         Label='AMProd  '
         nComp = 6
         nOrdOp = 2
         Call Allocate_Auxiliary()
         Call dcopy_(nComp,[Work(ipAMP  )],0,CoorO(1  ),3)
         Call dcopy_(nComp,[Work(ipAMP+1)],0,CoorO(1+1),3)
         Call dcopy_(nComp,[Work(ipAMP+2)],0,CoorO(1+2),3)
         Call dcopy_(3,Work(ipAMP),1,Ccoor,1)
C Symmetry labels iSymX  for operator d/dx, etc.
C Symmetry labels iSymLx for operator Lx, etc.
C Characters iChOx for operator Lx, etc.
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iChOx = iChBas(3) + iChBas(4)
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iChOy = iChBas(4) + iChBas(2)
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
         iChOz = iChBas(2) + iChBas(3)

C Symmetry labels and characters of products. Let G be the full
C  molecular point group, and Gsub=subgroup of G=stabilizer of
C gauge origin. The totally symmetric irrep of Gsub can be
C decomposed into irreps of G.
C Then symmetry label=packed array of bits, one for each irrep
C of G. The bit is set, if that irrep is included in the
C decomposition of the totally symmetric irrep of Gsub.
         OperI(1  )=MltLbl(iSymLx,iSymLx,nIrrep)
         OperI(1+1)=MltLbl(iSymLx,iSymLy,nIrrep)
         OperI(1+2)=MltLbl(iSymLx,iSymLz,nIrrep)
         OperI(1+3)=MltLbl(iSymLy,iSymLy,nIrrep)
         OperI(1+4)=MltLbl(iSymLy,iSymLz,nIrrep)
         OperI(1+5)=MltLbl(iSymLz,iSymLz,nIrrep)
         OperC(1  )=0
         OperC(1+1)=iEOr(iChOx,iChOy)
         OperC(1+2)=iEOr(iChOx,iChOz)
         OperC(1+3)=0
         OperC(1+4)=iEOr(iChOy,iChOz)
         OperC(1+5)=0
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(AMPInt,AMPMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
      End If
************************************************************************
************************************************************************
*17)                                                                   *
*                                                                      *
*     Contact term integrals                                           *
*                                                                      *
************************************************************************
************************************************************************
      ixyz=1
      iSymX = 2**IrrFnc(ixyz)
      ixyz=2
      iSymY = 2**IrrFnc(ixyz)
      ixyz=4
      iSymZ = 2**IrrFnc(ixyz)
      ixyz=3
      iSymXY = 2**IrrFnc(ixyz)
      ixyz=5
      iSymXZ = 2**IrrFnc(ixyz)
      ixyz=6
      iSymYZ = 2**IrrFnc(ixyz)
      ixyz=7
      iSyXYZ = 2**IrrFnc(ixyz)
      PLabel=' '
      rHrmt=One
*
      mCnt=0
      If (nOrdEF.eq.2) mCnt=nEF
*
      nComp = 1
      nOrdOp = 1
      Do iCnt = 1, mCnt
         Write (Label,'(A,I5)') 'Cnt',iCnt
         Call dcopy_(3,Work(ipEF+(iCnt-1)*3),1,Ccoor,1)
         Call Allocate_Auxiliary()
*
         iSymR(0) = 1
         If (Ccoor(1).ne.Zero) iSymR(0) = iOr(iSymR(0),iSymX)
         If (Ccoor(2).ne.Zero) iSymR(0) = iOr(iSymR(0),iSymY)
         If (Ccoor(3).ne.Zero) iSymR(0) = iOr(iSymR(0),iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymR(0) = iOr(iSymR(0),iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymR(0) = iOr(iSymR(0),iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymR(0) = iOr(iSymR(0),iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymR(0) = iOr(iSymR(0),iSyXYZ)
*
         OperI(1) = iSymR(0)
         OperC(1) = 0
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(CntInt,CntMem,Label,ipList,OperI,nComp,
     &              Ccoor,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*                                                                      *
************************************************************************
*                                                                      *
*        For picture-change corrected integrals.
*
         If (DKroll.and.Primitive_Pass) Then
            Write (Label,'(A,I2)') 'pCp   ', iCnt
            PLabel='CntInt'
            Call FZero(Nuc,nComp)
            Call OneEl(PXPInt,PXPMem,Label,ipList,OperI,nComp,
     &                 CCoor,nOrdOp+2,Nuc,rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
         End If
*
         Call Deallocate_Auxiliary()
      End Do
*
* For a properties calculation, read the CNT values saved in a temp file
* and write the sum through Add_Info
*
      If (PrPrt.and.mCnt.gt.0) Then
        SumEl=Zero
        SumNuc=Zero
*       Read and sum the values
        LuTmp=10
        Call DaName(LuTmp,'TMPPRP')
        iDisk=0
        Do iCnt=1,mCnt
          Call dDaFile(LuTmp,2,PtEl,1,iDisk)
          Call dDaFile(LuTmp,2,PtNuc,1,iDisk)
          SumEl=SumEl+PtEl
          SumNuc=SumNuc+PtNuc
        End Do
        Call DaClos(LuTmp)
*       set the tolerance according to the total number of centers
*       (assuming error scales with sqrt(mCnt))
        iTol=5
        iTol=iTol-NInt(Half*Log10(Dble(mCnt)))
        Write (label,'(a,a)') 'CNT','   el'
        Call Add_Info(label,SumEl,1,iTol)
        Write (label,'(a,a)') 'CNT','  nuc'
        Call Add_Info(label,SumNuc,1,iTol)
      End If
*
************************************************************************
************************************************************************
*18)                                                                   *
*     Gradient of overlap integrals with respect to the magnetic field *
*                                                                      *
************************************************************************
************************************************************************
      If (GIAO.and..Not.Primitive_Pass) Then
      PLabel=' '
      rHrmt=-One
      nB=3
      iLow = 0
      mMltpl=0 ! Do only overlap.
C     mMltpl=-1 ! Do only overlap.
      Do iMltpl = iLow, mMltpl
         Write (Label,'(A,I2)') 'dMP/dB', iMltpl
         mComp = (iMltpl+1)*(iMltpl+2)/2
         nComp = mComp * nB
         Call DCopy_(3,Coor_MpM(1,iMltpl+1),1,Ccoor,1)
         Call Allocate_Auxiliary()
*
         iComp=0
         Do ix = iMltpl, 0, -1
*
*---------- Pick up which irrep each of the cartesian components of
*           the operator belongs to. If the operator is associated
*           with a center other than origin add the total symmetric
*           irrep.
*
            If (Mod(ix,2).eq.0) Then
               iSymX=1
            Else
               ixyz=1
               iSymX=2**IrrFnc(ixyz)
               If (Ccoor(1).ne.Zero) iSymX = iOr(iSymX,1)
            End If
            Do iy = iMltpl-ix, 0, -1
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
*
*------------- Multiply cartesian components to generate which irreps
*              the current element of the multipole moment operator
*              belong to. The lowest significant bit if set indicate that
*              a particular irrep is included.
*
               iTemp = MltLbl(iSymX,MltLbl(iSymY,iSymZ,nIrrep),nIrrep)
*
               iChO = Mod(ix,2)*iChBas(2)
     &              + Mod(iy,2)*iChBas(3)
     &              + Mod(iz,2)*iChBas(4)
*
*------------- Now combine with the character of the first derivative
*              with respect to the magnetic field.
*
               iSymRx=2**IrrFnc(1)
               iSymRy=2**IrrFnc(2)
               iSymRz=2**IrrFnc(4)
*
               iB = 1
               iChOx = Mod(ix  ,2)*iChBas(2)
     &               + Mod(iy+1,2)*iChBas(3)
     &               + Mod(iz+1,2)*iChBas(4)
               OperC(1+(iB-1)*mComp+iComp) = iChOx
               iSymBx = MltLbl(iSymRy,iSymRz,nIrrep)
               OperI(1+(iB-1)*mComp+iComp) =
     &            MltLbl(iTemp,iSymBx,nIrrep)
               Call DCopy_(3,Coor_MPM(1,iMltpl+1),1,
     &                      CoorO(1+((iB-1)*mComp+iComp)*3),1)
*
               iB = 2
               iChOy = Mod(ix+1,2)*iChBas(2)
     &               + Mod(iy  ,2)*iChBas(3)
     &               + Mod(iz+1,2)*iChBas(4)
               OperC(1+(iB-1)*mComp+iComp) = iChOy
               iSymBy = MltLbl(iSymRz,iSymRx,nIrrep)
               OperI(1+(iB-1)*mComp+iComp) =
     &            MltLbl(iTemp,iSymBy,nIrrep)
               Call DCopy_(3,Coor_MPM(1,iMltpl+1),1,
     &                      CoorO(1+((iB-1)*mComp+iComp)*3),1)
*
               iB = 3
               iChOz = Mod(ix+1,2)*iChBas(2)
     &               + Mod(iy+1,2)*iChBas(3)
     &               + Mod(iz  ,2)*iChBas(4)
               OperC(1+(iB-1)*mComp+iComp) = iChOz
               iSymBz = MltLbl(iSymRx,iSymRy,nIrrep)
               OperI(1+(iB-1)*mComp+iComp) =
     &            MltLbl(iTemp,iSymBz,nIrrep)
               Call DCopy_(3,Coor_MPM(1,iMltpl+1),1,
     &                      CoorO(1+((iB-1)*mComp+iComp)*3),1)
*
               iComp = iComp + 1
            End Do
         End Do
*
*        Zero nuclear contribution.
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(MltInt_GIAO,MltMem_GIAO,
     &              Label,ipList,OperI,nComp,
     &              CoorO,iMltpl,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
*
      End Do
*
      End If
************************************************************************
************************************************************************
*19)                                                                   *
*     Atomic mean field integrals for spin-orbit calculations          *
*                                                                      *
************************************************************************
************************************************************************
      If (lAMFI.and..Not.Prprt.and..Not.Primitive_Pass) Then
         PLabel=' '
         rHrmt=-One
         Label='AMFI    '
         nComp = 3
         Call Allocate_Auxiliary()
         OperI(1  )=2**IrrFnc(6)
         OperC(1  )=iChBas(7)
         OperI(1+1)=2**IrrFnc(5)
         OperC(1+1)=iChBas(6)
         OperI(1+2)=2**IrrFnc(3)
         OperC(1+2)=iChBas(4)

* BP - Turn off AMFI integrals for certain atom types
*      as requested by the PAMF keyword
c         write(6,*) "nPAMFI:", nPAMFI
c         write(6,*) "iPAMFI:", iPAMFI(1:nPAMFI)

         Do iAtm=1,nCnttp
           iAtmNr2(iAtm) = iAtmNr(iAtm)
           Charge2(iAtm) = Charge(iAtm)

c           write(6,*) "iAtmNr2(iAtm)",iAtm, iAtmNr2(iAtm)
c           write(6,*) "Charge(iAtm)", iAtm, Charge(iAtm)

           do iPAM=1,nPAMFI
             if(iAtmNr(iAtm).EQ.iPAMFI(iPAM)) then
               write(6,*) "Disabling AMFI for atom type ",iAtmNr(iAtm)
               iAtmNr2(iAtm) = 0
               Charge2(iAtm) = 0.0d0
             end if
           end do
         End do

         Call Gen_RelPointers(-(Info-1))
         Call Drv_AMFI(Label,ipList,OperI,nComp,rHrmt,
     &                 OperC, iAtmNr2, Charge2,DInf,nDInf)
         Call Gen_RelPointers(Info-1)

         Call Deallocate_Auxiliary()
      End If
************************************************************************
************************************************************************
*KAMAL)                                                                *
*              GEN1INT                                                 *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
************************************************************************
      If (lPSOI.and..Not.Prprt.and..Not.Primitive_Pass) Then
#ifdef _GEN1INT_
         PLabel=' '
         rHrmt=-One
         nComp = 3
         Call Get_nAtoms_All(nAtoms)
         Do iCnt=1, nAtoms
         Write (Label,'(A,I2)') 'PSOI  ', iCnt
         nPSOI=iCnt
         nOrdOp = 1
         Call Allocate_Auxiliary()

         iComp=0
         do iComp= 1,nComp
c FIXME ipPSO is uninitialized
c        Call DCopy_(3,Work(ipPSO),1,CoorO(1+(iComp-1)*3),1)
         Call SysAbendMsg('Drv1El',
     &        'Faulty code (undefined Work index).',
     &        'Please correct it or contact the developer.')
         enddo

         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
*
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iChOx = iChBas(3) + iChBas(4)
         OperI(1  ) = iSymLx
         OperC(1  ) = iChOx
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iChOy = iChBas(4) + iChBas(2)
         OperI(1+1) = iSymLy
         OperC(1+1) = iChOy
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
         iChOz = iChBas(2) + iChBas(3)
         OperI(1+2) = iSymLz
         OperC(1+2) = iChOz

*        Zero nuclear contribution
         Call DCopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(PSOInt,PSOMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
         enddo
!1555      Continue
          !Call PrMtrx(Label,lOper,nComp,ip)
#else
         Call WarningMessage(2,
     &   'Drv1El: NO Gen1int interface available!')
         Call Abend()
#endif
       End If   ! lOAM
************************************************************************
************************************************************************
*20)                                                                   *
*     The gradient of the kinetic energy and the nuclear attraction    *
*     energy with respect to the magnetic field.                       *
*                                                                      *
************************************************************************
************************************************************************
      If (GIAO.and..Not.Primitive_Pass) Then
         PLabel=' '
         rHrmt=-One
         nOrdOp = 0
         nComp = 3
         Call Allocate_Auxiliary()
         Call dcopy_(3*nComp,[Zero],0,CoorO,1)
         ixyz=1
         OperI(1  ) = 2**IrrFnc(ixyz)
         OperC(1  ) = iChBas(2)
         ixyz=2
         OperI(1+1) = 2**IrrFnc(ixyz)
         OperC(1+1) = iChBas(3)
         ixyz=4
         OperI(1+2) = 2**IrrFnc(ixyz)
         OperC(1+2) = iChBas(4)
*
         Call dcopy_(3,[Zero],0,Nuc,1)
*
         Label='dT/dB   '
         Call OneEl(KneInt_GIAO,KneMem_GIAO,Label,
     &              ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         nOrdOp = 1
         Label='dV/dB   '
         Call OneEl( NAInt_GIAO, NAMem_GIAO,Label,
     &              ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Auxiliary()
*
*        Differentiate the generalized kinetic energy operator with respect
*        to the magnetic moment at the centers.
*
         nOrdOp = 1
         nComp = 3
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         ixyz=3
         iSymXY = 2**IrrFnc(ixyz)
         ixyz=5
         iSymXZ = 2**IrrFnc(ixyz)
         ixyz=6
         iSymYZ = 2**IrrFnc(ixyz)
         ixyz=7
         iSyXYZ = 2**IrrFnc(ixyz)
*
         iEF = 0
         Do  iCnttp = 1, nCnttp
            Do iCnt = 1, nCntr(iCnttp)
               iEF=iEF+1
               Write (Label,'(A,I2)') 'dT/dmu',iEF
               Call dcopy_(3,Work(ipCntr(iCnttp)+(iCnt-1)*3),1,Ccoor,1)
               Call Allocate_Auxiliary()
               iSymC = 1
               If (Ccoor(1).ne.Zero) iSymC = iOr(iSymC,iSymX)
               If (Ccoor(2).ne.Zero) iSymC = iOr(iSymC,iSymY)
               If (Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSymZ)
               If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &            iSymC = iOr(iSymC,iSymXY)
               If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &            iSymC = iOr(iSymC,iSymXZ)
               If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &            iSymC = iOr(iSymC,iSymYZ)
               If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &             Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSyXYZ)
*
               iComp=0
               Do ix = nOrdOp, 0, -1
                  Do iy = nOrdOp-ix, 0, -1
                     iComp = iComp + 1
                     iz = nOrdOp-ix-iy
                     ixyz=0
                     If (Mod(ix,2).ne.0) ixyz=iOr(ixyz,1)
                     If (Mod(iy,2).ne.0) ixyz=iOr(ixyz,2)
                     If (Mod(iz,2).ne.0) ixyz=iOr(ixyz,4)
                     iSym = 2**IrrFnc(ixyz)
                     If (Ccoor(iComp).ne.Zero ) iSym = iOr(iSym,1)
                     OperI(1+(iComp-1)) = MltLbl(iSymC,iSym,nIrrep)
                     OperC(1+(iComp-1)) = iChBas(iComp+1)
                     Call dcopy_(3,Ccoor,1,CoorO(1+(iComp-1)*3),1)
                  End Do
               End Do
*
*              Call EFNuc(CoorO,Chrg,Centr,kCentr,
*    &                    Nuc,nOrdOp)
*

               Call OneEl(dTdmu_Int,dTdmu_Mem,Label,ipList,
     &                    OperI,nComp,
     &                    CoorO,nOrdOp,Nuc,rHrmt,
     &                    OperC,dum,1,dum,idum,0,0,
     &                    dum,1,0)
*
               Call Deallocate_Auxiliary()
            End Do
         End Do
      End If
************************************************************************
************************************************************************
*                                                                      *
* 21) Atomic Fock matrix                                               *
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
         Call dcopy_(3,[Zero],0,CoorO,1)
         OperI(1) = 1
         OperC(1) = iChBas(1)
*
         Label='FckInt  '
         Call Drv_Fck(Label,ipList,OperI,nComp,
     &                CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &                dum,1,dum,idum,0,0,
     &                dum,1,0,DInf,nDInf)
*
         Call Deallocate_Auxiliary()
      End If
      Call Gen_RelPointers(Info-1)
*                                                                      *
************************************************************************
*                                                                      *
*
*     Deallocate memory for property calculation.
*
      If (Prprt) Then
         If (Short) Call mma_deallocate(Den)
         Call mma_deallocate(Vec)
         Call mma_deallocate(Occ)
         Call CollapseOutput(0,'   Molecular properties:')
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Is this the second pass in a relativistic calculation?
*     In that case: close ONEINT, open ONEREL, read the pVp
*     integrals, close ONEREL, re-open ONEINT, calculate the
*     DK integrals in the contracted basis set, transform and
*     replace the one-electron integrals on ONEINT.
*
      If (DKroll.and..Not.Primitive_Pass) then
      If (BSS) Then
            Call BSSint
         Else
            Call DKRelInt_DP
         End If
      End If
************************************************************************
************************************************************************
*20)                                                                   *
*     Produce information for the NEMO interface.                      *
*                                                                      *
************************************************************************
************************************************************************
      If (NEMO) Then
         If (Primitive_Pass) Then
*
*           Compute p-matrix and put it temporarily on ONEREL.
            PLabel=' '
            rHrmt=One
            nComp=3
            nOrdOp = 0
            Call Allocate_Auxiliary()
            Do iComp = 1, nComp
               Call dcopy_(3,[Zero],0,CoorO(1+(iComp-1)*3),1)
               OperI(1+(iComp-1)) = 1
               OperC(1+(iComp-1)) = iChBas(1)
            End Do
*
            Label='P_matrix'
            Call OneEl(P_Int,P_Mem,Label,ipList,OperI,nComp,
     &                 CoorO,nOrdOp,[Zero,Zero,Zero],rHrmt,OperC,
     &                 dum,1,dum,idum,0,0,
     &                 dum,1,0)
            Call Deallocate_Auxiliary()
         Else
*
*-----------Assemble transformation matrix between the contracted and
*           the primitive basis.
*
            Call NEMO_Opt1()
*
         End If
      End If
************************************************************************
************************************************************************
*21)                                                                   *
*     Fragment AIEMP integrals: projection and 2-electron interaction  *
*                               integrals contracted with the          *
*                               fragment's density matrices            *
*                                                                      *
************************************************************************
************************************************************************
      If(lFAIEMP.and..not.Primitive_Pass) Then
* projection integrals
        PLabel=' '
        rHrmt=One
        nComp=1
        nOrdOp = 0
        Call Allocate_Auxiliary()
        Call dcopy_(3,[Zero],0,CoorO,1)
        OperI(1) = 1
        OperC(1) = iChBas(1)
        Label='FragProj'
        Call OneEl(FragPInt,FragPMem,Label,ipList,OperI,nComp,
     &             CoorO,nOrdOp,[Zero],rHrmt,OperC,
     &             dum,1,dum,idum,0,0,
     &             dum,1,0)
        Call Deallocate_Auxiliary()
* add the results to the one-electron hamiltonian
        iOpt = 0
        iRC = -1
        Call mma_allocate(FragP,n2Tri(1)+4,label='FragP')
        Call RdOne(iRC,iOpt,Label,1,FragP,lOper)
        If (iRC.ne.0) Then
           Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
           Call Quit(_RC_IO_ERROR_READ_)
        End If
        Label = 'OneHam  '
        iRC = -1
        Call mma_allocate(OneHam,n2Tri(1)+4,label='OneHam')
        Call RdOne(iRC,iOpt,Label,1,OneHam,lOper)
        If (iRC.ne.0) Then
           Call WarningMessage(2,
     &                  'Drv1El: Error reading ONEINT;'
     &                //'Label='//Label)
           Call Quit(_RC_IO_ERROR_READ_)
        End If
        Call DaXpY_(n2Tri(1)+4,One,FragP,1,OneHam,1)
        iRC = -1
        Call WrOne(iRC,iOpt,Label,1,OneHam,lOper)
        If (iRC.ne.0) Then
           Call WarningMessage(2,
     &                  'Drv1El: Error writing ONEINT;'
     &                //'Label='//Label)
           Call Quit(_RC_IO_ERROR_WRITE_)
        End If
        iRC = -1
        Label = 'OneHam 0'
        Call WrOne(iRC,iOpt,Label,1,OneHam,lOper)
        If (iRC.ne.0) Then
           Call WarningMessage(2,
     &                  'Drv1El: Error writing ONEINT;'
     &                //'Label='//Label)
           Call Quit(_RC_IO_ERROR_WRITE_)
        End If
        Call mma_deallocate(FragP)
        Call mma_deallocate(OneHam)
* 2-electron interaction integrals (are added to the one-electron
* hamiltonian locally)
        Call Drv2El_FAIEMP()
      End If
*                                                                      *
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
      End Subroutine Drv1el
