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
* Copyright (C) 1991,1995, Roland Lindh                                *
*               Markus P. Fuelscher                                    *
************************************************************************
      SubRoutine Drvh1(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
* Object: driver for computation of gradient with respect to the one-  *
*         electron hamiltonian and the overlap matrix. The former will *
*         be contracted with the "variational" first order density     *
*         matrix and the latter will be contracted with the generalized*
*         Fock matrix.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
*                                                                      *
* modified by M.P. Fuelscher                                           *
* make use of the RELAX file                                           *
*                                                                      *
*             Modified to Sperical Well, External Field, and SRO       *
*             gradients, April '95. R. Lindh                           *
*             Modified to Self Consistent Reaction Fields, May '95     *
************************************************************************
      use PCM_arrays, only: PCM_SQ
      use External_Centers
      use Basis_Info, only: nCnttp, dbsc, nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External OvrGrd, KneGrd, NAGrd, PrjGrd, M1Grd, M2Grd, SROGrd,
     &         WelGrd, XFdGrd, RFGrd, PCMGrd, PPGrd, COSGrd, FragPGrd
      External OvrMmG, KneMmG, NAMmG, PrjMmG, M1MmG, M2MmG, SROMmG,
     &         WelMmg, XFdMmg, RFMmg, PCMMmg, PPMmG, FragPMmG
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "wldata.fh"
#include "rctfld.fh"
      Character Method*8, Label*80
      Real*8 Grad(nGrad), Temp(nGrad)
      parameter (lforce=20)
      real*8 force(lforce)
      common /finfld/force
      character*30 fldname
      External MltGrd,MltMmG
      Real*8, Allocatable:: Coor(:,:), Fock(:), D_Var(:)
      Integer, Allocatable:: lOper(:), lOperf(:)
CAOM>
      Logical DiffOp, lECP, lPP, lFAIEMP
*
*-----Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
*...  Prologue
      iRout = 131
      iPrint = nPrint(iRout)
      Call CWTime(TCpu1,TWall1)
      Call StatusLine(' Alaska:',' Computing 1-electron gradients')
*
*---- Allocate memory for density and Fock matrices
*
      nFock = 0
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nFock = nFock + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
      lECP=.False.
      lPP =.False.
      lFAIEMP =.False.
      Do i = 1, nCnttp
         lECP = lECP .or. dbsc(i)%ECP
         lPP  = lPP  .or. dbsc(i)%nPP.ne.0
         lFAIEMP = LFAIEMP .or. dbsc(i)%Frag
      End Do
*
*
*...  Get the method label
*     print *,' Read Method label'
      Call Get_cArray('Relax Method',Method,8)
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*     print *,' Read density matrix'

      Call mma_allocate(D_Var,nDens,Label='D_Var')
      Call Get_D1ao_Var(D_Var,nDens)
      If (iPrint.ge.99) then
         Write(6,*) 'variational 1st order density matrix'
         ii=1
         Do iIrrep = 0, nIrrep - 1
            Write(Label,*) 'symmetry block',iIrrep
            Call TriPrt(Label,' ',D_Var(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*
*...  Read the generalized Fock matrix
*...  Fock matrix in AO/SO basis
*     print *,' Read Fock matrix'
      If (.Not.HF_Force) Then
         Call mma_allocate(Fock,nDens,Label='Fock')
         Call Get_Fock_Occ(Fock,nDens)
         If (iPrint.ge.99) then
            Write(6,*) 'generalized Fock matrix'
            ii=1
            Do iIrrep = 0, nIrrep - 1
               Write(Label,*) 'symmetry block',iIrrep
               Call TriPrt(Label,' ',Fock(ii),nBas(iIrrep))
               ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
            End Do
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     nOrdOp: order/rank of the operator
*     lOper: lOper of each component of the operator
*
      nOrdOp=0
      nComp = nElem(nOrdOp)
      Call mma_allocate(Coor,3,nComp,Label='Coor')
      Call mma_allocate(lOper,nComp,Label='lOper')
      Coor(:,:)=Zero
      lOper(1)=1
      If (HF_Force) Go To 1003
************************************************************************
*1)                                                                    *
*     Trace the generalized Fock matrix with the gradient of the       *
*     overlap matrix.                                                  *
*                                                                      *
************************************************************************
*
      DiffOp = .False.
      Label  = ' The Renormalization Contribution'
      Call OneEl_g(OvrGrd,OvrMmG,Temp,nGrad,DiffOp,Coor,
     &             Fock,nFock,lOper,nComp,nOrdOp,Label)
      Call DaXpY_(nGrad,-One,Temp,1,Grad,1)
*
************************************************************************
*2)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the kinetic energy integrals.                        *
*                                                                      *
************************************************************************
*
      DiffOp = .False.
      Label  = ' The Kinetic Energy Contribution'
      Call OneEl_g(KneGrd,KneMmG,Temp,nGrad,DiffOp,Coor,
     &             D_Var,nDens,lOper,nComp,nOrdOp,Label)
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
CAOM<
C... Check finite field operators...
      nextfld=0
 1100 call GetNextFfield(nextfld,fldname,nOrdOpf,ncmp,force,lforce)
      if(nextfld.ne.0) then
        if(nOrdOpf.eq.0) then
          if(fldname(1:2).ne.'OV') then
            Call WarningMessage(2,'Error in Drvh1')
            Write(6,*) 'Finite field gradients only for MLTP'
            Call Quit_OnUserError()
          endif
        endif
        nCompf = nElem(nOrdOpf)
        Label=fldname
        if(ncompf.ne.ncmp) then
           Call WarningMessage(2,'Error in Drvh1')
           Write(6,*) 'Wrong number of components in FF grad'
           Call Quit_OnUserError()
        endif
        Call mma_allocate(lOperf,nCompf,Label='lOperf')
        lOperf(:)=1
        DiffOp=.False.
        if(nOrdOpf.gt.0) DiffOp=.True.
        Call OneEl_g(MltGrd,MltMmG,Temp,nGrad,DiffOp,Coor,
     &               D_Var,nDens,lOperf,nCompf,nOrdOpf,Label)
        Call MltGrdNuc(Temp,nGrad,nOrdOpf)
        Call mma_deallocate(lOperf)
        Call DaXpY_(nGrad,-One,Temp,1,Grad,1)
        goto 1100
      endif
CAOM>
*
************************************************************************
*3)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the nuclear attraction integrals.                    *
*                                                                      *
************************************************************************
*
 1003 Continue
      DiffOp = .True.
      Label = ' The Nuclear Attraction Contribution'
      Call OneEl_g(NAGrd,NAMmG,Temp,nGrad,DiffOp,Coor,
     &             D_Var,nDens,lOper,nComp,nOrdOp,Label)
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      If (HF_Force) Then
         If (lECP.or.nWel.ne.0.or.Allocated(XF).or.lRF) Then
            Call WarningMessage(2,'Error in Drvh1')
            Write (6,*) 'HF forces not implemented yet for this case!'
            Call Quit_OnUserError()
         End If
         Go To 1099
      End If
*
************************************************************************
*4)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the ECP integrals.                                   *
*                                                                      *
************************************************************************
*
      If (lECP) Then
         DiffOp = .True.
         Label = ' The Projection Operator Contribution'
         Call OneEl_g(PrjGrd,PrjMmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M1 Operator Contribution'
         Call OneEl_g( M1Grd, M1MmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M2 Operator Contribution'
         Call OneEl_g( M2Grd, M2MmG,Temp,nGrad,DiffOp,Coor,
     &                D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The SR Operator Contribution'
         Call OneEl_g(SROGrd,SROMmG,Temp,nGrad,DiffOp,Coor,
     &              D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
      If (lPP) Then
         Label = ' The Pseudo Potential Contribution'
         Call OneEl_g(PPGrd,PPMmG,Temp,nGrad,DiffOp,Coor,
     &              D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
*
************************************************************************
*5)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the Spherical well integrals.                        *
*                                                                      *
************************************************************************
      DiffOp = .True.
      Do iWel = 1, nWel
         r0   = Wel_Info(1,iWel)
         ExpB = Wel_Info(2,iWel)
         Label = ' The Spherical Well Contribution'
         Call OneEl_g(WelGrd,WelMmG,Temp,nGrad,DiffOp,Coor,
     &              D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Fact = Wel_Info(3,iWel)
         Call DaXpY_(nGrad,Fact,Temp,1,Grad,1)
      End Do
************************************************************************
*6)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the external field integrals.                        *
*                                                                      *
************************************************************************
      If (Allocated(XF)) Then
         DiffOp = .True.
         Label = ' The External Field Contribution'
         Call OneEl_g(XFdGrd,XFdMmG,Temp,nGrad,DiffOp,Coor,
     &              D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
*
************************************************************************
*7)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the reaction field integrals.                        *
*                                                                      *
************************************************************************
*
      If (lRF.and..Not.lLangevin.and..Not.PCM) Then
         Call mma_deallocate(lOper)
         Call mma_deallocate(Coor)
*
*------- The Kirkwood model
*
         nOrdOp=lMax
         nComp=(lMax+1)*(lMax+2)*(lMax+3)/6
         Call mma_allocate(lOper,nComp,Label='lOper')
*
*------- Store permutation symmetry of components of the EF
*
         iComp = 1
         Do iMltpl = 0, lMax
            Do ix = iMltpl, 0, -1
               If (Mod(ix,2).eq.0) Then
*                 iSymX=1
               Else
                  ixyz=1
*                 iSymX=2**IrrFnc(ixyz)
               End If
               Do iy = iMltpl-ix, 0, -1
                  If (Mod(iy,2).eq.0) Then
*                    iSymY=1
                  Else
                     ixyz=2
*                    iSymY=2**IrrFnc(ixyz)
                  End If
                  iz = iMltpl-ix-iy
                  If (Mod(iz,2).eq.0) Then
*                    iSymZ=1
                  Else
                     ixyz=4
*                    iSymZ=2**IrrFnc(ixyz)
                  End If
*                 lOper(iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ))
*-----------------Compute only total symmetric contributions
                  lOper(iComp) = 1
                  iComp = iComp + 1
               End Do
            End Do
         End Do
         Call mma_allocate(Coor,3,nComp,Label='Coor')
         Coor(:,:)=Zero
         DiffOp = .True.
         Label = ' The Electronic Reaction Field Contribution'
         Call OneEl_g(RFGrd,RFMmG,Temp,nGrad,DiffOp,Coor,
     &              D_Var,nDens,lOper,nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
      Else If (lRF.and.PCM) Then
         iCOSMO=0
*------- The PCM / COSMO model
*
         If (iCOSMO.le.0) Then
            iPrint=15
            Call DScal_(nTs*2,One/DBLE(nIrrep),PCM_SQ,1)
         End If
         lOper(1) = 1
         DiffOp = .True.
         If (iCOSMO.gt.0) Then
            Call dzero(Temp,ngrad)
            Label= ' The Electronic Reaction Field Contribution (COSMO)'
            Call OneEl_g(COSGrd,PCMMmG,Temp,nGrad,DiffOp,Coor,
     &                   D_Var,nDens,lOper,nComp,nOrdOp,
     &                   Label)
            If (iPrint.ge.15) Then
               Label=' Reaction Field (COSMO) Contribution'
               Call PrGrad(Label,Temp,nGrad,ChDisp,5)
            End If
         Else
            Label = ' The Electronic Reaction Field Contribution (PCM)'
            Call OneEl_g(PCMGrd,PCMMmG,Temp,nGrad,DiffOp,Coor,
     &                   D_Var,nDens,lOper,nComp,nOrdOp,
     &                   Label)
            If (iPrint.ge.15) Then
               Label=' Reaction Field (PCM) Contribution'
               Call PrGrad(Label,Temp,nGrad,ChDisp,5)
            End If
         End If

         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
         If (iCOSMO.eq.0) Call DScal_(nTs*2,DBLE(nIrrep),PCM_SQ,1)
*
      End If
*
************************************************************************
*8)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the FAIEMP integrals.                                *
*                                                                      *
************************************************************************
*
      If(lFAIEMP) Then
        DiffOp = .True.
        Label = ' The FAIEMP Projection Operator Contribution'
        Call OneEl_g(FragPGrd,FragPMmG,Temp,nGrad,DiffOp,Coor,
     &               D_Var,nDens,lOper,nComp,nOrdOp,Label)
        Call DaXpY_(nGrad,One,Temp,1,Grad,1)
        Call DrvG_FAIEMP(Grad,Temp,nGrad)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(lOper)
      Call mma_deallocate(Coor)
*                                                                      *
************************************************************************
1099  Continue
*                                                                      *
*...  Epilogue, end
*
      If (.Not.HF_Force) Call mma_deallocate(Fock)
      Call mma_deallocate(D_Var)
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      Return
      End
