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
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OneEl                                                   *
*              QExit                                                   *
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
      Implicit Real*8 (A-H,O-Z)
      External OvrGrd, KneGrd, NAGrd, PrjGrd, M1Grd, M2Grd, SROGrd,
     &         WelGrd, XFdGrd, RFGrd, PCMGrd, PPGrd, COSGrd, FragPGrd
      External OvrMmG, KneMmG, NAMmG, PrjMmG, M1MmG, M2MmG, SROMmG,
     &         WelMmg, XFdMmg, RFMmg, PCMMmg, PPMmG, FragPMmG
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
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
CAOM>
      Logical DiffOp
*
*-----Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
*...  Prologue
      iRout = 131
      iPrint = nPrint(iRout)
      Call CWTime(TCpu1,TWall1)
      Call qEnter('Drvh1')
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
*
*
*...  Get the method label
*     print *,' Read Method label'
      Call Get_cArray('Relax Method',Method,8)
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*     print *,' Read density matrix'
      Call Get_D1ao_Var(ipD_var,Length)
      If ( length.ne.nDens ) Then
         Call WarningMessage(2,'Error in Drvh1')
         Write (6,*) 'Drvh1: length.ne.nDens'
         Write (6,*) 'length,nDens=',length,nDens
         Call Abend()
      End If
      If (iPrint.ge.99) then
         Write(6,*) 'variational 1st order density matrix'
         ii=ipD_Var
         Do iIrrep = 0, nIrrep - 1
            Write(Label,*) 'symmetry block',iIrrep
            Call TriPrt(Label,' ',Work(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*
*...  Read the generalized Fock matrix
*...  Fock matrix in AO/SO basis
*     print *,' Read Fock matrix'
      If (.Not.HF_Force) Then
         Call Get_Fock_Occ(ipFock,Length)
         If ( length.ne.nDens ) Then
            Call WarningMessage(2,'Error in Drvh1')
            Write (6,*) 'Drvh1: length.ne.nDens'
            Write (6,*) 'length,nDens=',length,nDens
            Call Abend()
         End If
!         iprint=100
         If (iPrint.ge.99) then
            Write(6,*) 'generalized Fock matrix'
            ii=ipFock
            Do iIrrep = 0, nIrrep - 1
               Write(Label,*) 'symmetry block',iIrrep
               Call TriPrt(Label,' ',Work(ii),nBas(iIrrep))
               ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
            End Do
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     nOrdOp: order/rank of the operator
*     Work(ip1): lOper of each component of the operator
*
      nOrdOp=0
      nComp = nElem(nOrdOp)
      Call GetMem('Coor','Allo','Real',ipC,3*nComp)
      Call GetMem('lOper','Allo','Inte',ip1,nComp)
      call dcopy_(nComp*3,Zero,0,Work(ipC),1)
      iWork(ip1) = 1
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
      Call OneEl_g(OvrGrd,OvrMmG,Temp,nGrad,DiffOp,Work(ipC),
     &           Work(ipFock),nFock,iWork(ip1),nComp,nOrdOp,Label)
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
      Call OneEl_g(KneGrd,KneMmG,Temp,nGrad,DiffOp,Work(ipC),
     &           Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
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
        Call GetMem('lOperf','Allo','Inte',ip1f,nCompf)
        do ii=1,ncompf
          iwork(ip1f+ii-1)=1
        enddo
        DiffOp=.False.
        if(nOrdOpf.gt.0) DiffOp=.True.
        Call OneEl_g(MltGrd,MltMmG,Temp,nGrad,DiffOp,Work(ipC),
     &           Work(ipD_Var),nDens,iWork(ip1f),nCompf,nOrdOpf,Label)
        Call MltGrdNuc(Temp,nGrad,nOrdOpf)
        Call GetMem('lOperf','Free','Inte',ip1f,nCompf)
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
      Call OneEl_g(NAGrd,NAMmG,Temp,nGrad,DiffOp,Work(ipC),
     &           Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      If (HF_Force) Then
         If (lECP.or.nWel.ne.0.or.lXF.or.lRF) Then
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
         Call OneEl_g(PrjGrd,PrjMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M1 Operator Contribution'
         Call OneEl_g( M1Grd, M1MmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The M2 Operator Contribution'
         Call OneEl_g( M2Grd, M2MmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
         Label = ' The SR Operator Contribution'
         Call OneEl_g(SROGrd,SROMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
      End If
      If (lPP) Then
         Label = ' The Pseudo Potential Contribution'
         Call OneEl_g(PPGrd,PPMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
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
      Do iWel = 0, nWel-1
         r0   = Work(ipWel+iWel*3  )
         ExpB = Work(ipWel+iWel*3+1)
         Label = ' The Spherical Well Contribution'
         Call OneEl_g(WelGrd,WelMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Fact = Work(ipWel+iWel*3+2)
         Call DaXpY_(nGrad,Fact,Temp,1,Grad,1)
      End Do
************************************************************************
*6)                                                                    *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the external field integrals.                        *
*                                                                      *
************************************************************************
      If (lXF) Then
         DiffOp = .True.
         Label = ' The External Field Contribution'
         Call OneEl_g(XFdGrd,XFdMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
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
         Call GetMem('lOper','Free','Inte',ip1,nComp)
         Call GetMem('Coor','Free','Real',ipC,3*nComp)
*
*------- The Kirkwood model
*
         nOrdOp=lMax
         nComp=(lMax+1)*(lMax+2)*(lMax+3)/6
         Call GetMem('lOper','Allo','Inte',ip1,nComp)
*
*------- Store permutation symmetry of components of the EF
*
         iComp = 0
         Do iMltpl = 0, lMax
            Do ix = iMltpl, 0, -1
               If (Mod(ix,2).eq.0) Then
                  iSymX=1
               Else
                  ixyz=1
                  iSymX=2**IrrFnc(ixyz)
               End If
               Do iy = iMltpl-ix, 0, -1
                  If (Mod(iy,2).eq.0) Then
                     iSymY=1
                  Else
                     ixyz=2
                     iSymY=2**IrrFnc(ixyz)
                  End If
                  iz = iMltpl-ix-iy
                  If (Mod(iz,2).eq.0) Then
                     iSymZ=1
                  Else
                     ixyz=4
                     iSymZ=2**IrrFnc(ixyz)
                  End If
*                 iWork(ip1+iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ,
*    &                                      nIrrep),nIrrep)
*-----------------Compute only total symmetric contributions
                  iWork(ip1+iComp) = 1
                  iComp = iComp + 1
               End Do
            End Do
         End Do
         Call GetMem('Coor','Allo','Real',ipC,3*nComp)
         call dcopy_(nComp*3,Zero,0,Work(ipC),1)
         DiffOp = .True.
         Label = ' The Electronic Reaction Field Contribution'
         Call OneEl_g(RFGrd,RFMmG,Temp,nGrad,DiffOp,Work(ipC),
     &              Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*
      Else If (lRF.and.PCM) Then
         iCOSMO=0
*------- The PCM / COSMO model
*
         If (iCOSMO.le.0) Then
            iPrint=15
            Call DScal_(nTs*2,One/DBLE(nIrrep),Work(ip_Q),1)
         End If
*        iWork(ip1) = 255
         iWork(ip1) = 1
         DiffOp = .True.
         If (iCOSMO.gt.0) Then
            Call dzero(Temp,ngrad)
            Label= ' The Electronic Reaction Field Contribution (COSMO)'
            Call OneEl_g(COSGrd,PCMMmG,Temp,nGrad,DiffOp,Work(ipC),
     &                   Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,
     &                   Label)
            If (iPrint.ge.15) Then
               Label=' Reaction Field (COSMO) Contribution'
               Call PrGrad(Label,Temp,nGrad,lIrrep,ChDisp,5)
            End If
         Else
            Label = ' The Electronic Reaction Field Contribution (PCM)'
            Call OneEl_g(PCMGrd,PCMMmG,Temp,nGrad,DiffOp,Work(ipC),
     &                   Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,
     &                   Label)
            If (iPrint.ge.15) Then
               Label=' Reaction Field (PCM) Contribution'
               Call PrGrad(Label,Temp,nGrad,lIrrep,ChDisp,5)
            End If
         End If

         Call DaXpY_(nGrad,One,Temp,1,Grad,1)
         If (iCOSMO.eq.0) Call DScal_(nTs*2,DBLE(nIrrep),Work(ip_Q),1)
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
        Call OneEl_g(FragPGrd,FragPMmG,Temp,nGrad,DiffOp,Work(ipC),
     &               Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
        Call DaXpY_(nGrad,One,Temp,1,Grad,1)
        Call DrvG_FAIEMP(Grad,Temp,nGrad)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('lOper','Free','Inte',ip1,nComp)
      Call GetMem('Coor','Free','Real',ipC,3*nComp)
*                                                                      *
************************************************************************
1099  Continue
*                                                                      *
*...  Epilogue, end
*
      If (.Not.HF_Force) Call GetMem('Fock','Free','Real',ipFock,nFock)
      Call GetMem('D0  ','Free','Real',ipD_Var,nDens)
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      Call qExit('Drvh1')
      Return
      End
