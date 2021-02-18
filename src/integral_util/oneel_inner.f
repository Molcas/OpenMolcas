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
* Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      Subroutine OneEl_Inner
     &                 (Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,
     &                  nOrdOp,rHrmt,iChO,
     &                  opmol,opnuc,ipad,iopadr,idirect,isyop,
     &                  iStabO,nStabO,nIC,
     &                  PtChrg,nGrid,iAddPot,Array,LenTot)
************************************************************************
*                                                                      *
* Object: to compute the one-electron integrals. The method employed at*
*         this point is not necessarily the fastest. However, the total*
*         time for the computation of integrals will depend on the time*
*         spent in computing the two-electron integrals.               *
*         The memory at this point is assumed to be large enough to do *
*         the computation in core.                                     *
*         The data is structured with respect to four indices, two (my *
*         ny or i j) refer to primitives or basis functions and two (a *
*         b) refer to the components of the cartesian or spherical     *
*         harmonic gaussians.                                          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*             Modified for Hermite-Gauss quadrature November '90       *
*             Modified for Rys quadrature November '90                 *
*             Modified for multipole moments November '90              *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for general kernel routines January  91         *
*             Modified for nonsymmetrical operators February  91       *
*             Modified for better symmetry treatement October  93      *
*             Modified loop structure April 99                         *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info, only: dbsc
      use Sizes_of_Seward, only: S
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm, Rsv_Tsk
C     Logical Addpot
#include "real.fh"
#include "rmat_option.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
#include "property_label.fh"
      Real*8 Array(LenTot)
      Real*8, Dimension(:), Allocatable :: Zeta, ZI, Kappa, PCoor,
     &                                     SOInt, FArray, Scrtch,
     &                                     ScrSph, Kern
      Integer, Dimension(:,:), Allocatable :: Ind_ij
      Real*8 CCoor(3,nComp), PtChrg(nGrid)
      dimension opmol(*),opnuc(*),iopadr(nComp,*)
      Character Label*8
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      Logical Do_PGamma, Rsv_Tsk
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 112
      iPrint = nPrint(iRout)
      RMat_type_integrals=.False.
      Do_PGamma = .True.
*
*-----Auxiliary memory allocation.
*
      Call mma_allocate(Zeta,S%m2Max,label='Zeta')
      Call mma_allocate(ZI,S%m2Max,label='ZI')
      Call mma_allocate(Kappa,S%m2Max,label='Kappa')
      call mma_allocate(PCoor,S%m2Max*3,label='PCoor')
*                                                                      *
************************************************************************
*                                                                      *
      Call Nr_Shells(nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      Call mma_allocate(Ind_ij,2,nskal*(nSkal+1)/2,label='Ind_ij')
      nijS = 0
      is = 0
      js = 0
      Do I = 1,nSkal*(nSkal+1)/2
         nijS = nijS + 1
         js = js + 1
         If (jS .gt. iS) Then
            iS = jS
            jS = 1
         End If
         Ind_ij(1,nijS)=iS
         Ind_ij(2,nijS)=jS
      End Do
      Call Init_Tsk(id_Tsk,nijS)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate scratch for the integral evaluation.
*
      lFinal=1
      lScrt1=1
      lScrt2=1
      MemKrn=1
      Do ijS = 1, nijS
         iS=Ind_ij(1,ijS)
         jS=Ind_ij(2,ijS)
         iPrim=iSD(5,iS)
         jPrim=iSD(5,jS)
         iBas=iSD(3,iS)
         jBas=iSD(3,jS)
         iAng=iSD(1,iS)
         jAng=iSD(1,jS)
*
         mFinal=nIC*iPrim*jPrim*nElem(iAng)*nElem(jAng)
         lFinal=Max(lFinal,mFinal)
*
         If (Label(1:4).eq.'PSOI') Cycle
         mScrt1=nIC*Max(iPrim,jBas)*Max(iBas,jPrim)
     &         *nElem(iAng)*nElem(jAng)
         lScrt1=Max(mScrt1,lScrt1)
*
         mScrt2=nIC*iBas*jBas*nElem(iAng)*nElem(jAng)
         lScrt2=Max(mScrt2,lScrt2)
*
         Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
*
         If (PLabel.ne.' ') Then
            la0=iAng
            lb0=jAng
            MemAux= 1 + 3*nElem(la0)*nElem(lb0+1)*nIC
            la1=la0
            lb1=lb0+1
            MemBux= 1 + 3*nElem(la1+1)*nElem(lb1)*nIC
            If (la1.ne.0) MemBux=MemBux+3*nElem(la1-1)*nElem(lb1)*nIC
            If (lb0.ne.0) Then
               lb1=lb0-1
               MemAux=MemAux+3*nElem(la0)*nElem(lb0-1)*nIC
               MemCux=1+3*nElem(la1+1)*nElem(lb1)*nIC
               If (la1.ne.0) MemCux=MemCux+3*nElem(la1-1)*nElem(lb1)*nIC
            Else
               MemCux=0
            End If
            MemAux = MemAux + Max(MemBux,MemCux)
            MemKer = MemKer + MemAux
         End If
*
         MemKrn=Max(MemKer*iPrim*jPrim,MemKrn)
      End Do
*
      Call mma_Allocate(FArray,lFinal,label='Final')
      Call mma_allocate(Scrtch,lScrt1,label='Scrtch')
      Call mma_allocate(ScrSph,lScrt2,label='ScrSph')
      Call mma_allocate(Kern,MemKrn,label='Kern')
*                                                                      *
************************************************************************
*                                                                      *
*     big loop over individual tasks, distributed over individual nodes
      ijSh = 0
 10   Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_Tsk(id_Tsk,ijSh)) Go To 11
      iS = Ind_ij(1,ijSh)
      jS = Ind_ij(2,ijSh)

      iCmp  =iSD(2,iS)
      iBas  =iSD(3,iS)
      iAO   =iSD(7,iS)
      iShell=iSD(11,iS)
      iCnttp=iSD(13,iS)

      jCmp  =iSD(2,jS)
      jBas  =iSD(3,jS)
      jAO   =iSD(7,jS)
      jShell=iSD(11,jS)
      jCnttp=iSD(13,jS)

      nSO=0
      Do iComp = 1, nComp
         iSmLbl=lOper(iComp)
         nSO=nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
      End Do
      If (iPrint.ge.29) Write (6,*) ' nSO=',nSO
*
*     Do not compute matrix elements in which electronic and
*     muonic basis sets are mixed.
*
      If (nSO.gt.0 .AND.
     &   dbsc(iCnttp)%fMass.eq.dbsc(jCnttp)%fMass
     &   ) Then
         l_SOInt=iBas*jBas*nSO
         Call mma_allocate(SOInt,l_SOInt,label='SOInt')
         SOInt(:)=Zero
         ipSO=1
         Call OneEl_IJ(iS,jS,iPrint,Do_PGamma,
     &                 Zeta,ZI,Kappa,PCoor,
     &                 Kernel,KrnlMm,Label,lOper,nComp,CCoor,
     &                 nOrdOp,iChO,
     &                 iStabO,nStabO,nIC,
     &                 PtChrg,nGrid,iAddPot,SOInt,l_SOInt,
     &                 FArray,lFinal,Scrtch,lScrt1,
     &                 ScrSph,lScrt2,Kern,MemKrn)
         iSOBlk = ipSO
         Do iComp = 1, nComp
            iSmLbl=lOper(iComp)
            If (n2Tri(iSmLbl).ne.0) Then
               mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            Else
               mSO=0
            End If
*
*           Special trick for integrals over electromagnetic field
*           radiation integrals.
*
            rHrmt_Save=rHrmt

            If (Label(1:5).eq.'EMFR '.or.Label(1:5).eq.'TMOM ') Then
               If (MOD((iComp+5),6).lt.3) Then
                  rHrmt= One
               Else
                  rHrmt=-One
               End If
            End If
*           Write (*,*) 'Label,iComp,rHrmt=',Label,iComp,rHrmt
            If (mSO.ne.0) Then
               Call SOSctt(SOInt(iSOBlk),iBas,jBas,mSO,Array(ip(iComp)),
     &                     n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,
     &                     jShell,iAO,jAO,nComp,Label,lOper,rHrmt)
               iSOBlk = iSOBlk + mSO*iBas*jBas
            End If
            rHrmt=rHrmt_Save
         End Do
         Call mma_deallocate(SOInt)
      End If
      Goto 10
   11 Continue
      Call Free_Tsk(id_Tsk)
      Do iComp = 1, nComp
         iSmLbl=lOper(iComp)
         Call GADSum(Array(ip(iComp)),n2Tri(iSmLbl))
      End Do
*
      Call mma_deallocate(Kern)
      Call mma_deallocate(ScrSph)
      Call mma_deallocate(Scrtch)
      Call mma_deallocate(FArray)
      Call mma_deallocate(Ind_ij)
      Call mma_deallocate(PCoor)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(ZI)
      Call mma_deallocate(Zeta)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(opmol)
         Call Unused_real_array(opnuc)
         Call Unused_integer(ipad)
         Call Unused_integer_array(iopadr)
         Call Unused_integer(idirect)
         Call Unused_integer(isyop)
      End If
      End
