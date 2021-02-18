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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
*               1994, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Dot1El2(Kernel,KrnlMm,Hess,nGrad,DiffOp,CCoor,
     &                 FD,nordop)
************************************************************************
*                                                                      *
* Object: to compute gradients of the one electron integrals.          *
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
*             Modified for general kernel routines January '91         *
*             Modified for nonsymmetrical operators February '91       *
*             Modified for gradients October '91                       *
*             Modified for Hessians by AB   Dec '94                    *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iOper
      use Sizes_of_Seward, only: S
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), Ccoor(3), FD(*),
     &       RB(3), Hess(nGrad)
      Integer iDCRR(0:7), iDCRT(0:7), iStabM(0:7),iCoM(0:7,0:7),
     &           nOp(2),
     &           iStabO(0:7),IndGrd(2,3,3,0:7)
      Logical AeqB, EQ, DiffOp
      Real*8, Allocatable:: Zeta(:), ZI(:), Kappa(:), PCoor(:,:),
     &                      Kern(:), Scrt1(:), Scrt2(:), DAO(:),
     &                      DSOpr(:), DSO(:)
      Logical, External :: TF
*
*     Statement function
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      call dcopy_(nGrad,[Zero],0,Hess,1)
*
*     Auxiliary memory allocation.
*
      Call mma_allocate(Zeta,S%m2Max,Label='Zeta')
      Call mma_allocate(ZI,S%m2Max,Label='ZI')
      Call mma_allocate(Kappa,S%m2Max,Label='Kappa')
      Call mma_allocate(PCoor,S%m2Max,3,Label='PCoor')
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal)
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*     Double loop over shells.
*
      nTasks = nSkal*(nSkal+1)/2
      iS = 0
      jS = 0
      Do ijS = 1, nTasks
         jS = jS + 1
         If (jS.gt.iS) Then
            iS = jS
            jS = 1
         End If
C     Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
C        Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*       Call kernel routine to get memory requirement.
*
        Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
        MemKrn=MemKer*S%m2Max
        Call mma_allocate(Kern,MemKrn,Label='Kern')
*
*       Allocate memory for the final integrals, all in the
*       primitive basis.
*
*
*       Scratch area for contraction step
*
        nScrt1 =  S%MaxPrm(iAng)*S%MaxPrm(jAng) *
     &           nElem(iAng)*nElem(jAng)
        Call mma_allocate(Scrt1,nScrt1,Label='Scrt1')
*
*       Scratch area for the transformation to spherical gaussians
*
        nScrt2=S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
        Call mma_allocate(Scrt2,nScrt2,Label='Scrt2')
*
        nDAO=iPrim*jPrim*nElem(iAng)*nElem(jAng)
        Call mma_allocate(DAO,nDAO,Label='DAO')
*
*       At this point we can compute Zeta.
*
          Call ZXia(Zeta,ZI,
     &              iPrim,jPrim,Shells(iShll)%Exp,
     &                          Shells(jShll)%Exp)
*
            AeqB = iS.eq.jS
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,
     &                     dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
            If (.Not.DiffOp .and. nDCRR.eq.1 .and. EQ(A,B)) Go To 131
*
*-----------Find the stabilizer for A and B
*
            Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,
     &                 dc(mdcj)%iStab,dc(mdcj)%nStab,
     &                 iStabM,nStabM)
*
*          Generate all possible (left) CoSet
*          To the stabil. of A and B
*
            Do 433 i = 0, nIrrep-1
               Do 434 j = 0, nStabM-1
                  iCoM(i,j) = iEor(iOper(i),iStabM(j))
 434            Continue
 433         Continue
*           Order the Coset so we will have the unique ones first
            nMax = 1
            Do 435 j = 1, nIrrep-1
*              Check uniqueness
               Do 436 i = 0, nMax - 1
                  Do 437 ielem = 0, nStabM-1
                     If (iCoM(i,1).eq.iCoM(j,ielem))
     &                  Go To 435
 437               Continue
 436            Continue
*              Move unique CoSet
               nMax = nMax + 1
               Do 438 ielem = 0, nStabM-1
                  iTmp = iCoM(nMax-1,ielem)
                  iCoM(nMax-1,ielem) = iCoM(j,ielem)
                  iCoM(j,ielem) = iTmp
 438            Continue
               If (nMax.eq.nIrrep/nStabM) Go To 439
 435         Continue
 439         Continue
*
*           Allocate memory for the elements of the Fock or 1st order
*           denisty matrix which are associated with the current shell
*           pair.
*
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            If (nSO.eq.0) Go To 131
            Call mma_allocate(DSOpr,nSO*iPrim*jPrim,Label='DSOpr')
            DSOpr(:)=Zero
            Call mma_allocate(DSO,nSO*iPrim*jPrim,Label='DSO')
            DSO(:)=Zero
*
*           Gather the elements from 1st order density / Fock matrix.
*
            Call SOGthr(DSO,iBas,jBas,nSO,FD,
     &                  n2Tri(iSmLbl),iSmLbl,
     &                  iCmp,jCmp,iShell,jShell,
     &                  AeqB,iAO,jAO)
*
*           Project the Fock/1st order density matrix in AO
*           basis on to the primitive basis.
*
*
*           Transform IJ,AB to J,ABi
            Call DGEMM_('T','T',
     &                  jBas*nSO,iPrim,iBas,
     &                  1.0d0,DSO,iBas,
     &                        Shells(iShll)%pCff,iPrim,
     &                  0.0d0,DSOpr,jBas*nSO)
*           Transform J,ABi to AB,ij
            Call DGEMM_('T','T',
     &                  nSO*iPrim,jPrim,jBas,
     &                  1.0d0,DSOpr,jBas,
     &                        Shells(jShll)%pCff,jPrim,
     &                  0.0d0,DSO,nSO*iPrim)
*           Transpose to ij,AB
            Call DGeTmO(DSO,nSO,nSO,iPrim*jPrim,DSOpr,
     &                  iPrim*jPrim)
            Call mma_deallocate(DSO)
*
*
*           Loops over symmetry operations.
*
            nOp(1) = NrOpr(0)
            if(jBas.lt.-999999) write(6,*) 'gcc overoptimization',nDCRR
            Do 140 lDCRR = 0, nDCRR-1
               Call OA(iDCRR(lDCRR),B,RB)
               nOp(2) = NrOpr(iDCRR(lDCRR))
               If (EQ(A,RB).and. (.Not.DiffOp)) Go To 140
*
*
               lloper=1
               Call SOS(iStabO,nStabO,llOper)
               Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
*--------------Compute normalization factor due the DCR symmetrization
*              of the two basis functions and the operator.
*
               iuv = dc(mdci)%nStab*dc(mdcj)%nStab
               FactNd = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LmbdT)
               If (MolWgh.eq.1) Then
                  FactNd = FactNd * DBLE(nIrrep)**2 / DBLE(iuv)
               Else If (MolWgh.eq.2) Then
                  FactNd = sqrt(DBLE(iuv))*nStabO/DBLE(nIrrep*LmbdT)
               End If
*
*
*--------------Desymmetrize the matrix with which we will
*              contracte the trace.
*
               Call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                     iShell,jShell,iShll,jShll,
     &                     iAO,jAO,DAO,iPrim,jPrim,
     &                     DSOpr,nSO,nOp,FactNd)
*
*--------------Project the spherical harmonic space onto the
*              cartesian space.
*
               kk = nElem(iAng)*nElem(jAng)
               If (Shells(iShll)%Transf.or.Shells(jShll)%Transf) Then
*
*-----------------ij,AB --> AB,ij
                  Call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,
     &                        iCmp*jCmp,Scrt1,iCmp*jCmp)
*-----------------AB,ij --> ij,ab
                  Call SphCar(Scrt1,iCmp*jCmp,iPrim*jPrim,
     &                        Scrt2,nScr2,
     &                        RSph(ipSph(iAng)),iAng,
     &                        Shells(iShll)%Transf,
     &                        Shells(iShll)%Prjct,
     &                        RSph(ipSph(jAng)),jAng,
     &                        Shells(jShll)%Transf,
     &                        Shells(jShll)%Prjct,
     &                        DAO,kk)
               End If
*
*--------------Compute kappa and P.
*
               Call Setup1(Shells(iShll)%Exp,iPrim,
     &                     Shells(jShll)%Exp,jPrim,
     &                     A,RB,Kappa,PCoor,ZI)
*
*
               Call Icopy(18*nirrep,[0],0,IndGrd,1)
               kk=0
               Do jIrrep=0,nirrep-1
               Do Jcar=1,3
                iirrep=irrfnc(2**(jcar-1))
                If (iirrep.eq.jirrep) Then
                jj=0
                Do i=0,jirrep-1
                 jj=ldisp(i)+jj
                End do
                nDisp = IndDsp(mdci,jIrrep)-jj
                Do iCar=1,3
                 iComp = 2**(iCar-1)
                 If ( TF(mdci,jIrrep,iComp)) Then
                  ndisp=ndisp+1
                  IndGrd(1,icar,jcar,jIrrep) = kk+nDisp
                 End If
                end do
                kk=kk+ldisp(jirrep)
               End If
               End Do
               End Do
*
               kk=0
               Do jIrrep=0,nirrep-1
               Do Jcar=1,3
                iirrep=irrfnc(2**(jcar-1))
                If (iirrep.eq.jirrep) then
                jj=0
                Do i=0,jirrep-1
                 jj=ldisp(i)+jj
                End do
                nDisp = IndDsp(mdcj,jIrrep)-jj
                Do iCar=1,3
                 iComp = 2**(iCar-1)
                 If ( TF(mdcj,jIrrep,iComp)) Then
                  ndisp=ndisp+1
                  IndGrd(2,icar,jcar,jIrrep) = kk+nDisp
                 End If
                end do
                kk=kk+ldisp(jirrep)
               End If
               End Do
               End Do

*
*--------------Compute gradients of the primitive integrals and
*              trace the result.
*


               Call Kernel(Shells(iShll)%Exp,iPrim,
     &                     Shells(jShll)%Exp,jPrim,
     &                     Zeta,ZI,
     &                     Kappa,Pcoor,
     &                     iPrim*jPrim,
     &                     iAng,jAng,A,RB,nOrder,Kern,
     &                     MemKer,Ccoor,
     &                     nOrdOp,Hess,
     &                     indgrd,DAO,
     &                     mdci,mdcj,nOp,
     &                     iStabM,nStabM)
*
 140        Continue
*
            Call mma_deallocate(DSOpr)
 131        Continue
         Call mma_deallocate(DAO)
         Call mma_deallocate(Scrt2)
         Call mma_deallocate(Scrt1)
         Call mma_deallocate(Kern)
*
C        End Do
C     End Do
      End Do
*
      Call Free_iSD()
*
      Call mma_deallocate(PCoor)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(ZI)
      Call mma_deallocate(Zeta)
*
      Return
      End
