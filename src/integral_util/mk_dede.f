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
* Copyright (C) 1990-1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine mk_DeDe(FD,nFD,mFD,ipOffD,nOffD,ipDeDe,ipD00,MaxDe,
     &                   mDeDe,mIndij,Special_NoSym,DFT_Storage,
     &                   DInf,nDInf,DeDe,nDeDe)
************************************************************************
*                                                                      *
* Object: to decontract, desymmetrize the 1st order density matrix.    *
*         The memory at this point is assumed to be large enough to do *
*         the computation in core.                                     *
*         The data is structured with respect to four indices, two (my *
*         ny or i j) refer to primitives or basis functions and two (a *
*         b) refer to the components of the cartesian or spherical     *
*         harmonic gaussians.                                          *
*                                                                      *
*         The indices are here ordered canonically!!!                  *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy    (ESSL)                                         *
*              DGEMM_   (ESSL)                                         *
*              DGeTMO   (ESSL)                                         *
*              DaXpY    (ESSL)                                         *
*              SOGthr                                                  *
*              Desym1                                                  *
*              DScal    (ESSL)                                         *
*              TriPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*             Modified for Hermite-Gauss quadrature November '90       *
*             Modified for Rys quadrature November '90                 *
*             Modified for multipole moments November '90              *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for general kernel routines January '91         *
*             Modified for nonsymmetrical operators February '91       *
*             Modified for gradients October '91                       *
*             Modified to process 1st order density matrices, Dec. '92 *
************************************************************************
      use Real_Spherical
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "lundio.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 DInf(*), DeDe(nDeDe)
      Real*8, Dimension (:), Allocatable :: Scrt, DAO, DSOp, DSOc, DSO
      Real*8 FD(nFD,mFD)
      Character ChOper(0:7)*3
      Integer    iDCRR(0:7), nOp(2), ipOffD(2+mFD,nOffD)
      Logical AeqB, Special_NoSym, DFT_Storage
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 112
      iPrint = nPrint(iRout)
      Call CWTime(TCpu1,TWall1)
C     Call QEnter('DeDe')
*
      If (iPrint.ge.99) Then
         Write (6,*)
         Write (6,*) ' Differential 1st order density matrix'
         iFD = 1
         Do iIrrep = 0, nIrrep - 1
            Write (6,*)
            Write (6,*) 'iIrrep=',iIrrep
            Do jFD = 1, mFD
               Write (6,*) 'jFD=',jFD
               Write (6,*)
               Call TriPrt(' Diagonal block',' ',
     &                     FD(iFD,jFD),nBas(iIrrep))
            End Do
            iFD = iFD + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*
*     ipD00:
*     MaxDCR: max number of possible pairs
*     MaxDe:
*
      mIndij = 0
      iIrrep = 0
      jOffD = 0
      Inc=3
      If (mFD.eq.2) Inc=4
      Call ICopy(nOffD,[ipD00],0,ipOffD(1,1),Inc)
      If (mFD.eq.2) Call ICopy(nOffD,[ipD00],0,ipOffD(4,1),Inc)
      Call ICopy(nOffD,[MaxDCR],0,ipOffD(2,1),Inc)
      Call ICopy(nOffD,[MaxDe],0,ipOffD(3,1),Inc)
*                                                                      *
************************************************************************
*                                                                      *
      Call Nr_Shells(nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*     i.e. (dd), (dp), (pp), etc. This is the ab index.
*
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iCff   = iSD( 4,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
*
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jCff   = iSD( 4,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            ijShll = iTri(iShell,jShell)
*
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            If (nSO.eq.0) Go To 131
*                                                                      *
************************************************************************
*                                                                      *
            If (iPrint.ge.19) Then
               Write (6,*) 'iS,jS=',iS,jS
               Write (6,'(A,A,A,A,A)')
     &            ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
            End If
*
*---------- Scratch area for contraction step
*
            nScr1 =  MaxPrm(iAng)*MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
*
*---------- Scratch area for the transformation to spherical gaussians
*
            nScr2=MaxPrm(iAng)*MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
*
            Call mma_allocate(DAO,Max(iBas*jBas,iPrim*jPrim)*iCmp*jCmp,
     &                        label='DAO')
*
            AeqB = iS .eq. jS
*
*-----------Find the DCR for A and B
*
            Call ICopy(nIrrep,iOper,1,iDCRR,1)
            nDCRR=nIrrep
            LmbdR=1
            If (iPrint.ge.49) Write (6,'(10A)')
     &         ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
*
*-----------Compute normalization factor due the DCR symmetrization
*           of the two basis functions and the operator.
*
            iuv = nStab(mdci)*nStab(mdcj)
            FactNd = DBLE(iuv)/DBLE(nIrrep * LmbdR)
            If (MolWgh.eq.1) Then
               FactNd = FactNd * DBLE(nIrrep)/DBLE(iuv)
            Else If (MolWgh.eq.2) Then
               FactNd = sqrt(DBLE(iuv))*DBLE(nIrrep)/DBLE(LmbdR)
            End If
*
*-----------Allocate memory for the elements of the Fock or 1st order
*           density matrix which are associated with the current shell
*           pair.
*
            Call mma_allocate(DSOp,nSO*iPrim*jPrim,label='DSOp')
            Call mma_allocate(DSOc,nSO*iBas *jBas,label='DSOc')
            Call mma_allocate(DSO,nSO*iPrim*jPrim,label='DSO')
            Call mma_allocate(Scrt,Max(iBas*jBas,iPrim*jPrim),
     &                        label='Scrt')
*
*-----------Introduce canonical order for the contracted basis
*
            If (iShell.ge.jShell) Then
              iSh   = iShell
              jSh   = jShell
              iAOi  = iAO
              jAOj  = jAO
              iBasi = iBas
              jBasj = jBas
              iPrimi= iPrim
              jPrimj= jPrim
              iCffi = iCff
              jCffj = jCff
              iAngi = iAng
              jAngj = jAng
              iCmpi = iCmp
              jCmpj = jCmp
              iShlli= iShll
              jShllj= jShll
            Else
              iSh   = jShell
              jSh   = iShell
              iAOi  = jAO
              jAOj  = iAO
              iBasi = jBas
              jBasj = iBas
              iPrimi= jPrim
              jPrimj= iPrim
              iCffi = jCff
              jCffj = iCff
              iAngi = jAng
              jAngj = iAng
              iCmpi = jCmp
              jCmpj = iCmp
              iShlli= jShll
              jShllj= iShll
            End If
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over the density and the spin-density (optional)
*
            Do iFD = 1,  mFD
               If (iPrint.ge.99) Then
                  If (iFD.eq.1) Then
                     Write (6,*) 'Processing the density'
                  Else
                     Write (6,*) 'Processing the spin-density'
                  End If
               End If
*                                                                      *
************************************************************************
*                                                                      *
*-----------Gather the elements from 1st order density / Fock matrix.
*
            Call SOGthr(DSOc,iBasi,jBasj,nSO,FD(1,iFD),
     &                  n2Tri(iSmLbl),iSmLbl,
     &                  iCmpi,jCmpj,iSh,jSh,AeqB,iAOi,jAOj)
*                                                                      *
************************************************************************
*                                                                      *
*-----------Project the Fock/1st order density matrix in AO
*           basis on to the primitive basis.
*
            If (iPrint.ge.99) Then
               Call RecPrt(' Left side contraction',' ',
     &                     DInf(iCffi),iPrimi,iBasi)
               Call RecPrt(' Right side contraction',' ',
     &                     DInf(jCffj),jPrimj,jBasj)
            End If
*
*-----------Transform IJ,AB to J,ABi
            Call DGEMM_('T','T',
     &                  jBasj*nSO,iPrimi,iBasi,
     &                  1.0d0,DSOc,iBasi,
     &                  DInf(iCffi),iPrimi,
     &                  0.0d0,DSOp,jBasj*nSO)
*-----------Transform J,ABi to AB,ij
            Call DGEMM_('T','T',
     &                  nSO*iPrimi,jPrimj,jBasj,
     &                  1.0d0,DSOp,jBasj,
     &                  DInf(jCffj),jPrimj,
     &                  0.0d0,DSO,nSO*iPrimi)
*-----------Transpose to ij,AB
            Call DGeTmO(DSO,nSO,nSO,iPrimi*jPrimj,DSOp,
     &                  iPrimi*jPrimj)
*
            If (iPrint.ge.99) Call
     &         RecPrt(' Decontracted 1st order density/Fock matrix',' ',
     &                        DSOp,iPrimi*jPrimj,nSO)
*                                                                      *
************************************************************************
*                                                                      *
*-----------Loops over symmetry operations.
*
            nOp(1) = NrOpr(0,iOper,nIrrep)
            If (iFD.eq.1) Then
*
*------------- Store away pointer to the block of density info
*
               ipOffD(1,ijShll) = jOffD + ipDeDe
               ipOffD(2,ijShll) = nDCRR
               If (nIrrep.eq.1.and.Special_NoSym) Then
                  ipOffD(3,ijShll) = iCmp*jCmp
     &                             + iPrim*jPrim + 1
               Else
                  ipOffD(3,ijShll) = iCmp*jCmp*(iBas*jBas+1)
     &                             + iPrim*jPrim + 1
               End If
            Else
*
*------------- Store away pointer to the block of spin-density info
*
               ipOffD(4,ijShll) = jOffD + ipDeDe
            End If
            mIndij = mIndij + (nIrrep/nStab(mdci))*
     &                        (nIrrep/nStab(mdcj))
            If (iPrint.ge.99) Then
               Write (6,*) ' ipDeDe+jOffD,nDCRR,iCmp*jCmp*iBas*jBas=',
     &                       ipDeDe+jOffD,nDCRR,iCmp*jCmp*iBas*jBas
            End If
            Do lDCRR = 0, nDCRR-1
               nOp(2) = NrOpr(iDCRR(lDCRR),iOper,nIrrep)
*                                                                      *
************************************************************************
*                                                                      *
*--------------Desymmetrize the 1st order density matrix(contracted).
*
               Call Desym1(iSmLbl,iAngi,jAngj,iCmpi,jCmpj,
     &                     iSh,jSh,iShlli,jShllj,
     &                     DAO,iBasi,jBasj,
     &                     DSOc,nSO,nOp,FactNd,Scrt)
*
*--------------Store away result
*
               Temp=Zero
               ipStart=ipDeDe+jOffD
               If (DFT_Storage) Then
*
*---------------- Storage format for numerical integration
*
*                 D(iBas*iCmp,jBas*jCmp)
*
                  Call ResortD(DAO,DeDe(ipStart),
     &                         iBas,iCmp,jBas,jCmp)
                  jOffD = jOffD + iBas*iCmp*jBas*jCmp
*
               Else
*
*---------------- Storage format for direct SCF
*
*                 D(iBas*jBas+1,iCmp*jCmp)
*
                  jpDAO = 1
                  Do ijCmp = 1, iCmp*jCmp
                     If (nIrrep.ne.1.or..Not.Special_NoSym) Then
                        call dcopy_(iBas*jBas,DAO(jpDAO),1,
     &                                       DeDe(ipDeDe+jOffD),1)
                        jOffD = jOffD + iBas*jBas
                     End If
*--------------------Find the largest density for this angular
*                    combination
                     iHigh = iDAMax_(iBas*jBas,DAO(jpDAO),1)
                     DeDe(ipDeDe+jOffD) = Abs(DAO(jpDAO+iHigh-1))
                     If (Temp.lt.Abs(DAO(jpDAO+iHigh-1))) Then
                        Temp=Abs(DAO(jpDAO+iHigh-1))
                     End If
                     jOffD = jOffD + 1
                     jpDAO = jpDAO + iBas*jBas
                  End Do
                  If ( (nIrrep.ne.1.or..Not.Special_NoSym) .and.
     &                iPrint.ge.99) Call RecPrt(' DAO(+AMax)',' ',
     &               DeDe(ipStart),iBas*jBas+1,iCmp*jCmp)
               End If
*
               If (DFT_Storage) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*--------------Desymmetrize the 1st order density matrix(primitive).
*
               Call Desym1(iSmLbl,iAngi,jAngj,iCmpi,jCmpj,
     &                     iSh,jSh,iShlli,jShllj,
     &                     DAO,iPrimi,jPrimj,
     &                     DSOp,nSO,nOp,FactNd,Scrt)
*
*--------------Change order so it follows what is used in TwoEl
*
               iFLOP1 = iBas*iPrim*jPrim + iBas*jPrim*jBas
               iFLOP2 = iPrim*jPrim*jBas + iBas*iPrim*jBas
               ipStart=ipDeDe+jOffD
               Do j = 1, jPrimj
                  Do i = 1, iPrimi
                     ij = (j-1)*iPrimi + i
                     iHigh = IDAMax_(iCmpi*jCmpj,DAO(ij),
     &                       iPrimi*jPrimj)-1
                     DeDe(ipDeDe+jOffD) =
     &                       Abs(DAO(ij+iHigh*iPrimi*jPrimj))
                     jOffD = jOffD + 1
                  End Do
               End Do
*--------------Find the overall largest density
               DeDe(ipDeDe+jOffD) = Temp
               jOffD = jOffD + 1
               If (iPrint.ge.99) Then
                  Call RecPrt(' D,prim',' ',
     &               DeDe(ipStart),iPrimi,jPrimj)
                  Write (6,*) ' Max(DAO)=',Temp
               End If
*
 99            Continue
*                                                                      *
************************************************************************
*                                                                      *
            End Do  ! lDCRR
*                                                                      *
************************************************************************
*                                                                      *
            End Do  ! iFD
*                                                                      *
************************************************************************
*                                                                      *
            Call mma_deallocate(Scrt)
            Call mma_deallocate(DSO)
            Call mma_deallocate(DSOc)
            Call mma_deallocate(DSOp)
            Call mma_deallocate(DAO)
 131        Continue
         End Do
      End Do
      mDeDe = jOffD
*
      If (mDeDe.ne.nDeDe) Then
         Write (6,*) 'DeDe:  mDeDe =', mDeDe,' nDeDe =', nDeDe
         Call ErrTra
         Call Abend
      End If
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(2,TCpu2-TCpu1,TWall2-TWall1)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nDInf)
      End
