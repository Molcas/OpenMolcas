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
* Copyright (C) 2008, Roland Lindh                                     *
************************************************************************
      Subroutine ReNorm2(iCnttp)
      use Wrj12
*
      Call ICopy(4*8,[0],0,iOffA,1)
      Do ire_do = 1, 2
*
         Call ReNorm2_Internal(iCnttp)
*
      End Do
*
      Return
      End
      Subroutine ReNorm2_Internal(iCnttp)
************************************************************************
*                                                                      *
*    Objective: Orthonormalize parts of the auxiliary basis set.       *
*                                                                      *
* Called from: Mk_RICD_Shells                                          *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theor. Chemi., Lund Univ., Sweden.*
*                                                                      *
*             Modified to transform the auxiliary basis to a true      *
*             Cholesky basis set while on TACC 2008 conference in      *
*             Songjiang District, Shanghai, China, 23-27 Sept. 2008.   *
*                                                                      *
************************************************************************
      use SOAO_Info, only: iAOtSO, nSOInf
      use Real_Spherical
      use Basis_Info
      use Sizes_of_Seward, only: S
      Implicit Real*8 (A-H,O-Z)
      External Integral_RI_2
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "print.fh"
#include "status.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8, Allocatable :: TInt_c(:), ADiag(:)
      Logical In_Core
      Interface
            Subroutine Drv2El_Atomic_NoSym(Integral_RI_2,
     &                                     ThrAO,iCnttp,jCnttp,
     &                                     TInt_c,nTInt_c,
     &                                     In_Core,ADiag,Lu_A,
     &                                     ijS_req,Keep_Shell)
            External Integral_RI_2
            Real*8 ThrAO
            Integer iCnttp, jCnttp, nTInt_c, Lu_A,ijS_req, Keep_Shell
            Real*8, Allocatable :: TInt_c(:), ADiag(:)
            Logical In_Core
            End Subroutine
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Let us now Gram-Schmidt orthonormalize the auxiliary basis for
*     better numerics and balance.
*
*     Update kOffAO and lOffAO to include the auxiliary basis too.
*
      Call Setup_OffAO()
*
*     Set up transformation matrix from Cartesian to real spherical
*     harmonics.
*
      Call Sphere(S%iAngMx)
*
      Call Flip_Flop(.False.) ! Contracted mode.
*
      Thr_CB=Max(1.0D-14,Thrshld_CD*1.0D-10)
      ThrAO=Zero
*
*     Do iCnttp = 1, nCnttp
*        Skip the dummy shell
         If (iCnttp.eq.iCnttp_dummy) Go To 2222
*        skip non-auxiliary basis sets
         If (.Not.dbsc(iCnttp)%Aux) Go To 2222
*                                                                      *
************************************************************************
*                                                                      *
*        Define some parameters to facilitate the atomic calculation
*
         S%nShlls= dbsc(iCnttp)%nVal
         nTest = dbsc(iCnttp)%nVal-1
*
*        Define AOtSO
*
         iAO = 0
         iSO = 0
         nSO=0
         Do iAng = 0, nTest
            iShll_ = dbsc(iCnttp)%iVal + iAng
            nCmp = (iAng+1)*(iAng+2)/2
            If (Shells(iShll_)%Prjct) nCmp = 2*iAng+1
            iSO = 0
            If (Shells(iShll_)%nBasis_C*Shells(iShll_)%nExp==0) Cycle
            Do iCmp = 1, nCmp
               iAO = iAO + 1
               If (iAO>nSOInf) Then
                  Write (6,*) 'renorm2_internal: iAO>nSOInf'
                  Write (6,*) 'iAO=',iAO
                  Write (6,*) 'nSOInf=',nSOInf
                  Call Abend()
               End If
               iAOtSO(iAO,0) = iSO + 1
               iSO = iSO + Shells(iShll_)%nBasis
            End Do
            nSO=nSO+iSO
         End Do
*
         ijS_req=0
         Keep_Shell=iTabMx
         Do iAng = 0, nTest
            iShll = dbsc(iCnttp)%iVal + iAng
            nExpi = Shells(iShll)%nExp
            nBasisi=Shells(iShll)%nBasis
            If (nExpi*nBasisi.eq.0) Go To 2221
*
            nCmp = (iAng+1)*(iAng+2)/2
            If (Shells(iShll)%Prjct) nCmp = 2*iAng+1
*
            ijS_req=ijS_req+1
*
            Call Drv2El_Atomic_NoSym(Integral_RI_2,
     &                               ThrAO,iCnttp,iCnttp,
     &                               TInt_c,nTInt_c,
     &                               In_Core,ADiag,Lu_A,ijS_req,
     &                               Keep_Shell)
#ifdef _DEBUG_
            Call TriPrt('TInt_c',' ',TInt_c,nTInt_c)
#endif
*
            If (.NOT.In_Core) Then
               Call WarningMessage(2,'Error in ReNorm')
               Write (6,*) 'Out-of-core acCD not implemented!'
               Call Abend()
            End If
*
*           Produce the reduced set, in-place reduction.
*
            Call Allocate_Work(ipA,nTInt_c**2)
            ijT=0
            Do iBas = 1, nTInt_c
               Do jBas = 1, iBas
                  ijT=ijT+1
                  ijS=(jBas-1)*nTInt_c+iBas
                  jiS=(iBas-1)*nTInt_c+jBas
                  Work(ipA+ijS-1)=TInt_c(ijT)
                  Work(ipA+jiS-1)=TInt_c(ijT)
               End Do
            End Do
            Call mma_deallocate(TInt_c)
            ip_TInt_c=ipA
#ifdef _DEBUG_
            Call RecPrt('TInt_c',' ',Work(ip_TInt_c),nTInt_c,nTInt_c)
#endif
*
            ij=0
            iCmp=1
            jCmp=1
            Do jBas = 1, nBasisi
               j=(jCmp-1)*nBasisi+jBas
               Do iBas = 1, nBasisi
                  i=(iCmp-1)*nBasisi+iBas
                  ijF=(j-1)*nBasisi*nCmp+i
                  ij=ij+1
*
                  Work(ij+ip_TInt_c-1)=Work(ijF+ip_TInt_c-1)
*
               End Do
            End Do
#ifdef _DEBUG_
            Call RecPrt('TInt_c(r)','(5G20.10)',
     &                  Work(ip_TInt_c),nBasisi,nBasisi)
#endif
*
            Call GetMem('ADiag','Allo','Real',ipADiag,nBasisi)
            Call GetMem('iADiag','Allo','Inte',ipiADiag,nBasisi)
*
            iSeed=77
            Lu_A=IsFreeUnit(iSeed)
            Call DaName_MF_WA(Lu_A,'AMat09')
*
            iDisk=0
            Call dDaFile(Lu_A,1,Work(ip_TInt_c),nBasisi**2,iDisk)
*
            iSeed=iSeed+1
            Lu_Q=IsFreeUnit(iSeed)
            Call DaName_MF_WA(Lu_Q,'QMat09')
*
            call dcopy_(nBasisi,Work(ip_TInt_c),nBasisi+1,
     &                 Work(ipADiag),1)
*
            Call CD_AInv_(nBasisi,m,Work(ipADiag),iWork(ipiADiag),
     &                    Lu_A,Lu_Q,Thr_CB)
*
            Call GetMem('iADiag','Free','Inte',ipiADiag,n)
            Call GetMem('ADiag','Free','Real',ipADiag,n)
            Call Free_Work(ip_TInt_c)
*
*           Transform the contraction coefficients according to the
*           Cholesky vectors.
*
            Call Allocate_Work(ipTmp,nBasisi*nExpi)
            Call Allocate_Work(ipQVec,nBasisi**2)
            Call FZero(Work(ipQVec),nBasisi**2)
*
            iDisk=0
            Call dDaFile(Lu_Q,2,Work(ipQVec),nBasisi*m,iDisk)
            Call DaEras(Lu_Q)
#ifdef _DEBUG_
            Call RecPrt('QVec',' ',Work(ipQVec),nBasisi,m)
#endif
*
            iOff=0
            Do iCase = 1, 2
               call dcopy_(nExpi*nBasisi,
     &                     Shells(iShll)%Cff_c(1,1,iCase),1,
     &                     Work(ipTmp),1)
#ifdef _DEBUG_
               Call RecPrt('Coeff(old)',' ',
     &                     Shells(iShll)%Cff_c(1,1,iCase),
     &                     nExpi,nBasisi)
#endif
               Call DGEMM_('N','N',
     &                    nExpi,nBasisi,nBasisi,
     &                    1.0D0,Work(ipTmp),nExpi,
     &                          Work(ipQVec),nBasisi,
     &                    0.0D0,Shells(iShll)%Cff_c(1,1,iCase),
     &                          nExpi)
#ifdef _DEBUG_
               Call RecPrt('Coeff(new)',' ',
     &                     Shells(iShll)%Cff_c(1,1,iCase),
     &                     nExpi,nBasisi)
#endif
            End Do
*
            Call Free_Work(ipQVec)
            Call Free_Work(ipTmp)
*
 2221       Continue
         End Do
*
 2222    Continue
*     End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

