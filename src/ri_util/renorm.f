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
      Subroutine ReNorm()
      use Wrj12
*
      Call ICopy(4*8,[0],0,iOffA,1)
      Do ire_do = 1, 2
*
         Call ReNorm_()
*
      End Do
*
      Return
      End
      Subroutine ReNorm_()
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
      use Real_Spherical
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
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('ReNorm')
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
      If (Allocated(RSph)) Call mma_deallocate(RSph)
      If (Allocated(ipSph)) Call mma_deallocate(ipSph)
      Call Sphere(iAngMx)
*
      Call Flip_Flop(.False.) ! Contracted mode.
*
      Thr_CB=Max(1.0D-14,Thrshld_CD*1.0D-10)
      ThrAO=Zero
*
      Do iCnttp = 1, nCnttp
*        Skip the dummy shell
         If (iCnttp.eq.iCnttp_dummy) Go To 2222
*        skip non-auxiliary basis sets
         If (.Not.AuxCnttp(iCnttp)) Go To 2222
*                                                                      *
************************************************************************
*                                                                      *
*        Define some parameters to facilitate the atomic calculation
*
         nShlls= nVal_Shells(iCnttp)
         nTest= nVal_Shells(iCnttp)-1
*
*        Define AOtSO
*
         iAO = 0
         iSO = 0
         nSO=0
         Do iAng = 0, nTest
            iShll_ = ipVal(iCnttp) + iAng
            nCmp = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll_)) nCmp = 2*iAng+1
            iSO = 0
            If (nBasis_Cntrct(iShll_).ne.0 .and.
     &          nExp(iShll_).ne.0) Then
               Do iCmp = 1, nCmp
                  iAO = iAO + 1
                  iAOtSO(iAO,0) = iSO + 1
                  Do iCont = 1, nBasis(iShll_)
                     iSO = iSO + 1
                  End Do
               End Do
            End If
            nSO=nSO+iSO
         End Do
*
         ijS_req=0
         Keep_Shell=iTabMx
         Do iAng = 0, nTest
            iShll = ipVal(iCnttp) + iAng
            If (nExp(iShll)*nBasis(iShll).eq.0) Go To 2221
*
            nCmp = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll)) nCmp = 2*iAng+1
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
            Do jBas = 1, nBasis(iShll)
               j=(jCmp-1)*nBasis(iShll)+jBas
               Do iBas = 1, nBasis(iShll)
                  i=(iCmp-1)*nBasis(iShll)+iBas
                  ijF=(j-1)*nBasis(iShll)*nCmp+i
                  ij=ij+1
*
                  Work(ij+ip_TInt_c-1)=Work(ijF+ip_TInt_c-1)
*
               End Do
            End Do
#ifdef _DEBUG_
            Call RecPrt('TInt_c(r)','(5G20.10)',
     &                  Work(ip_TInt_c),nBasis(iShll),
     &                  nBasis(iShll))
#endif
*
            Call GetMem('ADiag','Allo','Real',ipADiag,nBasis(iShll))
            Call GetMem('iADiag','Allo','Inte',ipiADiag,nBasis(iShll))
*
            iSeed=77
            Lu_A=IsFreeUnit(iSeed)
            Call DaName_MF_WA(Lu_A,'AMat09')
*
            iDisk=0
            Call dDaFile(Lu_A,1,Work(ip_TInt_c),nBasis(iShll)**2,iDisk)
*
            iSeed=iSeed+1
            Lu_Q=IsFreeUnit(iSeed)
            Call DaName_MF_WA(Lu_Q,'QMat09')
*
            call dcopy_(nBasis(iShll),Work(ip_TInt_c),nBasis(iShll)+1,
     &                 Work(ipADiag),1)
*
            Call CD_AInv_(nBasis(iShll),m,Work(ipADiag),iWork(ipiADiag),
     &                    Lu_A,Lu_Q,Thr_CB)
*
            Call GetMem('iADiag','Free','Inte',ipiADiag,n)
            Call GetMem('ADiag','Free','Real',ipADiag,n)
            Call Free_Work(ip_TInt_c)
*
*           Transform the contraction coefficients according to the
*           Cholesky vectors.
*
            Call Allocate_Work(ipTmp,nBasis(iShll)*nExp(iShll))
            Call Allocate_Work(ipQVec,nBasis(iShll)**2)
            Call FZero(Work(ipQVec),nBasis(iShll)**2)
*
            iDisk=0
            Call dDaFile(Lu_Q,2,Work(ipQVec),nBasis(iShll)*m,iDisk)
            Call DaEras(Lu_Q)
#ifdef _DEBUG_
            Call RecPrt('QVec',' ',Work(ipQVec),nBasis(iShll),m)
#endif
*
            iOff=0
            Do iCase = 1, 2
               call dcopy_(nExp(iShll)*nBasis(iShll),
     &                    Work(ipCff(iShll)+iOff),1,
     &                    Work(ipTmp),1)
#ifdef _DEBUG_
               Call RecPrt('Coeff(old)',' ',Work(ipCff(iShll)+ioff),
     &                     nExp(iShll),nBasis(iShll))
#endif
               Call DGEMM_('N','N',
     &                    nExp(iShll),nBasis(iShll),nBasis(iShll),
     &                    1.0D0,Work(ipTmp),nExp(iShll),
     &                          Work(ipQVec),nBasis(iShll),
     &                    0.0D0,Work(ipCff(iShll)+iOff),nExp(iShll))
#ifdef _DEBUG_
               Call RecPrt('Coeff(new)',' ',Work(ipCff(iShll)+iOff),
     &                     nExp(iShll),nBasis(iShll))
#endif
               iOff = iOff + nBasis(iShll)*nExp(iShll)
            End Do
*
            Call Free_Work(ipQVec)
            Call Free_Work(ipTmp)
*
 2221       Continue
         End Do
*
 2222    Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('ReNorm')
      Return
      End
      Subroutine ReNorm2(iCnttp)
      use Wrj12
*
      Call ICopy(4*8,[0],0,iOffA,1)
      Do ire_do = 1, 2
*
         Call ReNorm2_(iCnttp)
*
      End Do
*
      Return
      End
      Subroutine ReNorm2_(iCnttp)
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
      use Real_Spherical
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
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('ReNorm')
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
      If (Allocated(RSph)) Call mma_deallocate(RSph)
      If (Allocated(ipSph)) Call mma_deallocate(ipSph)
      Call Sphere(iAngMx)
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
         If (.Not.AuxCnttp(iCnttp)) Go To 2222
*                                                                      *
************************************************************************
*                                                                      *
*        Define some parameters to facilitate the atomic calculation
*
         nShlls= nVal_Shells(iCnttp)
         nTest= nVal_Shells(iCnttp)-1
*
*        Define AOtSO
*
         iAO = 0
         iSO = 0
         nSO=0
         Do iAng = 0, nTest
            iShll_ = ipVal(iCnttp) + iAng
            nCmp = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll_)) nCmp = 2*iAng+1
            iSO = 0
            If (nBasis_Cntrct(iShll_).ne.0 .and.
     &          nExp(iShll_).ne.0) Then
               Do iCmp = 1, nCmp
                  iAO = iAO + 1
                  iAOtSO(iAO,0) = iSO + 1
                  Do iCont = 1, nBasis(iShll_)
                     iSO = iSO + 1
                  End Do
               End Do
            End If
            nSO=nSO+iSO
         End Do
*
         ijS_req=0
         Keep_Shell=iTabMx
         Do iAng = 0, nTest
            iShll = ipVal(iCnttp) + iAng
            If (nExp(iShll)*nBasis(iShll).eq.0) Go To 2221
*
            nCmp = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll)) nCmp = 2*iAng+1
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
            Do jBas = 1, nBasis(iShll)
               j=(jCmp-1)*nBasis(iShll)+jBas
               Do iBas = 1, nBasis(iShll)
                  i=(iCmp-1)*nBasis(iShll)+iBas
                  ijF=(j-1)*nBasis(iShll)*nCmp+i
                  ij=ij+1
*
                  Work(ij+ip_TInt_c-1)=Work(ijF+ip_TInt_c-1)
*
               End Do
            End Do
#ifdef _DEBUG_
            Call RecPrt('TInt_c(r)','(5G20.10)',
     &                  Work(ip_TInt_c),nBasis(iShll),
     &                  nBasis(iShll))
#endif
*
            Call GetMem('ADiag','Allo','Real',ipADiag,nBasis(iShll))
            Call GetMem('iADiag','Allo','Inte',ipiADiag,nBasis(iShll))
*
            iSeed=77
            Lu_A=IsFreeUnit(iSeed)
            Call DaName_MF_WA(Lu_A,'AMat09')
*
            iDisk=0
            Call dDaFile(Lu_A,1,Work(ip_TInt_c),nBasis(iShll)**2,iDisk)
*
            iSeed=iSeed+1
            Lu_Q=IsFreeUnit(iSeed)
            Call DaName_MF_WA(Lu_Q,'QMat09')
*
            call dcopy_(nBasis(iShll),Work(ip_TInt_c),nBasis(iShll)+1,
     &                 Work(ipADiag),1)
*
            Call CD_AInv_(nBasis(iShll),m,Work(ipADiag),iWork(ipiADiag),
     &                    Lu_A,Lu_Q,Thr_CB)
*
            Call GetMem('iADiag','Free','Inte',ipiADiag,n)
            Call GetMem('ADiag','Free','Real',ipADiag,n)
            Call Free_Work(ip_TInt_c)
*
*           Transform the contraction coefficients according to the
*           Cholesky vectors.
*
            Call Allocate_Work(ipTmp,nBasis(iShll)*nExp(iShll))
            Call Allocate_Work(ipQVec,nBasis(iShll)**2)
            Call FZero(Work(ipQVec),nBasis(iShll)**2)
*
            iDisk=0
            Call dDaFile(Lu_Q,2,Work(ipQVec),nBasis(iShll)*m,iDisk)
            Call DaEras(Lu_Q)
#ifdef _DEBUG_
            Call RecPrt('QVec',' ',Work(ipQVec),nBasis(iShll),m)
#endif
*
            iOff=0
            Do iCase = 1, 2
               call dcopy_(nExp(iShll)*nBasis(iShll),
     &                    Work(ipCff(iShll)+iOff),1,
     &                    Work(ipTmp),1)
#ifdef _DEBUG_
               Call RecPrt('Coeff(old)',' ',Work(ipCff(iShll)+ioff),
     &                     nExp(iShll),nBasis(iShll))
#endif
               Call DGEMM_('N','N',
     &                    nExp(iShll),nBasis(iShll),nBasis(iShll),
     &                    1.0D0,Work(ipTmp),nExp(iShll),
     &                          Work(ipQVec),nBasis(iShll),
     &                    0.0D0,Work(ipCff(iShll)+iOff),nExp(iShll))
#ifdef _DEBUG_
               Call RecPrt('Coeff(new)',' ',Work(ipCff(iShll)+iOff),
     &                     nExp(iShll),nBasis(iShll))
#endif
               iOff = iOff + nBasis(iShll)*nExp(iShll)
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
      Call qExit('ReNorm')
      Return
      End

