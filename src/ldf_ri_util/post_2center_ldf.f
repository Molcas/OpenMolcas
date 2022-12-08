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
* Copyright (C) 1990,1991,1993,1998,2005, Roland Lindh                 *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Post_2Center_LDF(A_Diag,AB,MaxCntr,Lu_AB,Local_A,
     &                            SO2C,nSO_Aux)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals.                          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March 1990                                               *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August 1991                        *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. 1993         *
*             Modified driver. Jan. 1998                               *
*             Modified to 2-center ERIs for RI June 2005               *
************************************************************************
      use Basis_Info, only: nBas_Aux
      use RI_glob, only: iOffA, Lu_A, Lu_Q, nChV
      use Gateway_global, only: force_out_of_core
      use RICD_Info, only: Thrshld_CD
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "nsd.fh"
      Character Name_Q*6
      Integer nQvec(0:7)
      Real*8 :: A_Diag(*)
      Real*8, Allocatable :: Local_A(:,:)
      Integer, Allocatable :: SO2C(:), AB(:,:)

      Real*8, Allocatable :: Scr(:)
      Integer, Allocatable :: iDiag(:), SO2lO(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the number of AB blocks
*
      nAB=MaxCntr*(MaxCntr+1)/2
      Call mma_allocate(AB,2,nAB,Label='AB')
      AB(:,:)=0
*
*     Find the size of the largest AB block!
*
      Max_AA=0
      Do iCenter = 1, MaxCntr
         nCenter=0
         Do i = 1, nSO_Aux
            jCenter=SO2C(i)
            If (jCenter.eq.iCenter) nCenter=nCenter+1
         End Do
         Max_AA = Max(Max_AA,nCenter)
      End Do
      nLocal_A=(2*Max_AA)**2
      Call mma_allocate(Local_A,nLocal_A,2,Label='Local_A')
      Call mma_allocate(SO2lO,nSO_Aux,Label='SO2lO')
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Some setup for the generation of the Q-vectors
*
      nB=0
      nXZ=0
      nQm_Tot=0
      lScr=0
      Do iIrrep = 0, nIrrep-1
         iOffA(3,iIrrep)=nB
         mB=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) mB = mB - 1
         nB=nB+mB
         nQm_Tot=nQm_Tot + mB**2
         nXZ = Max(nXZ,mB)
      End Do
*
      nA_Diag=nB
      lScr=3*nXZ
      Call mma_allocate(iDiag,nA_Diag,Label='iDiag')
      Call mma_allocate(Scr,lScr,Label='Scr')
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*     Process the atomic and diatomic blocks of the A matrix
*
      iAdr_AB=0
      Lu_AB=22
      Call DaName_MF_WA(Lu_AB,'ABLOCKS')
      Do iCenter = 1, MaxCntr
         niSO=0
         nlO=0
         Do iSO = 1, nSO_Aux
            kCenter = SO2C(iSO)
            If (kCenter.eq.iCenter) Then
               niSO=niSO+1
               nlO = nlO + 1
               SO2lO(iSO)=nlO
            End If
         End Do

         Do jCenter = 1, iCenter
            ijCenter = iCenter*(iCenter-1)/2 + jCenter
*
*           Process the atomic blocks
*
            iAddr=0
            iIrrep=0
            If (jCenter.eq.iCenter) Then
*
               Local_A(1:niSO**2,1)=Zero
*
               Do iSO = 1, nSO_Aux
                  Call dDaFile(Lu_A(iIrrep),2,Scr,
     &                         nSO_Aux,iAddr)
C                 Call RecPrt('Scr','(6G23.15)',Scr,1,nSO_Aux)
                  kCenter = SO2C(iSO)
                  If (kCenter.eq.iCenter) Then
                      ilO=SO2lO(iSO)
                      Do jSO = 1, iSO
                         lCenter = SO2C(jSO)
                         If (lCenter.eq.iCenter) Then
                            AElement=Scr(jSO)
                            jlO=SO2lO(jSO)
                            ij = (ilO-1)*nlO + jlO
                            ji = (jlO-1)*nlO + ilO
                            AElement=Scr(jSO)
                            Local_A(ij,1)=AElement
                            Local_A(ji,1)=AElement
                         End If
                      End Do
                  End If
               End Do
C              Call RecPrt('Local_A','(6G23.15)',
C    &                     Local_A(:,1),nlO,nlO)
*
               Call CD_AInv(Local_A(:,1),nlO,Local_A(:,2),Thrshld_CD)
*
C              Call RecPrt('Local_AInv','(6G23.15)',
C    &                     Local_A(:,2),nlO,nlO)
*
*              Store the address of the matrix, size and write the
*              matrix to *              disk.
*
               AB(1,ijCenter)=iAdr_AB
               AB(2,ijCenter)=nlO
               Call dDaFile(Lu_AB,1,Local_A(:,2),nlO**2,iAdr_AB)
*
            Else
*
*              Process the diatomic blocks
*
               njSO=0
               mlO = 0
               Do jSO = 1, nSO_Aux
                  lCenter = SO2C(jSO)
                  If (lCenter.eq.jCenter) Then
                     njSO=njSO+1
                     mlO = mlO + 1
                     SO2lO(jSO)=mlO+nlO
                  End If
               End Do
*
               Local_A(1:(niSO+njSO)**2,1)=Zero
*
               Do iSO = 1, nSO_Aux
                  Call dDaFile(Lu_A(iIrrep),2,Scr,nSO_Aux,iAddr)
                  kCenter = SO2C(iSO)
                  If (kCenter.eq.iCenter.or.
     &                kCenter.eq.jCenter) Then
                      ilO=SO2lO(iSO)
                      Do jSO = 1, iSO
                         lCenter = SO2C(jSO)
                         If (lCenter.eq.iCenter.or.
     &                       lCenter.eq.jCenter) Then
                            AElement=Scr(jSO)
                            jlO=SO2lO(jSO)
                            ij = (ilO-1)*(nlO+mlO) + jlO
                            ji = (jlO-1)*(nlO+mlO) + ilO
                            AElement=Scr(jSO)
                            Local_A(ij,1)=AElement
                            Local_A(ji,1)=AElement
                         End If
                      End Do
                  End If
               End Do
*
               Call CD_AInv(Local_A(:,1),nlO+mlO,Local_A(:,2),
     &                      Thrshld_CD)
*
*              Store the address of the matrix, size and write the
*               matrix to *              disk.
*
               AB(1,ijCenter)=iAdr_AB
               AB(2,ijCenter)=nlO+mlO
               Call dDaFile(Lu_AB,1,Local_A(:,2),(nlO+mlO)**2,iAdr_AB)
*
            End If
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Fill in the lower part of the A matrix as it is stored on disk.
*
      Call mma_maxDBLE(MaxMem)
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB = nB - 1 ! subtract dummay af
         Call Square_A(Lu_A(iIrrep),nB,MaxMem,Force_Out_of_Core)
      End Do
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      ThrQ=1.0D-12
      ichk=0
      is=1
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB = nB - 1 ! subtract dummy aux. func.
         If (nB.gt.0) Then
*
            iSeed=55+iIrrep
            Lu_Q(iIrrep)=IsFreeUnit(iSeed)
            Write(Name_Q,'(A4,I2.2)') 'QVec',iIrrep
            Call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)
*
            Call get_pivot_idx(A_Diag(is),nB,nQvec(iIrrep),
     &                         Lu_A(iIrrep),Lu_Q(iIrrep),
     &                         iDiag(is),Scr,lScr,ThrQ)
            nChV(iIrrep)=nQvec(iIrrep)
            ichk=ichk+Min(1,nB-nQvec(iIrrep))
            is=is+nB
         End If
         Call DaClos(Lu_A(iIrrep))
      End Do
*
      If (ichk.ne.0) Then
         Write(6,*)
         Write(6,*)'Detected lin. dependences in the auxiliary basis.'
         Write(6,'(A,8I6)')
     & ' # of AuxBas before l. d. removal: ',
     &   nBas_Aux(0)-1,(nBas_Aux(i),i=1,nIrrep-1)
         Write(6,'(A,8I6)')
     & ' # of AuxBas after  l. d. removal: ',(nQvec(i),i=0,nIrrep-1)
         Write(6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Scr)
      Call mma_deallocate(iDiag)
      Call mma_deallocate(SO2lO)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
