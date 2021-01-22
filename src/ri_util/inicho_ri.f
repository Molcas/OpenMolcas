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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  IniCho_RI
*
*> @brief
*>   Initialize Cholesky environment for RI calculations
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Initialize Cholesky environment for RI calculations.
*>
*> @note
*> Needs a call to ::SetUp_Ints with indexation turned on
*>
*> @param[in] nSkal    The number of shells (excl. aux. basis)
*> @param[in] nVec_Aux Number of aux. basis vectors per irrep
*> @param[in] nIrrep   Number of irreps
*> @param[in] iTOffs   Offset vector
*> @param[in] iShij    Index vector of shell pairs
*> @param[in] nShij    Number of shell pairs
************************************************************************
      SubRoutine IniCho_RI(nSkal,nVec_Aux,nIrrep,iTOffs,iShij,nShij)
      Use Para_Info, Only: Is_Real_Par
      use ChoArr, only: iSP2F
      use ChoSwp, only: InfRed, InfVec
      Implicit None
      Integer nSkal, nIrrep, nShij
      Integer nVec_Aux(0:nIrrep-1)
      Integer iTOffs(3,nIrrep)
      Integer iShij(2,nShij)
#include "cholesky.fh"
#include "choprint.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#if defined (_MOLCAS_MPP_)
#include "choglob.fh"
#endif

      Logical SetDefaultsOnly, Skip_PreScreen, Alloc_Bkm
      Integer iDummy, LuOut, iSym, iVec, ijS
      Integer iTri, i, j

      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

C     Set defaults for those parameters that can normally be changed
C     through user input to the Cholesky decomposition.
C     -----------------------------------------------------------------

      SetDefaultsOnly = .True.
      iDummy = -1
      LuOut = 6
      Call Cho_Inp(SetDefaultsOnly,iDummy,LuOut)

C     Reset Cholesky Threshold for RI
C     -------------------------------
      Call Get_thrc_RI(ThrCom)

C     Reset parallel config.
C     ----------------------

      CHO_FAKE_PAR = .False.
      Call Cho_ParConf(CHO_FAKE_PAR)

C     Set run mode to "external" (should be irrelevant for RI).
C     ---------------------------------------------------------

      RUN_MODE = RUN_EXTERNAL

C     Silence the Cholesky routines.
C     ------------------------------

      iPrint = 0

C     Set number of shells (excl. aux. basis) in cholesky.fh
C     -------------------------------------------------------

      nShell = nSkal

C     To avoid unnecessary allocations of shell-pair-to-reduced-set
C     maps, set decomposition algorithm to 1 ("one-step") and the Seward
C     interface to "1" (full shell quadruple storage). Both values are,
C     obviously, irrelevant for RI. In parallel runs, use default values
C     to avoid warnings being printed in Cho_P_Check.
C     ------------------------------------------------------------------

      If (Is_Real_Par()) Then
         Cho_DecAlg = 4
         IfcSew = 2
      Else
         Cho_DecAlg = 1
         IfcSew = 1
      End If

C     Change MaxRed to 1 (all vectors have identical dimension, namely
C     full => only 1 reduced set).
C     ----------------------------------------------------------------

      MaxRed = 1

C     Set MaxVec to the largest number of vectors (= number of linearly
C     independent auxiliary basis functions). In this way we avoid
C     allocating more memory for InfVec than needed.
C     ------------------------------------------------------------------

      MaxVec = nVec_Aux(0)
      Do iSym = 1,nIrrep-1
         MaxVec = max(MaxVec,nVec_Aux(iSym))
      End Do

C     Other initializations. Most importantly, allocate InfRed and
C     InfVec arrays (defined in choswp.f90).
C     We skip diagonal prescreening, as it has already been done.
C     Instead, allocate and set the mapping from reduced to full shell
C     pairs here.
C     ----------------------------------------------------------------

      nnShl = nShij
      Call mma_allocate(iSP2F,nnShl,Label='iSP2F')
      Do ijS = 1,nnShl
         iSP2F(ijS) = iTri(iShij(1,ijS),iShij(2,ijS))
      End Do
      Skip_PreScreen = .True.
      Alloc_Bkm = .False.
      Call Cho_Init(Skip_PreScreen,Alloc_Bkm)

C     Set number of vectors equal to the number of lin. indep. auxiliary
C     basis functions.
C     ------------------------------------------------------------------

      Do iSym = 1,nSym
         NumCho(iSym) = nVec_Aux(iSym-1)
      End Do
#if defined (_MOLCAS_MPP_)
      If (Is_Real_Par()) Then
         Call iCopy(nSym,NumCho,1,NumCho_G,1)
         Call iZero(myNumCho,nSym)
      End If
#endif

C     Do allocations that are normally done during or after the
C     computation of the diagonal (since the dimension of the 1st
C     reduced set is unknown until the screened diagonal is known).
C     -------------------------------------------------------------

      Call IniCho_RI_Xtras(iTOffs,nIrrep,iShij,nShij)

C     Set start disk addresses.
C     -------------------------

      XnPass = 0 ! it should be zeroed in Cho_Inp, but just in case.
      Call Cho_SetAddr(InfRed,InfVec,MaxRed,MaxVec,SIZE(InfVec,2),nSym)

C     Set vector info.
C     Parent diagonal is set equal to the vector number, parent pass
C     (i.e. reduced set) to 1.
C     --------------------------------------------------------------

      Do iSym = 1,nSym
         Do iVec = 1,NumCho(iSym)
            Call Cho_SetVecInf(iVec,iSym,iVec,1,1)
         End Do
      End Do
*
      Return
      End
      SubRoutine IniCho_RI_Xtras(iTOffs,nIrrep,iShij,nShij)
      use ChoArr, only: iRS2F, nDimRS
      use ChoSwp, only: nnBstRSh, iiBstRSh
      use ChoSwp, only:   IndRSh,   IndRSh_Hidden
      use ChoSwp, only:   IndRed,   IndRed_Hidden
      Implicit None
      Integer nIrrep, nShij
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Logical DoDummy

      Integer iiBst(8), nnBst(8), iTOffs(3,nIrrep), iShij(2,nShij)
      Integer iSym, iCount, nnBstT

      Integer i

C     Define max. dimensions and offsets of the symmetry blocks of the
C     integrals matrix.
C     ----------------------------------------------------------------

      iCount = 0
      Do iSym = 1,nSym
         iiBst(iSym) = iCount
         nnBst(iSym) = iTOffs(3,iSym)
         iCount = iCount + nnBst(iSym)
      End Do

C     Set dimensions of reduced sets equal to full dimension.
C     -------------------------------------------------------

      Do i = 1,3
         nnBstT=0
         Do iSym = 1,nSym
            iiBstR(iSym,i) = iiBst(iSym)
            nnBstR(iSym,i) = nnBst(iSym)
            nnBstT=nnBstT+nnBstR(iSym,i)
         End Do
         nnBstRT(i) = nnBstT
      End Do
      mmBstRT = nnBstRT(1)

C     Allocate index arrays for reduced sets: IndRed and IndRsh.
C     ----------------------------------------------------------

      Call mma_allocate(IndRed_Hidden,nnBstRT(1),3,
     &                  Label='IndRed_Hidden')
      IndRed => IndRed_Hidden
      Call mma_allocate(IndRSh_Hidden,nnBstRT(1),Label='IndRSh_Hidden')
      IndRSh => IndRSh_Hidden

C     Allocate iScr array used by reading routines.
C     ---------------------------------------------

      DoDummy = .False.
      Call Cho_Allo_iScr(DoDummy)

C     Initialize reduced set dimensions used for reading vectors.
C     (Note: here they are all the same - there is one reduced sets!)
C     ---------------------------------------------------------------

      Do i = 1,MaxRed
         Do iSym = 1,nSym
            nDimRS(iSym,i) = nnBstR(iSym,1)
         End Do
      End Do

C     Allocate and set mapping array from 1st reduced set to full
C     storage.
C     -----------------------------------------------------------

      Call mma_allocate(iRS2F,2,nnBstRT(1),Label='iRS2F')

C     Set index arrays corresponding to full storage:
C     iiBstRSh, nnBstRSh, IndRed, IndRSh, and iRS2F.
C     -----------------------------------------------

      Call SetChoIndx_RI(iiBstRSh,nnBstRSh,
     &                   IndRed,IndRsh,iRS2F,
     &                   nSym,nnShl,nnBstRT(1),3,2,iShij,nShij)

      Return
      End
      SubRoutine SetChoIndx_RI(iiBstRSh,nnBstRSh,IndRed,IndRsh,iRS2F,
     &                         I_nSym,I_nnShl,I_mmBstRT,I_3,I_2,
     &                         iShij,nShij)
      use ChoArr, only: iSP2F, iBasSh, nBasSh, nBstSh
      Implicit Real*8 (a-h,o-z)
      Integer iiBstRSh(I_nSym,I_nnShl,I_3), nnBstRSh(I_nSym,I_nnShl,I_3)
      Integer IndRed(I_mmBstRT,I_3), IndRsh(I_mmBstRT)
      Integer iRS2F(I_2,I_mmBstRT), iShij(2,nShij)
#include "choorb.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"

      Integer  Cho_iSAOSh
      External Cho_iSAOSh

      Integer iRS(8)

      MulD2h(i,j)=iEor(i-1,j-1)+1
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

C     nnBstRSh(iSym,iSh_ij,1) = #elements in compound sym. iSym of
C                               shell-pair ab in 1st reduced set.
C     IndRSh(jRS): shell-pair to which element jRS of first reduced set
C                  belongs.
C     IndRed(jRS,1): address (without symmetry) in shell-pair of element
C                    jRS of first reduced set.
C     ------------------------------------------------------------------

      Call iCopy(nSym*nnShl,[0],0,nnBstRSh(1,1,1),1)
      Call iCopy(nSym,iiBstR(1,1),1,iRS,1)
      Do iSh_ij= 1,nShij
         iShla=iShij(1,iSh_ij)
         iShlb=iShij(2,iSh_ij)
         iShlab=iTri(iShla,iShlb)
C        Write (*,*) 'iSh_ij,iShlab,iShla,iShlb=',
C    &               iSh_ij,iShlab,iShla,iShlb
         If (iShlab .ne. iSP2F(iSh_ij)) Then
            Call SysAbendMsg('SetChoIndx_RI','SP2F setup error',' ')
         End If
*
         If (iShla.gt.iShlb) Then
*
*        code for shell a > shell b
*
            Do iSymb = 1,nSym
               Do ibb = 1,nBasSh(iSymb,iShlb)
                  ib = iBasSh(iSymb,iShlb) + ibb
                  Do iSyma = 1,nSym
                     iSym = MulD2h(iSyma,iSymb)
                     Do iaa = 1,nBasSh(iSyma,iShla)
                        ia = iBasSh(iSyma,iShla) + iaa
                        iab = nBstSh(iShla)*(ib-1) + ia
                        nnBstRSh(iSym,iSh_ij,1) =
     &                     nnBstRSh(iSym,iSh_ij,1) + 1
                        iRS(iSym) = iRS(iSym) + 1
                        IndRSh(iRS(iSym)) = iShlab
                        IndRed(iRS(iSym),1) = iab
                     End Do
                  End Do
               End Do
            End Do
*
         Else
*
*            code for shell a = shell b follows
*
            Do ia = 1,nBstSh(iShla)
               iSyma = Cho_iSAOSh(ia,iShla)
               Do ib = 1,ia
                  iab = iTri(ia,ib)
                  iSymb = Cho_iSAOSh(ib,iShlb)
                  iSym = MulD2h(iSyma,iSymb)
                  nnBstRSh(iSym,iSh_ij,1) = nnBstRSh(iSym,iSh_ij,1) + 1
                  iRS(iSym) = iRS(iSym) + 1
                  IndRSh(iRS(iSym)) = iShlab
                  IndRed(iRS(iSym),1) = iab
               End Do
            End Do
*
         End If
      End Do   ! iSh_ij

C     Check.
C     ------

      nErr = 0
      Do iSym = 1,nSym
         iCount = nnBstRSh(iSym,1,1)
         Do iSh_ij = 2,nnShl
            iCount = iCount + nnBstRSh(iSym,iSh_ij,1)
         End Do
         If (iCount .ne. nnBstR(iSym,1)) Then
            nErr = nErr + 1
         End If
      End Do
      If (nErr .ne. 0) Then
         Call SysAbendMsg('SetChoIndx_RI','Setup error',
     &                    'iCount vs. nnBstR')
      End If
      Do iSym = 1,nSym
         If ((iRS(iSym)-iiBstR(iSym,1)) .ne. nnBstR(iSym,1)) Then
            nErr = nErr + 1
         End If
      End Do
      If (nErr .ne. 0) Then
         Call SysAbendMsg('SetChoIndx_RI','Setup error','ShP RS1 count')
      End If

C     iiBstRSh(iSym,iSh_ij,1) = offset to elements in compound sym. iSym
C                               of shell-pair ab in 1st reduced set.
C     ------------------------------------------------------------------

      Do iSym = 1,nSym
         iiBstRSh(iSym,1,1) = 0
         Do iSh_ij = 2,nnShl
            iiBstRSh(iSym,iSh_ij,1) = iiBstRSh(iSym,iSh_ij-1,1)
     &                              + nnBstRSh(iSym,iSh_ij-1,1)
         End Do
      End Do

C     Check.
C     ------

      nErr = 0
      Do iSym = 1,nSym
         Do iSh_ij = 1,nnShl
            jRS1 = iiBstR(iSym,1) + iiBstRSh(iSym,iSh_ij,1) + 1
            jRS2 = jRS1 + nnBstRSh(iSym,iSh_ij,1) - 1
            Do jRS = jRS1,jRS2
               If (IndRSh(jRS) .ne. iSP2F(iSh_ij)) Then
                  nErr = nErr + 1
               End If
            End Do
         End Do
      End Do
      If (nErr .ne. 0) Then
         Call SysAbendMsg('SetChoIndx_RI','Setup error','IndRSh')
      End If

C     Copy index arrays to "locations" 2 and 3.
C     Note: IndRed here returns the index in 1st reduced set.
C     -------------------------------------------------------

      Do i = 2,3
         Do jRS = 1,nnBstRT(1)
            IndRed(jRS,i) = jRS
         End Do
         Call iCopy(nSym*nnShl,iiBstRSh(1,1,1),1,iiBstRSh(1,1,i),1)
         Call iCopy(nSym*nnShl,nnBstRSh(1,1,1),1,nnBstRSh(1,1,i),1)
      End Do

      Call Cho_RStoF(iRS2F,2,nnBstRT(1),1)

      Return
      End
************************************************************
*
************************************************************
      Subroutine Get_thrc_RI(Thr_CD)
      use RICD_Info, only: Thrshld_CD
      Real*8 Thr_CD

      Thr_CD = Thrshld_CD

      Return
      End
