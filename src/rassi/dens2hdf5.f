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
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
      Module Dens2HDF5

      Integer, Allocatable :: IdxState(:,:)
      Real*8, Parameter, Private :: Thrs=1.0D-14

      Contains

************************************************************************
*  UpdateIdx
*
*> @brief
*>   Update index of TDMs to save to HDF5
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Update the table of indices of states for the TDMs that will be saved
*> on the HDF5 file. If SUBSets keyword is not used, all TDMs are saved,
*> otherwise only those for the selected transitions. Since the TDMs are
*> saved in SF basis, for a SO calculation we need to figure out which
*> SF TDMs contribute to the desired transitions.
*>
*> @param[in] IndexE  SF states sorted by energy
*> @param[in] nSS     number of SO states
*> @param[in] nSS     number of SO states
*> @param[in] USOR    SO coefficients in SF basis (real part)
*> @param[in] USOI    SO coefficients in SF basis (imaginary part)
*> @param[in] MapSt   map of SF states expanded by multiplicity
************************************************************************
      Subroutine UpdateIdx(IndexE, nSS, USOR, USOI, MapSt)
      Implicit None
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer, Dimension(nState), Intent(In) :: IndexE
      Integer, Intent(In) :: nSS
      Real*8, Intent(In), Optional :: USOR(nSS,nSS),USOI(nSS,nSS)
      Integer, Intent(In), Optional :: MapSt(nSS)
      Integer :: i,j,i_,j_,iState,jState,iSS,jSS,iEnd,jStart
      Real*8 :: f1,f2
      If (.not.Allocated(IdxState)) Then
        Call mma_Allocate(IdxState,nState,nState,Label='IdxState')
        Call iCopy(nState**2,[0],0,IdxState,1)
      End If
      If (ReduceLoop) Then
        iEnd=LoopDivide
        jStart=LoopDivide+1
      Else
        iEnd=nState
        If (nSS.gt.0) iEnd=nSS
        jStart=1
      End If
!     list states for/between which we want to store the density matrices
      If (nSS.gt.0) Then
        Do iSS=1,nSS
          Do jSS=1,nSS
            If ((jSS.ne.iSS).and.
     &          ((jSS.lt.jStart).or.(iSS.gt.iEnd))) Cycle
            Do i=1,nSS
              iState=MapSt(i)
              Do j=1,i
                jState=MapSt(j)
                i_=Max(jState,iState)
                j_=Min(jState,iState)
                If (IdxState(i_,j_).gt.0) Cycle
                f1=(USOR(i,iSS)*USOR(j,jSS)+USOI(i,iSS)*USOI(j,jSS))**2
     &            +(USOR(i,iSS)*USOI(j,jSS)-USOI(i,iSS)*USOR(j,jSS))**2
                f2=(USOR(j,iSS)*USOR(i,jSS)+USOI(j,iSS)*USOI(i,jSS))**2
     &            +(USOR(j,iSS)*USOI(i,jSS)-USOI(j,iSS)*USOR(i,jSS))**2
!               this should be Thrs**2, but let's be looser with SO states
                If (Max(f1,f2).ge.Thrs) IdxState(i_,j_)=2
            End Do
            End Do
          End Do
        End Do
      Else
        Do i=1,nState
          iState=IndexE(i)
          IdxState(iState,iState)=1
          If (i.gt.iEnd) Cycle
          Do j=Max(i+1,jStart),nState
            jState=IndexE(j)
            i_=Max(jState,iState)
            j_=Min(jState,iState)
            IdxState(i_,j_)=1
          End Do
        End Do
      End If
      End Subroutine UpdateIdx

************************************************************************
*  StoreDens
*
*> @brief
*>   Store the density matrices to HDF5
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Save the selected state and transition density matrices to the HDF5
*> file. The TDMs were computed in input state basis, but will be stored
*> in SF eigen states, so a transformation is needed.
*>
*> @param[in] EigVec  coefficients of the SF eigen states in the input state basis
************************************************************************
#ifdef _HDF5_
      Subroutine StoreDens(EigVec)
      Use rassi_aux, Only : iDisk_TDM
      Use rassi_global_arrays, Only : JbNum
      Use mh5, Only: mh5_put_dset_array_real
      Implicit None
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "rassiwfn.fh"
      Real*8, Intent(In) :: EigVec(nState,nState)
      Integer :: iState,jState,k,l,nThisTDMZZ
      Integer :: Job1,Job2,iSym1,iSym2,iSy12,iDisk,iEmpty,iOpt,iGo
      Logical isZero(3)
      Real*8, Allocatable :: TDMZZ(:),TSDMZZ(:),WDMZZ(:),
     &                       TDMIJ(:),TSDMIJ(:),WDMIJ(:)
      Real*8 :: f1,f2
      If (.not.Allocated(IdxState)) Return
C Transform TDMs to SF eigenstates
      Call mma_Allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
      Call mma_Allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
      Call mma_Allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
      Call mma_Allocate(TDMIJ,nTDMZZ,Label='TDMIJ')
      Call mma_Allocate(TSDMIJ,nTDMZZ,Label='TSDMIJ')
      Call mma_Allocate(WDMIJ,nTDMZZ,Label='WDMIJ')
      Do iState=1,nState
        Do jState=1,iState
          If (IdxState(iState,jState).eq.0) Cycle
          Call dCopy_(nTDMZZ,[0.0D0],0,TDMIJ,1)
          Call dCopy_(nTDMZZ,[0.0D0],0,TSDMIJ,1)
          Call dCopy_(nTDMZZ,[0.0D0],0,WDMIJ,1)
          isZero=[.True.,.True.,.True.]
          Do k=1,nState
            Job1=JbNum(k)
            iSym1=Irrep(Job1)
            Do l=1,k
              f1=EigVec(k,iState)*EigVec(l,jState)
              f2=EigVec(l,iState)*EigVec(k,jState)
              If (Max(Abs(f1),Abs(f2)).lt.Thrs) Cycle
              Job2=JbNum(l)
              iSym2=Irrep(Job2)
              iSy12=Mul(iSym1,iSym2)
              iDisk=iDisk_TDM(k,l,1)
              iEmpty=iDisk_TDM(k,l,2)
              iOpt=2
              iGo=3
              If (IfSO) iGo=iGo+4
              Call dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,
     &                       LuTDM,iDisk,iEmpty,iOpt,iGo,k,l)
              If (IAnd(iEmpty,1).ne.0) Then
                isZero(1)=.False.
                If (Abs(f1).ge.Thrs)
     &            Call daXpY_(nTDMZZ,f1,TDMZZ,1,TDMIJ,1)
                If ((k.ne.l).and.(Abs(f2).ge.Thrs)) Then
                  Call Transpose_TDM(TDMZZ,iSy12)
                  Call daXpY_(nTDMZZ,f2,TDMZZ,1,TDMIJ,1)
                End If
              End If
              If (IAnd(iEmpty,2).ne.0) Then
                isZero(2)=.False.
                If (Abs(f1).ge.Thrs)
     &            Call daXpY_(nTDMZZ,f1,TSDMZZ,1,TSDMIJ,1)
                If ((k.ne.l).and.(Abs(f2).ge.Thrs)) Then
                  Call Transpose_TDM(TSDMZZ,iSy12)
                  Call daXpY_(nTDMZZ,f2,TSDMZZ,1,TSDMIJ,1)
                End If
              End If
              If (IFSO.and.(IAnd(iEmpty,4).ne.0)) Then
                isZero(3)=.False.
                If (Abs(f1).ge.Thrs)
     &            Call daXpY_(nTDMZZ,f1,WDMZZ,1,WDMIJ,1)
                If ((k.ne.l).and.(Abs(f2).ge.Thrs)) Then
                  Call Transpose_TDM(WDMZZ,iSy12)
                  Call daXpY_(nTDMZZ,f2,WDMZZ,1,WDMIJ,1)
                End If
              End If
            End Do
          End Do
          If (All(isZero)) Cycle
          nThisTDMZZ=0
          Do iSym1=1,nSym
            iSym2=Mul(iSy12,iSym1)
            nThisTDMZZ=nThisTDMZZ+NBASF(iSym1)*NBASF(iSym2)
          End Do
          If (.not.isZero(1)) Then
            call mh5_put_dset_array_real(wfn_sfs_tdm,
     &        TDMIJ,[nThisTDMZZ,1,1], [0,iState-1,jState-1])
          End If
          If (.not.isZero(2)) Then
          call mh5_put_dset_array_real(wfn_sfs_tsdm,
     &      TSDMIJ,[nThisTDMZZ,1,1], [0,iState-1,jState-1])
          End If
          If (IFSO.and.(.not.isZero(3))) Then
            call mh5_put_dset_array_real(wfn_sfs_wetdm,
     &        WDMIJ,[nThisTDMZZ,1,1], [0,iState-1,jState-1])
          End If
        End Do
      End Do
      Call mma_deAllocate(IdxState)
      Call mma_deAllocate(TDMZZ)
      Call mma_deAllocate(TSDMZZ)
      Call mma_deAllocate(WDMZZ)
      Call mma_deAllocate(TDMIJ)
      Call mma_deAllocate(TSDMIJ)
      Call mma_deAllocate(WDMIJ)
      End Subroutine StoreDens
#endif

************************************************************************
*  Transpose_TDM
*
*> @brief
*>   Transpose a transition density matrix in place.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Transpose a transition density matrix, stored in symmetry blocks,
*> replacing the original matrix. The matrices contain only the symmetry
*> blocks that match the total symmetry of the transition.
*>
*> @param[in,out] TDM       Transition density matrix
*> @param[in]     Symmetry  Symmetry of the transition
************************************************************************
      Subroutine Transpose_TDM(TDM,Symmetry)
      Implicit None
      Real*8, Intent(InOut) :: TDM(*)
      Integer, Intent(In) :: Symmetry
#include "rassi.fh"
#include "symmul.fh"
#include "stdalloc.fh"
      Integer :: iSym1,iSym2,nTot,i,j
      Integer :: iBlock(0:8)
      Real*8, Allocatable :: Tmp(:)
* Compute the location of all the stored symmetry blocks
      nTot=0
      iBlock(0)=0
      Do iSym1=1,nSym
        iSym2=Mul(Symmetry,iSym1)
        nTot=nTot+nBasF(iSym1)*nBasF(iSym2)
        iBlock(iSym1)=nTot
      End Do
* Make a copy so we can transpose in place
      Call mma_Allocate(Tmp,nTot,Label='Tmp')
      Call dCopy_(nTot,TDM,1,Tmp,1)
* Transpose symmetry block (a,b) onto symmetry block (b,a)
      Do iSym1=1,nSym
        iSym2=Mul(Symmetry,iSym1)
        Do i=1,nBasF(iSym2)
          Do j=1,nBasF(iSym1)
            TDM(iBlock(iSym2-1)+(j-1)*nBasF(iSym2)+i) =
     &      Tmp(iBlock(iSym1-1)+(i-1)*nBasF(iSym1)+j)
          End Do
        End Do
      End Do
      Call mma_deAllocate(Tmp)
      End Subroutine Transpose_TDM

      End Module Dens2HDF5
