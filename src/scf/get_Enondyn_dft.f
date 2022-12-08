************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Get_Enondyn_dft(nh1,Grad,nGrad,DFTFOCK)
      use SCF_Arrays, only: CMO
      use InfSCF
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8  Grad(nGrad)
      Character*4 DFTFOCK
#include "dcscf.fh"
      Real*8, Allocatable :: F_DFT(:,:), D_DS(:,:)
*
      Erest_xc=0.0d0
      Call mma_allocate(D_DS ,nBT,2,Label='D_DS ')
      D_DS(:,:)=Zero
*
      iOff=1
      jOff=1
      Do iSym=1,nSym
         Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nOcc(iSym,1),
     &                    1.0d0,CMO(iOff,1),nBas(iSym),
     &                          CMO(iOff,1),nBas(iSym),
     &                    0.0d0,D_DS(jOff:,1),nBas(iSym))
         Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nOcc(iSym,2),
     &                    1.0d0,CMO(iOff,2),nBas(iSym),
     &                          CMO(iOff,2),nBas(iSym),
     &                    0.0d0,D_DS(jOff:,2),nBas(iSym))
         Do j=1,nBas(iSym)
            Do i=1,j-1
               ji=j*(j-1)/2+i
               iDji=iOff-1+ji
               D_DS(iDji,1)=2.0d0*D_DS(iDji,1)
               D_DS(iDji,2)=2.0d0*D_DS(iDji,2)
            End Do
         End Do
         iOff=iOff+nBas(iSym)*nOrb(iSym)
         jOff=jOff+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
*
*----------------------------------------------------------------------*
      Call Get_Fmat_nondyn(D_DS(:,1),D_DS(:,2),nBT,.true.)
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
      Call mma_allocate(F_DFT,nBT,2,Label='F_DFT')
      Call Get_Exc_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,nBT,2,
     &                     KSDFT)
*----------------------------------------------------------------------*
*
      Call mma_deallocate(D_DS)
      Call mma_deallocate(F_DFT)
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Get_Exc_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,nBT,nD,
     &                       KSDFT)
      use nq_Info, only: Dens_I, Grad_I, Tau_I
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "debug.fh"
      Real*8  Grad(nGrad)
      Logical Do_MO,Do_TwoEl,Do_Grad
      Character*4 DFTFOCK
      Character*80  KSDFT
#include "dcscf.fh"
      Real*8 :: D_DS(nBT,nD), F_DFT(nBT,nD)
*
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
*     DFT functionals, compute integrals over the potential
*
      Func            =Zero
      Dens_I          =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
      Do_Grad=.false.
*
      nFckDim=2
*                                                                      *
************************************************************************
*                                                                      *
      Call Driver(KSDFT,Do_Grad,Func,Grad,nGrad,
     &            Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nFckDim,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Erest_xc=Erest_xc-Func
*
#ifdef _DEBUGPRINT_
      write(6,*) ' XC-part of energy-restoring term : ',-Func
      write(6,*)
      write(6,*) ' XC-potentials: (itri,F_alpha,F_beta)'
      write(6,*)
      Do i=1,nh1
        Write(6,'(i4,3f22.16)') i,F_DFT(i,1),F_DFT(i,2)
      End Do
#endif
*
      Return
      End
