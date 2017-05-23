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
* Copyright (C) 1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine LinSer()
************************************************************************
*                                                                      *
*     purpose: Line Search for QNR steps                               *
*              between dE and dE(1st order)                            *
*                                                                      *
*     input:                                                           *
*       Ovrlp   : overlap integrals                                    *
*       mBT     : size of overlap integrals, triangular storage        *
*       CMO     : Molecular Orbital coefficients                       *
*       mBB     : size of the CMO array                                *
*                                                                      *
*     output:                                                          *
*       Enew    : new actual SCF energy after Line Search              *
*       En1V    : variational 1el energy after Line Search             *
*       En2V    : variational 2el energy after Line Search             *
*       FstItr  : used for semidirect Fock matrix construction in      *
*                 Drv2El_dscf ...                                      *
*     In this Subroutine, the CMOs, Densities and TwoHam Matrices are  *
*     altered in the Micro Iteractions, hopefully towards lower SCF    *
*     energy values...                                                 *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*     calls to: RotMOs,DMat,PMat,GrdClc                                *
*               uses SubRoutines and Functions from Module lnklst.f    *
*               -linked list implementation to store series of vectors *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: corresponds to large extends to Jeppe s LinSe2          *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "file.fh"
#include "llists.fh"
#include "print.fh"
*
*---- declaration of local variables
*
      Integer jpgrd
      Real*8, Dimension(:,:), Allocatable:: Scr
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     declarations of functions
*
      Integer LstPtr
*                                                                      *
************************************************************************
*                                                                      *
*
*                                                                      *
************************************************************************
*                                                                      *
*     initialization stuff
*
      nD = iUHF + 1
      Call mma_allocate(Scr,nOV,nD,Label='Scr')
*                                                                      *
*     compute gradient of actual point g(x(n))
*
*     Call GrdClc('Lst',.true.)
*                                                                      *
************************************************************************
*                                                                      *
*     get Pointer to the gradient back from LList (only address)
      jpgrd=LstPtr(LuGrd,iter,LLGrad)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute
*
*             dg(n-1)= g(n) - g(n-1) : g(n)=Work(jpgrd)
*
*     and update the data on files and in memory.
*
*     (5) get g(n-1)
*
      Call GetNod(iter-1,LLGrad,inode)
      If (inode.eq.0) Then
         Write (6,*) 'inode.eq.0'
         Call Abend()
      End If
      Call iVPtr(LuGrd,Scr,nOV*nD,inode)
*
*     (6) compute dg(n-1)=g(n)-g(n-1)
*
      Call DaXpY_(nOV*nD,-One,Work(jpgrd),1,Scr,1)
      Call DScal_(nOV*nD,-One,Scr,1)
*
*     (7) put dg(n-1) on its LList. Will be used later for DIIS or RS-RFO
*
      Call PutVec(Scr,nOV*nD,LudGd,iter-1,MemRsv,'NOOP',LLdGrd)
*
*     and free temp. allocated memory
*
      Call mma_deallocate(Scr)
*
      Return
      End
*define _DEBUG_X_
#ifdef _DEBUG_X_
      SubRoutine LinSer(FstItr,nD,Ovrlp,mBT,CMO,mBB)
************************************************************************
*                                                                      *
*     purpose: Line Search for QNR steps                               *
*              between dE and dE(1st order)                            *
*                                                                      *
*     input:                                                           *
*       Ovrlp   : overlap integrals                                    *
*       mBT     : size of overlap integrals, triangular storage        *
*       CMO     : Molecular Orbital coefficients                       *
*       mBB     : size of the CMO array                                *
*                                                                      *
*     output:                                                          *
*       Enew    : new actual SCF energy after Line Search              *
*       En1V    : variational 1el energy after Line Search             *
*       En2V    : variational 2el energy after Line Search             *
*       FstItr  : used for semidirect Fock matrix construction in      *
*                 Drv2El_dscf ...                                      *
*     In this Subroutine, the CMOs, Densities and TwoHam Matrices are  *
*     altered in the Micro Iteractions, hopefully towards lower SCF    *
*     energy values...                                                 *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*     calls to: RotMOs,DMat,PMat,GrdClc                                *
*               uses SubRoutines and Functions from Module lnklst.f    *
*               -linked list implementation to store series of vectors *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: corresponds to large extends to Jeppe s LinSe2          *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
#include "file.fh"
#include "llists.fh"
#include "print.fh"
*
*     declaration procedure parameters
      Real*8 Enew,En1V,En2V
      Real*8 Ovrlp(mBT), CMO(mBB,nD)
      Logical FstItr
*
*---- declaration of local variables
*
      Integer jpgrd
      Real*8, Dimension(:,:), Allocatable:: Scr
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     If _DEBUG_X_ is defined this will activate the computation of
*     the gradient of the energy with respect to the matrix elements
*     if the anti-symmetric X matrix. The code use a four point
*     symmetric formula with an adjusted step depending on the analytic
*     value of the gradient. Accuracy of the computed gradients has been
*     observed to devite for values below 0.1D-6. This is a bit of a
*     worry!
*
*define _DEBUG_X_
#ifdef _DEBUG_X_
      Real*8, Dimension(:,:), Allocatable:: Grdnm, CMO_ref
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     declarations of functions
*
      Integer LstPtr
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUG_X_
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     This code is here to debug the option which use the
*     anti-symmetric X matrix to represent the orbital
*     rotations. Here we will recompute the analytically computed
*     gradient numerically.
*
      Call mma_allocate(Grdnm,nOV,nD,Label='Grdnm')
      Call mma_allocate(CMO_ref,nBO,nD,Label='CMO_ref')
*
*     save current CMOs...
*
      call dcopy_(nBO*nD,CMO,1,CMO_ref,1)
*
*     get g(n-1)
*
      Call GetNod(iter-1,LLGrad,inode)
      If (inode.eq.0) Then
         Write (6,*) 'inode.eq.0'
         Call Abend()
      End If
      Call iVPtr(LuGrd,Grdnm,nOV*nD,inode)
*
      Call FZero(Scr,nOV*nD)
      Write (6,'(A)') '  iOV   iD   Grdnm     numerical'
      Do iD = 1, nD
         Do iOV = 1, nOV
*
*           use symmetric two-point formula
*
            EDiff = 0.0D0
*define _TWO_POINT_
#ifdef _TWO_POINT_
            Delta =1.0D-4
            Do i = -1, 1, 2
               Scr(iOV,iD) =  DBLE(i)*Delta
               Call DCopy_(nBO*nD,CMO_ref,1,CMO,1)
               Call RotMOs(Scr,nOV,CMO,nBO,nD,Ovrlp,mBT)
               Call SCF_Energy(FstItr,En1V,En2V,Enew)
*
               EDiff = EDiff + DBLE(i)*ENew
               Scr(iOV,iD)=0.0D0
            End Do
*
            Write (6,'(2I5,2E10.3)')
     &            iOV, iD, Grdnm(iOV,iD), EDiff/(2.0D0*Delta)
#else
            Delta =5.0D-4
            E1=0.0D0
            E2=0.0D0
            E3=0.0D0
            E4=0.0D0
            If (Abs(Grdnm(iOV,iD)).lt.1.0D-09) Delta=5.0D-2
            If (Abs(Grdnm(iOV,iD)).lt.1.0D-12) Delta=5.0D-1
*           Write (*,*) 'Delta=',Delta
            Do i = -2, 2
               If (i.eq.0) Cycle
               Scr(iOV,iD) =  DBLE(i)*Delta
               Call DCopy_(nBO*nD,CMO_ref,1,CMO,1)
               Call RotMOs(Scr,nOV,CMO,nBO,nD,Ovrlp,mBT)
               Call SCF_Energy(FstItr,En1V,En2V,Enew)
*
               If (i.eq.-2) Then
                  E1= 1.0D0*Enew
               Else If(i.eq.-1) Then
                  E2=-1.0D0*Enew
               Else If(i.eq. 1) Then
                  E3= 1.0D0*Enew
               Else
                  E4=-1.0D0*Enew
               End If
               Scr(iOV,iD)=0.0D0
            End Do
*           Write (6,*) E1
*           Write (6,*) E2
*           Write (6,*) E3
*           Write (6,*) E4
            E1 = E1 + E4
            E2 = 8.0D0*(E2 + E3)
            EDIFF = E1 + E2
*
            Write (6,'(2I5,2E10.3)')
     &            iOV, iD, Grdnm(iOV,iD), EDiff/(12.0D0*Delta)
#endif
*
            Write (6,*)
            Delta=1.0D-3
            EDiff=0.0D0
            Do i = -1, 1
               Scr(iOV,iD) =  DBLE(i)*Delta
*
               Call DCopy_(nBO*nD,CMO_ref,1,CMO,1)
               Call RotMOs(Scr,nOV,CMO,nBO,nD,Ovrlp,mBT)
               Call SCF_Energy(FstItr,En1V,En2V,Enew)
*
               If (i.eq.0) Then
                  EDiff = EDiff - 2.0D0*Enew
               Else
                  EDiff = EDiff + Enew
               End If
*
               Scr(iOV,iD)=0.0D0
            End Do
            Write (6,'(2I5,1E10.3)') iOV, iD, EDiff/(Delta**2)
         End Do
      End Do
*
*     After completed work we abend the calculations.
*
      Call mma_deallocate(CMO_ref)
      Call mma_deallocate(Grdnm)
      Call mma_deallocate(Scr)
      Call Abend()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#else
*                                                                      *
************************************************************************
*                                                                      *
*     initialization stuff
*
      Call mma_allocate(Scr,nOV,nD,Label='Scr')
*                                                                      *
*     compute gradient of actual point g(x(n))
*
*     Call GrdClc('Lst',.true.)
*                                                                      *
************************************************************************
*                                                                      *
*     get Pointer to the gradient back from LList (only address)
      jpgrd=LstPtr(LuGrd,iter,LLGrad)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute
*
*             dg(n-1)= g(n) - g(n-1) : g(n)=Work(jpgrd)
*
*     and update the data on files and in memory.
*
*     (5) get g(n-1)
*
      Call GetNod(iter-1,LLGrad,inode)
      If (inode.eq.0) Then
         Write (6,*) 'inode.eq.0'
         Call Abend()
      End If
      Call iVPtr(LuGrd,Scr,nOV*nD,inode)
*
*     (6) compute dg(n-1)=g(n)-g(n-1)
*
      Call DaXpY_(nOV*nD,-One,Work(jpgrd),1,Scr,1)
      Call DScal_(nOV*nD,-One,Scr,1)
*
*     (7) put dg(n-1) on its LList. Will be used later for DIIS or RS-RFO
*
      Call PutVec(Scr,nOV*nD,LudGd,iter-1,MemRsv,'NOOP',LLdGrd)
*
*     and free temp. allocated memory
*
      Call mma_deallocate(Scr)
#endif
*
      Return
      End
#endif
