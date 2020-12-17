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
      Subroutine BMtrx(nsAtom,nInter,Coor,nIter,mTtAtm,nQQ,iIter,nWndw)
      Use Slapaf_Info, Only: Cx, ANr, Shift, qInt, KtB, BMx, Smmtrc,
     &                       Lbl
      Use Slapaf_Parameters, only: Curvilinear, Redundant, nDimBC,
     &                             User_Def, MaxItr, BSet, HSet,
     &                             rHidden, lOld, Numerical
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 Coor(3,nsAtom)
      External Get_SuperName
      Character(LEN=100) Get_SuperName
      Integer, Allocatable:: TabB(:,:), TabA(:,:,:), TabAI(:,:), AN(:)
      Real*8, Allocatable:: TR(:), TRNew(:), TROld(:), Scr2(:),
     &                      Vec(:,:), Coor2(:,:), EVal(:), Hss_X(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
        Subroutine Box(Coor,nsAtom,iANr,TabB,TabA,nBonds,
     &                nMax)
        Integer nsAtom
        Real*8 Coor(3,nsAtom)
        Integer iANr(nsAtom)
        Integer, Allocatable:: TabB(:,:), TabA(:,:,:)
        Integer nBonds, nMax
        End Subroutine Box
        Subroutine Hidden(mTtAtm,Coor,AN,nHidden)
        Integer mTtAtm
        Real*8, Allocatable:: Coor(:,:)
        Integer, Allocatable:: AN(:)
        Integer nHidden
        End Subroutine Hidden
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      iRout=133
      iPrint=nPrint(iRout)
*
      Lu=6
*                                                                      *
************************************************************************
*                                                                      *
*     iIter: point at the reference geometry, used for non-redundant
*            internal coordinate to define the K-matrix and to generate
*            the raw model Hessian and TR vectors.
*
      If (Numerical) Then
         iIter = 1               ! Numerical Hessian Computation
      Else If (iIter.eq.0) Then
         If (.Not.BSet) Then
            iIter = nIter-1      ! Compute cartesian Structure
         Else
            iIter = nIter        ! Normal Computation
         End If
      End If
*
#ifdef _DEBUGPRINT_
      Write (6,*) ' Actual structure from iteration',iIter
      Write (6,*) ' Last structure from iteration',nIter
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Find the translational and rotational eigenvectors for the
*     current structure.
*
      Call mma_allocate(TR,18*nsAtom,Label='TR')
      TR(:)=Zero
*
      Call TRPGen(nDimBC,nsAtom,Cx(1,1,iIter),mTR,.False.,TR)
*
      Call mma_allocate(TRnew,3*nsAtom*mTR,Label='TRNew')
      TRNew(:)=Zero
      i = 0
      Do ix = 1, 3*nsAtom
         iAtom = (ix+2)/3
         ixyz = ix - (iAtom-1)*3
         If (Smmtrc(ixyz,iAtom)) Then
            i = i + 1
            call dcopy_(mTR,TR(i),-nDimBC,TRNew(ix),3*nsAtom)
         End If
      End Do
      Call Put_dArray('TR',TRnew,3*nsAtom*mTR)
      Call mma_deallocate(TRnew)
*
*     Call RecPrt('TR',' ',TR,nDimBC,mTR)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(TabAI,2,mTtAtm,Label='TabAI')
      Call mma_allocate(Vec,3*mTtAtm,nDimBC,Label='Vec')
      Call mma_allocate(AN,mTtAtm,Label='AN')
      Call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')
*
*-----Generate Grand atoms list
*
      Call GenCoo(Cx(1,1,iIter),nsAtom,Coor2,mTtAtm,Vec,nDimBC,ANr,
     &            AN,TabAI)
*
*---- Are there some hidden frozen atoms ?
*
      nHidden = 0
      If (rHidden.ge.Two) Call Hidden(mTtAtm,Coor2,AN,nHidden)
*
*-----Generate bond list
*
      mTtAtm = mTtAtm+nHidden
      Call Box(Coor2,mTtAtm,AN,TabB,TabA,nBonds,nMax)
      mTtAtm = mTtAtm-nHidden
*                                                                      *
************************************************************************
*                                                                      *
*---- First compute the approximate Hessian in cartesians
*
*     OBSERVE if the analytic Hessian is available it will be used
*     rather than the model Hessian.
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the raw Cartesian Hessian
*
      Call mma_allocate(EVal,(3*mTtAtm)*(3*mTtAtm+1)/2,Label='EVal')
      Call mma_Allocate(Hss_X,(3*mTtAtm)**2,Label='Hss_X')
      Call mma_allocate(Scr2,(3*mTtAtm)**2,Label='Scr2')
*
      If (HSet.or..Not.(Curvilinear.or.User_Def))
     &   Call LNM(Coor2,mTtAtm,EVal,Hss_X,Scr2,Vec,nsAtom,nDimBC,AN,
     &            nIter,
     &            TabB,TabA,nBonds,nMax,nHidden)
*
      Call mma_deallocate(Scr2)
      Call mma_deallocate(Coor2)
*                                                                      *
************************************************************************
*                                                                      *
*     The internal coordinates are of either two sets.
*
*              i)  user supplied internal coordinates
*
*             or
*
*             ii)  automatic internal coordinates
*                    a) cartesian coordinates (lnm)
*                    b) non-redundant internal coordinates (nrc)
*                       1) Conventional
*                       2) Curvature weighted   (HWRS)
*
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (User_Def) Then
*                                                                      *
************************************************************************
*                                                                      *
         Call BMtrx_User_Defined(
     &                 nsAtom,nInter,Lbl,Coor,nDimBC,
     &                 nIter,
     &                 mTR,nQQ)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Curvilinear) Then
*                                                                      *
************************************************************************
*                                                                      *
         If (Redundant) Then
            Call WarningMessage(2,
     &           ' Bmtrx: Redundant option not implemented yet.')
            Call Abend()
         End If
*
*------- Re-generate the bonds if there were hidden atoms
*
         If (nHidden.ne.0) Then
            Call Box(Coor2,mTtAtm,AN,TabB,TabA,
     &               nBonds,nMax)
         End If
         Call BMtrx_Internal(
     &                 nsAtom,nDimBC,
     &                 nIter,
     &                 mTtAtm,
     &                 iIter,mTR,TR,TabAI,
     &                 TabA,TabB,nBonds,nMax,
     &                 iIter,nQQ,nWndw)
*
*------- Set the Labels for internal coordinates.
*
         Do i = 1, nInter
            Write (Lbl(i),'(A,I3.3,A)') 'nrc',i,'  '
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      Else   ! Cartesian coordinates
*                                                                      *
************************************************************************
*                                                                      *
         Call BMtrx_Cartesian(nsAtom,nInter,nDimBC,
     &                        nIter,
     &                        mTtAtm,mTR,TR,EVal,Hss_X,
     &                        nQQ,nWndw)
*
*------- Set the Labels for cartesian normal modes.
*
         Do i = 1, nQQ
            Write (Lbl(i),'(A,I3.3,A)') 'lnm',i,'  '
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If ((BSet.and.HSet.and..NOT.lOld)) Then
         Call Put_dArray('Hss_X',Hss_X,nDimBC**2)
         Call Put_dArray('KtB',KtB,nDimBC*nQQ)
         Call mma_deallocate(KtB)
      End If
      Call mma_deallocate(Hss_X)
      Call mma_deallocate(EVal)
      Call mma_deallocate(TabA)
      Call mma_deallocate(TabB)
      Call mma_deallocate(AN)
      Call mma_deallocate(Vec)
      Call mma_deallocate(TabAI)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the shift vector in the basis.
*
      If (Bset) Then
         Call mma_allocate(Shift,nQQ,MaxItr,Label='Shift')
         Shift(:,:)=Zero
         Call ShfANM(nQQ,nIter,qInt,Shift,iPrint)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Store B matrices to be used to transform the numerical
*     Hessian to the new basis as the optimization proceeds.
*
      If ((nIter.eq.1.and.BSet).and.
     &    (Get_SuperName().ne.'numerical_gradient')) Then

         Call Put_dArray('BMxOld',BMx,3*nsAtom*nQQ)

         If (mTR.ne.0) Then
            Call mma_allocate(TROld,3*nsAtom*mTR,Label='TROld')
            TROld(:)=Zero
#ifdef _DEBUGPRINT_
            Call RecPrt('TRVec',' ',TR,3*nsAtom,mTR)
#endif
            i = 0
            Do ix = 1, 3*nsAtom
               iAtom = (ix+2)/3
               ixyz = ix - (iAtom-1)*3
               If (Smmtrc(ixyz,iAtom)) Then
                  i = i + 1
                  call dcopy_(mTR,TR(i),-nDimBC,TROld(ix),3*nsAtom)
               End If
            End Do
            Call Put_dArray('TROld',TROld,3*nsAtom*mTR)
            Call mma_deallocate(TROld)
         End If
      End IF
*
*---- Print the B-matrix
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' The BMtrx',' ',BMx,3*nsAtonsAtom,nQQ)
#endif
      Call mma_deallocate(TR)
*                                                                      *
************************************************************************
*                                                                      *
*.... Print out the values of the internal coordinates
*
      If (iPrint.ge.99) Then
         Write (6,*)
         Write (6,*) ' Internal coordinates'
         Write (6,*)
         Write (6,'(1X,A,2X,F10.4)')
     &         (Lbl(iInter),qInt(iInter,nIter),iInter=1,nQQ)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
