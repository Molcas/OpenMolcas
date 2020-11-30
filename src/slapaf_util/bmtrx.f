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
      Subroutine BMtrx(nLines,nBVec,ipBMx,nAtom,nInter,
     &                 Lbl,Coor,nDim,dMass,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,
     &                 Gx,mTtAtm,iAnr,iOptH,User_Def,
     &                 nStab,jStab,Curvilinear,Numerical,
     &                 DDV_Schlegel,HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,mxdc,iCoSet,lOld,
     &                 rHidden,nFix,nQQ,iIter,Redundant,MaxItr,
     &                 nWndw)
      Use Slapaf_Info, Only: Cx, Shift, qInt
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom),  dMass(nAtom), Degen(3*nAtom),
     &       Gx(3*nAtom,nIter)
      Character Lbl(nInter)*8,Name(nAtom)*(LENIN)
      Integer   iAnr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Logical Smmtrc(3*nAtom), BSet, HSet, Redundant,
     &        User_Def, Curvilinear, Numerical, DDV_Schlegel,
     &        HWRS, Analytic_Hessian, PrQ, lOld
      External Get_SuperName
      Character(LEN=100) Get_SuperName
      Integer, Allocatable:: TabB(:,:), TabA(:,:,:), TabAI(:,:), AN(:)
      Integer, External:: ip_of_iWork
      Real*8, Allocatable:: TR(:), TRNew(:), TROld(:), Scr2(:),
     &                      Vec(:,:), Coor2(:,:), EVal(:), Hss_X(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
        Subroutine Box(Coor,nAtoms,iANr,iOptC,Schlegel,TabB,TabA,nBonds,
     &                nMax)
        Integer nAtoms
        Real*8 Coor(3,nAtoms)
        Integer iANr(nAtoms)
        Integer iOptC
        Logical Schlegel
        Integer, Allocatable:: TabB(:,:), TabA(:,:,:)
        Integer nBonds, nMax
        End Subroutine Box
        Subroutine Hidden(mTtAtm,Coor,AN,nHidden,rHidden)
        Integer mTtAtm
        Real*8, Allocatable:: Coor(:,:)
        Integer, Allocatable:: AN(:)
        Integer nHidden
        Real*8 rHidden
        End Subroutine Hidden
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
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
      Call mma_allocate(TR,18*nAtom,Label='TR')
      TR(:)=Zero
*
      Call TRPGen(nDim,nAtom,Cx(1,1,iIter),Degen,Smmtrc,mTR,dMass,
     &            .False.,TR)
*
      Call mma_allocate(TRnew,3*nAtom*mTR,Label='TRNew')
      TRNew(:)=Zero
      i = 0
      Do ix = 1, 3*nAtom
         If (Smmtrc(ix)) Then
            i = i + 1
            call dcopy_(mTR,TR(i),-nDim,TRNew(ix),3*nAtom)
         End If
      End Do
      Call Put_dArray('TR',TRnew,3*nAtom*mTR)
      Call mma_deallocate(TRnew)
*
*     Call RecPrt('TR',' ',TR,nDim,mTR)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(TabAI,2,mTtAtm,Label='TabAI')
      Call mma_allocate(Vec,3*mTtAtm,nDim,Label='Vec')
      Call mma_allocate(AN,mTtAtm,Label='AN')
      Call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')
*
*-----Generate Grand atoms list
*
      Call GenCoo(Cx(1,1,iIter),nAtom,Coor2,mTtAtm,Vec,Smmtrc,nDim,iAnr,
     &            AN,TabAI,Degen)
*
*---- Are there some hidden frozen atoms ?
*
      nHidden = 0
      If (rHidden.ge.Two) Call Hidden(mTtAtm,Coor2,AN,nHidden,rHidden)
*
*-----Generate bond list
*
      mTtAtm = mTtAtm+nHidden
      Call Box(Coor2,mTtAtm,AN,iOptC,ddV_Schlegel,TabB,TabA,nBonds,nMax)
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
     &   Call LNM(Coor2,mTtAtm,EVal,Hss_X,Scr2,Vec,nAtom,nDim,AN,Smmtrc,
     &            nIter,iOptH,Degen, DDV_Schlegel,Analytic_Hessian,
     &            iOptC,TabB,TabA,nBonds,nMax,nHidden)
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
     &                 nLines,nBVec,ipBMx,nAtom,nInter,Lbl,Coor,nDim,
     &                 dMass,Name,Smmtrc,Degen,BSet,HSet,nIter,Gx,
     &                 nStab,jStab,Numerical,Analytic_Hessian,
     &                 iOptC,mxdc,lOld,
     &                 nFix,mTR,ip_KtB,nQQ,Redundant,MaxItr)
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
            Call Box(Coor2,mTtAtm,AN,iOptC,ddV_Schlegel,TabB,TabA,
     &               nBonds,nMax)
         End If
         ip_TabA = ip_of_iWork(TabA(1,0,1))
         ip_TabB = ip_of_iWork(TabB(1,1))
         ip_TabAI= ip_of_iWork(TabAI(1,1))
         Call BMtrx_Internal(
     &                 ipBMx,nAtom,
     &                 nDim,dMass,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,
     &                 Gx,mTtAtm,iAnr,
     &                 nStab,jStab,Numerical,
     &                 HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,iCoSet,lOld,
     &                 iIter,mTR,TR,ip_TabAI,
     &                 ip_TabA,ip_TabB,nBonds,nMax,
     &                 iIter,ip_KtB,nQQ,MaxItr,nWndw)
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
         Call BMtrx_Cartesian(ipBMx,nAtom,nInter,nDim,Name,
     &                        Smmtrc,Degen,BSet,HSet,nIter,
     &                        Gx,mTtAtm,PrQ,lOld,mTR,TR,EVal,Hss_X,
     &                        ip_KtB,nQQ,Redundant,MaxItr,nWndw)
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
         Call Put_dArray('Hss_X',Hss_X,nDim**2)
         Call Put_dArray('KtB',Work(ip_KtB),nDim*nQQ)
         Call Free_Work(ip_KtB)
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

         Call Put_dArray('BMxOld',Work(ipBMx),3*nAtom*nQQ)

         If (mTR.ne.0) Then
            Call mma_allocate(TROld,3*nAtom*mTR,Label='TROld')
            TROld(:)=Zero
#ifdef _DEBUGPRINT_
            Call RecPrt('TRVec',' ',TR,3*nAtom,mTR)
#endif
            i = 0
            Do ix = 1, 3*nAtom
               If (Smmtrc(ix)) Then
                  i = i + 1
                  call dcopy_(mTR,TR(i),-nDim,TROld(ix),3*nAtom)
               End If
            End Do
            Call Put_dArray('TROld',TROld,3*nAtom*mTR)
            Call mma_deallocate(TROld)
         End If
      End IF
*
*---- Print the B-matrix
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' The BMtrx',' ',Work(ipBMx),3*nAtom,nQQ)
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
