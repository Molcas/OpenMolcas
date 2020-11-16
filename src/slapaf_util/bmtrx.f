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
     &                 ip_rInt,Lbl,Coor,nDim,dMass,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,ip_drInt,
     &                 ipShift,Gx,mTtAtm,iAnr,iOptH,User_Def,
     &                 nStab,jStab,Curvilinear,Numerical,
     &                 DDV_Schlegel,HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,mxdc,iCoSet,lOld,
     &                 rHidden,nFix,nQQ,iIter,Redundant,nqInt,MaxItr,
     &                 nWndw)
      Use Slapaf_Info, Only: Cx
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
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
      Character*100 Get_SuperName
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
      Call Allocate_Work(ipTR,18*nAtom)
      Call FZero(Work(ipTR),18*nAtom)
*
      Call TRPGen(nDim,nAtom,Cx(1,1,iIter),Degen,Smmtrc,mTR,dMass,
     &            .False.,Work(ipTR))
*
      Call Allocate_Work(ipTRnew,3*nAtom*mTR)
      Call FZero(Work(ipTRnew),3*nAtom*mTR)
      i = 0
      Do ix = 1, 3*nAtom
         If (Smmtrc(ix)) Then
            i = i + 1
            iOff = ipTRnew + ix - 1
            call dcopy_(mTR,Work(ipTR+i-1),-nDim,Work(iOff),3*nAtom)
         End If
      End Do
      Call Put_dArray('TR',Work(ipTRnew),3*nAtom*mTR)
      Call Free_Work(ipTRnew)
*
*     Call RecPrt('Work(ipTR)',' ',Work(ipTR),nDim,mTR)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('TabAI','Allo','Inte',ip_TabAI,2*mTtAtm)
      Call GetMem('Vect','Allo','Real',ipVec,3*mTtAtm*nDim)
      Call GetMem('AN','Allo','Inte',ipAN,mTtAtm)
      Call GetMem('Coor','Allo','Real',ipCoor,3*mTtAtm)
*
*-----Generate Grand atoms list
*
      Call GenCoo(Cx(1,1,iIter),nAtom,Work(ipCoor),mTtAtm,Work(ipVec),
     &            Smmtrc,nDim,iAnr,iWork(ipAN),iWork(ip_TabAI),Degen)
*
*---- Are there some hidden frozen atoms ?
*
      nHidden = 0
      ipMMKept = 0
      nMDstep = 0
      If (rHidden.ge.Two) Call Hidden(mTtAtm,ipCoor,ipAN,nHidden,
     &                                rHidden,nMDstep)
*
*-----Generate bond list
*
      ThrB=0.0D0  ! dummy
      mTtAtm = mTtAtm+nHidden
      Call Box(Work(ipCoor),mTtAtm,iWork(ipAN),iOptC,
     &         ddV_Schlegel,ip_TabB,ip_TabA,nBonds,nMax,ThrB)
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
      Call GetMem('EVal','Allo','Real',ipEVal,
     &            (3*mTtAtm)*(3*mTtAtm+1)/2)
      Call GetMem('scr1','Allo','Real',ip_Hss_X,(3*mTtAtm)**2)
      Call GetMem('scr2','Allo','Real',ipScr2,(3*mTtAtm)**2)
*
      If (HSet.or..Not.(Curvilinear.or.User_Def))
     &   Call LNM(Work(ipCoor),mTtAtm,Work(ipEVal),Work(ip_Hss_X),
     &            Work(ipScr2),Work(ipVec),nAtom,nDim,iWork(ipAN),
     &            Smmtrc,Gx,nIter,iOptH,Degen, DDV_Schlegel,
     &            Analytic_Hessian,iOptC,iWork(ip_TabB),iWork(ip_TabA),
     &            nBonds,nMax,nHidden,nMDstep,ipMMKept)
*
      Call GetMem('scr2','Free','Real',ipScr2,(3*mTtAtm)**2)
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
     &                 nLines,nBVec,ipBMx,nAtom,nInter,
     &                 ip_rInt,Lbl,Coor,nDim,dMass,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,ip_drInt,
     &                 Gx,mTtAtm,iAnr,
     &                 nStab,jStab,Numerical,
     &                 HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,mxdc,iCoSet,lOld,
     &                 nFix,mTR,ip_KtB,nQQ,Redundant,nqInt,MaxItr)
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
            Call Box(Work(ipCoor),mTtAtm,iWork(ipAN),iOptC,
     &               ddV_Schlegel,ip_TabB,ip_TabA,nBonds,nMax,ThrB)
         End If
         Call BMtrx_Internal(
     &                 nLines,ipBMx,nAtom,nInter,
     &                 ip_rInt,Coor,nDim,dMass,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,ip_drInt,
     &                 Gx,mTtAtm,iAnr,
     &                 nStab,jStab,Numerical,
     &                 HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,mxdc,iCoSet,lOld,
     &                 nFix,iIter,mTR,Work(ipTR),ip_TabAI,
     &                 ip_TabA,ip_TabB,nBonds,nMax,
     &                 iIter,ip_KtB,nQQ,Redundant,nqInt,MaxItr,nWndw)
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
         Call BMtrx_Cartesian(
     &                 nLines,ipBMx,nAtom,nInter,
     &                 ip_rInt,Coor,nDim,dMass,
     &                 Name,Smmtrc,
     &                 Degen,BSet,HSet,nIter,ip_drInt,
     &                 Gx,mTtAtm,iAnr,
     &                 nStab,jStab,Numerical,
     &                 HWRS,Analytic_Hessian,
     &                 iOptC,PrQ,mxdc,iCoSet,lOld,
     &                 nFix,mTR,Work(ipTR),ipEVal,ip_Hss_X,
     &                 ip_KtB,nQQ,Redundant,nqInt,MaxItr,nWndw)
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
         Call Put_dArray('Hss_X',Work(ip_Hss_X),nDim**2)
         Call Put_dArray('KtB',Work(ip_KtB),nDim*nQQ)
         Call Free_Work(ip_KtB)
      End If
      Call Free_Work(ip_Hss_X)
      Call GetMem('EVal','Free','Real',ipEVal,(3*mTtAtm)*(3*mTtAtm+1)/2)
      Call Free_iWork(ip_TabA)
      Call Free_iWork(ip_TabB)
      Call GetMem('Coor','Free','Real',ipCoor,3*mTtAtm)
      Call GetMem('AN','Free','Inte',ipAN,mTtAtm)
      Call GetMem('Vect','Free','Real',ipVec,3*mTtAtm*nDim)
      Call GetMem('TabAI','Free','Inte',ip_TabAI,2*mTtAtm)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the shift vector in the basis.
*
      If (Bset) Then
         Call GetMem('Shift','Allo','Real',ipShift,nQQ*MaxItr)
         Call FZero(Work(ipShift),nQQ*MaxItr)
         Call ShfANM(nQQ,nIter,Work(ip_rInt),Work(ipShift),iPrint)
      Else
         ipShift=ip_Dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Store B matrices to be used to transform the numerical
*     Hessian to the new basis as the optimization proceeds.
*
      If ((nIter.eq.1.and.BSet).and.
     &    (Get_SuperName().ne.'numerical_gradient')) Then
         Call Allocate_Work(ipBMxOld,3*nAtom*nQQ)
         Call FZero(Work(ipBmxOld),3*nAtom*nQQ)
         call dcopy_(3*nAtom*nQQ,Work(ipBMx),1,Work(ipBMxOld),1)
         Call Put_dArray('BMxOld',Work(ipBMxOld),3*nAtom*nQQ)
         Call Free_Work(ipBMxOld)
         If (mTR.ne.0) Then
            Call Allocate_Work(ipTROld,3*nAtom*mTR)
            Call FZero(Work(ipTROld),3*nAtom*mTR)
#ifdef _DEBUGPRINT_
            Call RecPrt('TRVec',' ',Work(ipTR),3*nAtom,mTR)
#endif
            i = 0
            Do ix = 1, 3*nAtom
               If (Smmtrc(ix)) Then
                  i = i + 1
                  iOff = ipTROld + ix - 1
                  call dcopy_(mTR,Work(ipTR+i-1),-nDim,
     &                           Work(iOff),3*nAtom)
               End If
            End Do
            Call Put_dArray('TROld',Work(ipTROld),3*nAtom*mTR)
            Call Free_Work(ipTROld)
         End If
      End IF
*
*---- Print the B-matrix
*
#ifdef _DEBUGPRINT_
      Call RecPrt(' The BMtrx',' ',Work(ipBMx),3*nAtom,nQQ)
#endif
      Call Free_Work(ipTR)
*                                                                      *
************************************************************************
*                                                                      *
*.... Print out the values of the internal coordinates
*
      If (iPrint.ge.99) Then
         ip = ip_rint + (nIter-1)*nInter -1
         Write (6,*)
         Write (6,*) ' Internal coordinates'
         Write (6,*)
         Write (6,'(1X,A,2X,F10.4)')
     &         (Lbl(iInter),Work(ip+iInter),iInter=1,nQQ)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
