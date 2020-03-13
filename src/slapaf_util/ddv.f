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
      subroutine ddV(Cart,nAtoms,Hess,iANr,Schlegel,iOptC,
     &               iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "sbs.fh"
      Real*8 Cart(3,nAtoms+nHidden),Hess((3*nAtoms)*(3*nAtoms+1)/2)
      Integer   iANr(nAtoms+nHidden), iTabBonds(3,nBonds),
     &          iTabAtoms(2,0:nMax,nAtoms+nHidden)
      Logical Schlegel
      Call QEnter('ddV')
*
*  Temporary big hessian
*
************************************************************************
*                                                                      *
*#define _DEBUG_
*                                                                      *
************************************************************************
      If (nHidden.gt.0) Then
         nTot = nAtoms+nHidden
         Call Allocate_Work(ipHBig,(3*nTot)*(3*nTot+1)/2)
*
* Temporary turn on the translational/rotational invariance
*
         iSBS = iEOr(iSBS,2**7)
         iSBS = iEOr(iSBS,2**8)
         Call ddV_(Cart,nTot,Work(ipHBig),iANr,Schlegel,iOptC,iTabBonds,
     &             iTabAtoms,nBonds,nMax,nHidden)
         iSBS = iOr(iSBS,2**7)
         iSBS = iOr(iSBS,2**8)
         Call dCopy_((3*nAtoms)*(3*nAtoms+1)/2,Work(ipHBig),1,Hess,1)
#ifdef _DEBUG_
         write(6,*) 'DDV: Improved Hessian'
         Call RecPrt('Coord (with hidden atoms):',' ',Cart,3,nTot)
         Call TriPrt('Hessian (hidden atoms):',' ',Work(ipHBig),3*nTot)
         Call TriPrt('Hessian (normal):',' ',Hess,3*nAtoms)
#endif
         Call Free_Work(ipHBig)
      Else
         Call ddV_(Cart,nAtoms,Hess,iANr,Schlegel,iOptC,iTabBonds,
     &             iTabAtoms,nBonds,nMax,nHidden)
      End If
      Call QExit('ddV')
      End
*
      Subroutine ddV_(Cart,nAtoms,Hess,iANr,Schlegel,iOptC,iTabBonds,
     &               iTabAtoms,nBonds,nMax,nHidden)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "sbs.fh"
#include "WrkSpc.fh"
      Real*8 Cart(3,nAtoms),rij(3),rjk(3),rkl(3),
     &       Hess((3*nAtoms)*(3*nAtoms+1)/2),si(3),sj(3),sk(3),
     &       sl(3),sm(3),x(2),y(2),z(2),
     &       xyz(3,4), C(3,4), Dum(1),
     &       ril(3), rik(3)
      Integer   iANr(nAtoms), iTabBonds(3,nBonds),
     &          iTabAtoms(2,0:nMax,nAtoms), iOper(0:7)
      Logical Schlegel, MinBas, Help, TransVar, RotVar, Torsion_Check,
     &        Invariant(3)
*
      Real*8 Trans(3), RotVec(3), RotMat(3,3)
*
#include "warnings.fh"
#define _FMIN_
#define _VDW_
#include "ddvdt.fh"
#include "ddvdt_RF.fh"
#define _SCHLEGEL_
#include "ddvdt_bond.fh"
#include "bondtypes.fh"
#include "ddvdt_bend.fh"
#include "ddvdt_trsn.fh"
#include "ddvdt_outofp.fh"
#include "constants.fh"
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Write (6,*) 'ddV: nBonds=',nBonds
      nqR=0
      nqB=0
      nqT=0
      nqO=0
      Do iAtom = 1, nAtoms
*
         nNeighbor_i = iTabAtoms(1,0,iAtom)
         Write (6,*) 'iAtom,nNeighbor=',iAtom,nNeighbor_i
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
      iRout=140
      iPrint=nPrint(iRout)
      Call QEnter('ddV_')
*
      f_const_min_=f_const_min * 1.0D-1
      f_const=0.0D0
*
      bohr=CONST_BOHR_RADIUS_IN_SI_ * 1.0D+10
      MinBas=.False.
      If (MinBas) Then
         Fact=1.3d0
      Else
         Fact=One
      End If
      rZero=1.0d-10
      n3=3*nAtoms
*
      call dcopy_((n3*(n3+1)/2),[Zero],0,Hess,1)
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian at start','(12f8.3)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      TransVar=iAnd(iSBS,2**7).eq. 2**7
      RotVar  =iAnd(iSBS,2**8).eq. 2**8
*
      Call Get_iScalar('NSYM',nSym)
      Call Get_iArray('Symmetry operations',iOper,nSym)
      Call Allocate_Work(ip_xMass,nAtoms)
      Call Get_Mass_All(Work(ip_xMass),nAtoms-nHidden)
      Do iAtom=nAtoms-nHidden+1,nAtoms
         Work(ip_xMass+iAtom-1)=rMass(iANr(iAtom))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Hessian for translational coordinates in case of a RF calculation.
*
      If (.Not.TransVar) Go To 777
*
      Do ixyz = 1, 3
         Invariant(ixyz)=.False.
         iTest=2**(ixyz-1)
         Do iSym = 0, nSym-1
            If (iOper(iSym).eq.iTest) Invariant(ixyz)=.True.
         End Do
      End Do
      If (Invariant(1).and.Invariant(2).and.Invariant(3)) Go To 777
*
      Fact=One
      If (.Not.RotVar) Fact=2.0D-2
*
      TMass=Zero
      Do iAtom = 1, nAtoms
         TMass=TMass+Work(ip_xMass+iAtom-1)
      End Do
      Do iAtom = 1, nAtoms
         f1=Work(ip_xMass+iAtom-1)/TMass
         Do jAtom = 1, iAtom-1
            f2=Work(ip_xMass+jAtom-1)/TMass
*
            f_const=Max(Trans_Const,f_const_Min_)
            gmm=Fact*f_const*f1*f2
*
            If (.Not.Invariant(1))
     &      Hess(LHR(1,iAtom,1,jAtom))=Hess(LHR(1,iAtom,1,jAtom))+gmm
            If (.Not.Invariant(2))
     &      Hess(LHR(2,iAtom,2,jAtom))=Hess(LHR(2,iAtom,2,jAtom))+gmm
            If (.Not.Invariant(3))
     &      Hess(LHR(3,iAtom,3,jAtom))=Hess(LHR(3,iAtom,3,jAtom))+gmm
*
         End Do
*
         f_const=Max(Trans_Const,f_const_Min_)
         gmm=Fact*f_const*f1*f1
*
         If (.Not.Invariant(1))
     &   Hess(LHR(1,iAtom,1,iAtom))=Hess(LHR(1,iAtom,1,iAtom))+gmm
         If (.Not.Invariant(1))
     &   Hess(LHR(2,iAtom,2,iAtom))=Hess(LHR(2,iAtom,2,iAtom))+gmm
         If (.Not.Invariant(1))
     &   Hess(LHR(3,iAtom,3,iAtom))=Hess(LHR(3,iAtom,3,iAtom))+gmm
*
      End Do
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian after Translation','(6f12.7)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Hessian for rotational coordinates
*
 777  Continue
      If (.Not.RotVar) Go To 778
*
      Do ixyz = 1, 3
         Invariant(ixyz)=.False.
         If (ixyz.eq.1) Then
            iTest=6
         Else If (ixyz.eq.2) Then
            iTest=5
         Else
            iTest=3
         End If
         Do iSym = 0, nSym-1
            If (iOper(iSym).eq.iTest) Invariant(ixyz)=.True.
         End Do
      End Do
      If (Invariant(1).and.Invariant(2).and.Invariant(3)) Go To 778
*
c     If (nAtoms.le.2) Then
c        Call WarningMessage(2,'Error in ddV')
c        Write (6,*)
c        Write (6,*) ' Warning!'
c        Write (6,*) ' Rotational internal coordinates not implemented'
c    &             //' for fewer than 3 atoms!'
c        Write (6,*) ' Add dummy atoms to your input and try again!'
c        Write (6,*)
c        Call Quit(_RC_GENERAL_ERROR_)
c     End If
      Call Allocate_Work(ip_Grad,3*3*nAtoms)
      Call Allocate_Work(ip_CurrXYZ,3*nAtoms)
      Do iAtom = 1, nAtoms
         If (iANr(iAtom).le.0) Then
            Work(ip_xMass-1+iAtom) = 1.0D-10
         End If
      End Do
      nOrder=1
      Call FZero(Trans,3)
      Call FZero(RotVec,3)
      call dcopy_(3*nAtoms,Cart,1,Work(ip_CurrXYZ),1)
      Call RotDer(nAtoms,Work(ip_xMass),Work(ip_CurrXYZ),Cart,
     &            Trans,RotAng,
     &            RotVec,RotMat,nOrder,Work(ip_Grad),dum,dum,dum)
      Call Free_Work(ip_CurrXYZ)
*
*
      Do iAtom = 1, nAtoms
         ix = (iAtom-1)*3 + 1
         iy = (iAtom-1)*3 + 2
         iz = (iAtom-1)*3 + 3
         dO1_dx1 =  Work( (ix-1)*3 +     ip_Grad)
         dO2_dx1 =  Work( (ix-1)*3 + 1 + ip_Grad)
         dO3_dx1 =  Work( (ix-1)*3 + 2 + ip_Grad)
         dO1_dy1 =  Work( (iy-1)*3 +     ip_Grad)
         dO2_dy1 =  Work( (iy-1)*3 + 1 + ip_Grad)
         dO3_dy1 =  Work( (iy-1)*3 + 2 + ip_Grad)
         dO1_dz1 =  Work( (iz-1)*3 +     ip_Grad)
         dO2_dz1 =  Work( (iz-1)*3 + 1 + ip_Grad)
         dO3_dz1 =  Work( (iz-1)*3 + 2 + ip_Grad)
         If (Invariant(1)) Then
            dO1_dx1 = Zero
            dO1_dy1 = Zero
            dO1_dz1 = Zero
         End If
         If (Invariant(2)) Then
            dO2_dx1 = Zero
            dO2_dy1 = Zero
            dO2_dz1 = Zero
         End If
         If (Invariant(3)) Then
            dO3_dx1 = Zero
            dO3_dy1 = Zero
            dO3_dz1 = Zero
         End If
         Do jAtom = 1, iAtom-1
            jx = (jAtom-1)*3 + 1
            jy = (jAtom-1)*3 + 2
            jz = (jAtom-1)*3 + 3
            dO1_dx2 =  Work( (jx-1)*3 +     ip_Grad)
            dO2_dx2 =  Work( (jx-1)*3 + 1 + ip_Grad)
            dO3_dx2 =  Work( (jx-1)*3 + 2 + ip_Grad)
            dO1_dy2 =  Work( (jy-1)*3 +     ip_Grad)
            dO2_dy2 =  Work( (jy-1)*3 + 1 + ip_Grad)
            dO3_dy2 =  Work( (jy-1)*3 + 2 + ip_Grad)
            dO1_dz2 =  Work( (jz-1)*3 +     ip_Grad)
            dO2_dz2 =  Work( (jz-1)*3 + 1 + ip_Grad)
            dO3_dz2 =  Work( (jz-1)*3 + 2 + ip_Grad)
            If (Invariant(1)) Then
               dO1_dx2 = Zero
               dO1_dy2 = Zero
               dO1_dz2 = Zero
            End If
            If (Invariant(2)) Then
               dO2_dx2 = Zero
               dO2_dy2 = Zero
               dO2_dz2 = Zero
            End If
            If (Invariant(3)) Then
               dO3_dx2 = Zero
               dO3_dy2 = Zero
               dO3_dz2 = Zero
            End If
*
            f_const=Max(Rot_Const,f_const_Min_)
            Hess(LHR(1,iAtom,1,jAtom))=Hess(LHR(1,iAtom,1,jAtom))
     &          + f_const  *(dO1_dx1*dO1_dx2
     &                      +dO2_dx1*dO2_dx2
     &                      +dO3_dx1*dO3_dx2)
            Hess(LHR(1,iAtom,2,jAtom))=Hess(LHR(1,iAtom,2,jAtom))
     &          + f_const  *(dO1_dx1*dO1_dy2
     &                      +dO2_dx1*dO2_dy2
     &                      +dO3_dx1*dO3_dy2)
            Hess(LHR(1,iAtom,3,jAtom))=Hess(LHR(1,iAtom,3,jAtom))
     &          + f_const  *(dO1_dx1*dO1_dz2
     &                      +dO2_dx1*dO2_dz2
     &                      +dO3_dx1*dO3_dz2)
            Hess(LHR(2,iAtom,1,jAtom))=Hess(LHR(2,iAtom,1,jAtom))
     &          + f_const  *(dO1_dy1*dO1_dx2
     &                      +dO2_dy1*dO2_dx2
     &                      +dO3_dy1*dO3_dx2)
            Hess(LHR(2,iAtom,2,jAtom))=Hess(LHR(2,iAtom,2,jAtom))
     &          + f_const  *(dO1_dy1*dO1_dy2
     &                      +dO2_dy1*dO2_dy2
     &                      +dO3_dy1*dO3_dy2)
            Hess(LHR(2,iAtom,3,jAtom))=Hess(LHR(2,iAtom,3,jAtom))
     &          + f_const  *(dO1_dy1*dO1_dz2
     &                      +dO2_dy1*dO2_dz2
     &                      +dO3_dy1*dO3_dz2)
            Hess(LHR(3,iAtom,1,jAtom))=Hess(LHR(3,iAtom,1,jAtom))
     &          + f_const  *(dO1_dz1*dO1_dx2
     &                      +dO2_dz1*dO2_dx2
     &                      +dO3_dz1*dO3_dx2)
            Hess(LHR(3,iAtom,2,jAtom))=Hess(LHR(3,iAtom,2,jAtom))
     &          + f_const  *(dO1_dz1*dO1_dy2
     &                      +dO2_dz1*dO2_dy2
     &                      +dO3_dz1*dO3_dy2)
            Hess(LHR(3,iAtom,3,jAtom))=Hess(LHR(3,iAtom,3,jAtom))
     &          + f_const  *(dO1_dz1*dO1_dz2
     &                      +dO2_dz1*dO2_dz2
     &                      +dO3_dz1*dO3_dz2)
*
         End Do
*
         Hess(LHR(1,iAtom,1,iAtom))=Hess(LHR(1,iAtom,1,iAtom))
     &          + f_const  *(dO1_dx1*dO1_dx1
     &                      +dO2_dx1*dO2_dx1
     &                      +dO3_dx1*dO3_dx1)
         Hess(LHR(2,iAtom,1,iAtom))=Hess(LHR(2,iAtom,1,iAtom))
     &          + f_const  *(dO1_dy1*dO1_dx1
     &                      +dO2_dy1*dO2_dx1
     &                      +dO3_dy1*dO3_dx1)
         Hess(LHR(2,iAtom,2,iAtom))=Hess(LHR(2,iAtom,2,iAtom))
     &          + f_const  *(dO1_dy1*dO1_dy1
     &                      +dO2_dy1*dO2_dy1
     &                      +dO3_dy1*dO3_dy1)
         Hess(LHR(3,iAtom,1,iAtom))=Hess(LHR(3,iAtom,1,iAtom))
     &          + f_const  *(dO1_dz1*dO1_dx1
     &                      +dO2_dz1*dO2_dx1
     &                      +dO3_dz1*dO3_dx1)
         Hess(LHR(3,iAtom,2,iAtom))=Hess(LHR(3,iAtom,2,iAtom))
     &          + f_const  *(dO1_dz1*dO1_dy1
     &                      +dO2_dz1*dO2_dy1
     &                      +dO3_dz1*dO3_dy1)
         Hess(LHR(3,iAtom,3,iAtom))=Hess(LHR(3,iAtom,3,iAtom))
     &          + f_const  *(dO1_dz1*dO1_dz1
     &                      +dO2_dz1*dO2_dz1
     &                      +dO3_dz1*dO3_dz1)
*
      End Do
      Call Free_Work(ip_Grad)
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian after Rotation','(12f12.7)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
 778  Continue
      Call Free_Work(ip_xMass)
*                                                                      *
************************************************************************
*                                                                      *
*     Hessian for tension
*
*
      Do iBond = 1, nBonds
         kAtom = iTabBonds(1,iBond)
         lAtom = iTabBonds(2,iBond)
         iBondType  = iTabBonds(3,iBond)
C        If (iBondType.gt.Magic_Bond) Go To 10
         kr=iTabRow(iANr(kAtom))
         lr=iTabRow(iANr(lAtom))
         Help = kr.gt.3.or.lr.gt.3
         xkl=Cart(1,kAtom)-Cart(1,lAtom)
         ykl=Cart(2,kAtom)-Cart(2,lAtom)
         zkl=Cart(3,kAtom)-Cart(3,lAtom)
         rkl2 = xkl**2 + ykl**2 + zkl**2
         r0=rAv(kr,lr)
         alpha=aAv(kr,lr)
*
         If (Schlegel.or.Help) Then
            Rab=Sqrt(rkl2)
            RabCov=CovRad(iANr(kAtom))+CovRad(iANr(lAtom))
            If ((kr.eq.1.and.lr.eq.1).or.Help) Then
               gmm=Fact*A_StrH(1)*Exp(-A_StrH(2)*(Rab-RabCov))
            Else
               ij=Max(kr,lr)*(Max(kr,lr)-1)/2 + Min(kr,lr)
               gmm=Fact*A_Str/(Rab-B_Str(ij))**3
            End If
         Else
            gmm=rkr     *Exp(alpha    *(r0    **2-rkl2))
            If (iAnd(iOptC,1024).eq.1024) Then
               r0_vdW= r_ref_vdW(kr,lr)
               g_vdW = rkr_vdW*Exp(-alpha_vdW*(r0_vdW-SQRT(rkl2))**2)
            Else
               g_vdW=0.0D0
            End If
            gmm = gmm + g_vdW
         End If
*
         f_const = Max(gmm,f_const_Min_)
#ifdef _DEBUG_
         nqR=nqR+1
         Write (6,*) 'ddV: bonds: kAtom,lAtom=',kAtom,LAtom
         Write (6,*) '          : Bondtype=',iBondType
*        Write (6,*) gmm/rkr, f_const, g_vdW
         Write (6,*) f_const
#endif
         Hxx=f_const*xkl*xkl/rkl2
         Hxy=f_const*xkl*ykl/rkl2
         Hxz=f_const*xkl*zkl/rkl2
         Hyy=f_const*ykl*ykl/rkl2
         Hyz=f_const*ykl*zkl/rkl2
         Hzz=f_const*zkl*zkl/rkl2
*
         Hess(LHR(1,kAtom,1,kAtom))=Hess(LHR(1,kAtom,1,kAtom))+Hxx
         Hess(LHR(2,kAtom,1,kAtom))=Hess(LHR(2,kAtom,1,kAtom))+Hxy
         Hess(LHR(2,kAtom,2,kAtom))=Hess(LHR(2,kAtom,2,kAtom))+Hyy
         Hess(LHR(3,kAtom,1,kAtom))=Hess(LHR(3,kAtom,1,kAtom))+Hxz
         Hess(LHR(3,kAtom,2,kAtom))=Hess(LHR(3,kAtom,2,kAtom))+Hyz
         Hess(LHR(3,kAtom,3,kAtom))=Hess(LHR(3,kAtom,3,kAtom))+Hzz
*
         Hess(LHR(1,kAtom,1,lAtom))=Hess(LHR(1,kAtom,1,lAtom))-Hxx
         Hess(LHR(1,kAtom,2,lAtom))=Hess(LHR(1,kAtom,2,lAtom))-Hxy
         Hess(LHR(1,kAtom,3,lAtom))=Hess(LHR(1,kAtom,3,lAtom))-Hxz
         Hess(LHR(2,kAtom,1,lAtom))=Hess(LHR(2,kAtom,1,lAtom))-Hxy
         Hess(LHR(2,kAtom,2,lAtom))=Hess(LHR(2,kAtom,2,lAtom))-Hyy
         Hess(LHR(2,kAtom,3,lAtom))=Hess(LHR(2,kAtom,3,lAtom))-Hyz
         Hess(LHR(3,kAtom,1,lAtom))=Hess(LHR(3,kAtom,1,lAtom))-Hxz
         Hess(LHR(3,kAtom,2,lAtom))=Hess(LHR(3,kAtom,2,lAtom))-Hyz
         Hess(LHR(3,kAtom,3,lAtom))=Hess(LHR(3,kAtom,3,lAtom))-Hzz
*
         Hess(LHR(1,lAtom,1,lAtom))=Hess(LHR(1,lAtom,1,lAtom))+Hxx
         Hess(LHR(2,lAtom,1,lAtom))=Hess(LHR(2,lAtom,1,lAtom))+Hxy
         Hess(LHR(2,lAtom,2,lAtom))=Hess(LHR(2,lAtom,2,lAtom))+Hyy
         Hess(LHR(3,lAtom,1,lAtom))=Hess(LHR(3,lAtom,1,lAtom))+Hxz
         Hess(LHR(3,lAtom,2,lAtom))=Hess(LHR(3,lAtom,2,lAtom))+Hyz
         Hess(LHR(3,lAtom,3,lAtom))=Hess(LHR(3,lAtom,3,lAtom))+Hzz
*
C10      Continue
      End Do
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian after tension','(12f12.7)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Hessian for bending
*
      If (nBonds.lt.2) Go To 999
      Do mAtom = 1, nAtoms
         mr=iTabRow(iANr(mAtom))
*
         nNeighbor=iTabAtoms(1,0,mAtom)
         If (nNeighbor.lt.2) Go To 20
         Do iNeighbor = 1, nNeighbor
            iAtom=iTabAtoms(1,iNeighbor,mAtom)
            iBond=iTabAtoms(2,iNeighbor,mAtom)
            iBondType  = iTabBonds(3,iBond)
            If (iBondType.gt.Magic_Bond) Go To 30
            ir=iTabRow(iANr(iAtom))
*
            xmi=(Cart(1,iAtom)-Cart(1,mAtom))
            ymi=(Cart(2,iAtom)-Cart(2,mAtom))
            zmi=(Cart(3,iAtom)-Cart(3,mAtom))
            rmi2 = xmi**2 + ymi**2 + zmi**2
            rmi=sqrt(rmi2)
            r0mi=rAv(mr,ir)
            ami=aAv(mr,ir)
*
            Do jNeighbor = 1, iNeighbor-1
               jAtom=iTabAtoms(1,jNeighbor,mAtom)
               jBond=iTabAtoms(2,jNeighbor,mAtom)
               jBondType  = iTabBonds(3,jBond)
               If (jBondType.gt.Magic_Bond) Go To 40
               jr=iTabRow(iANr(jAtom))
               Help=mr.gt.3.or.ir.gt.3.or.jr.gt.3
*
               xmj=(Cart(1,jAtom)-Cart(1,mAtom))
               ymj=(Cart(2,jAtom)-Cart(2,mAtom))
               zmj=(Cart(3,jAtom)-Cart(3,mAtom))
               rmj2 = xmj**2 + ymj**2 + zmj**2
               rmj=sqrt(rmj2)
               r0mj=rAv(mr,jr)
               amj=aAv(mr,jr)
               r=sqrt(rmj2+rmi2)
*
*------------- Test if zero angle
*
               Test_zero=xmi*xmj+ymi*ymj+zmi*zmj
               Test_zero=Test_zero/(rmi*rmj)
               If (Abs(Test_zero-One).lt.1.0D-12) Go To 40
*
               xij=(Cart(1,jAtom)-Cart(1,iAtom))
               yij=(Cart(2,jAtom)-Cart(2,iAtom))
               zij=(Cart(3,jAtom)-Cart(3,iAtom))
               rij2 = xij**2 + yij**2 + zij**2
               rrij=sqrt(rij2)
*
               If (Schlegel.or.Help) Then
                  Rab=rmi
                  RabCov=CovRad(iANr(iAtom))+CovRad(iANr(mAtom))
                  Rbc=rmj
                  RbcCov=CovRad(iANr(jAtom))+CovRad(iANr(mAtom))
                  If (ir.eq.1.or.jr.eq.1) Then
                     gij=Fact*A_Bend(1)
                  Else
                     gij=Fact*A_Bend(2)
                  End If
               Else
                  gim=exp(ami*(r0mi**2-rmi2))
                  gjm=exp(amj*(r0mj**2-rmj2))
                  If (iAnd(iOptC,1024).eq.1024) Then
                     r0_vdW_im= r_ref_vdW(ir,mr)
                     g_vdW_im = Exp(-alpha_vdW*
     &                              (r0_vdW_im-SQRT(rmi2))**2
     &                             )
                     r0_vdW_jm= r_ref_vdW(jr,mr)
                     g_vdW_jm = Exp(-alpha_vdW*
     &                              (r0_vdW_jm-SQRT(rmj2))**2
     &                             )
                  Else
                     g_vdW_im=0.0D0
                     g_vdW_jm=0.0D0
                  End If
                  g_vdW_im = g_vdW_im * rkr_vdW / rkr
                  g_vdW_jm = g_vdW_jm * rkr_vdW / rkr
                  gim = gim + Half*g_vdW_im
                  gjm = gjm + Half*g_vdW_jm
                  gij = rkf* gim*gjm
               End If
               rL2=(ymi*zmj-zmi*ymj)**2
     &            +(zmi*xmj-xmi*zmj)**2
     &            +(xmi*ymj-ymi*xmj)**2
chjw modified
            If(rL2.lt.1.d-14) then
              rL=Zero
            else
              rL=sqrt(rL2)
            end if
            gij=Max(gij,f_const_Min_)
#ifdef _DEBUG_
            Write (6,*) 'iAtom,mAtom,jAtom=',iAtom,mAtom,jAtom
            Write (6,*) 'gij=',gij
            Write (6,*) 'rmj=',rmj
            Write (6,*) 'rmi=',rmi
            Write (6,*) 'rrij=',rrij
#endif
*
            if ((rmj.gt.rZero).and.(rmi.gt.rZero).and.
     &                                (rrij.gt.rZero)) Then
#ifdef _DEBUG_
              nqB=nqB+1
#endif
              SinPhi=rL/(rmj*rmi)
              rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
              CosPhi=rmidotrmj/(rmj*rmi)
*
*-------------Non-linear case
*
              Thr_Line=Sin(Pi*25.0d0/180.0D0)
              If (nAtoms.eq.3) Thr_Line=rZero
              If (SinPhi.gt.Thr_Line) Then
                si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                sm(1)=-si(1)-sj(1)
                sm(2)=-si(2)-sj(2)
                sm(3)=-si(3)-sj(3)
                Do icoor=1,3
                   Do jCoor=1,3
                    If (mAtom.gt.iAtom) Then
                       Hess(LHR(icoor,mAtom,jcoor,iAtom))=
     &                        Hess(LHR(icoor,mAtom,jcoor,iAtom))
     &                        +gij*sm(icoor)*si(jcoor)
                    else
                      Hess(LHR(icoor,iAtom,jcoor,mAtom))=
     &                        Hess(LHR(icoor,iAtom,jcoor,mAtom))
     &                        +gij*si(icoor)*sm(jcoor)
                    End If
                    If (mAtom.gt.jAtom) Then
                        Hess(LHR(icoor,mAtom,jcoor,jAtom))=
     &                        Hess(LHR(icoor,mAtom,jcoor,jAtom))
     &                        +gij*sm(icoor)*sj(jcoor)
                    else
                      Hess(LHR(icoor,jAtom,jcoor,mAtom))=
     &                        Hess(LHR(icoor,jAtom,jcoor,mAtom))
     &                        +gij*sj(icoor)*sm(jcoor)
                    End If
                    If (iAtom.gt.jAtom) Then
                        Hess(LHR(icoor,iAtom,jcoor,jAtom))=
     &                        Hess(LHR(icoor,iAtom,jcoor,jAtom))
     &                        +gij*si(icoor)*sj(jcoor)
                     else
                        Hess(LHR(icoor,jAtom,jcoor,iAtom))=
     &                        Hess(LHR(icoor,jAtom,jcoor,iAtom))
     &                        +gij*sj(icoor)*si(jcoor)
                     End If
                   End Do
                End Do
                Do icoor=1,3
                  Do jCoor=1,icoor
                    Hess(LHR(icoor,iAtom,jcoor,iAtom))=
     &                        Hess(LHR(icoor,iAtom,jcoor,iAtom))
     &                        +gij*si(icoor)*si(jcoor)
                    Hess(LHR(icoor,mAtom,jcoor,mAtom))=
     &                        Hess(LHR(icoor,mAtom,jcoor,mAtom))
     &                        +gij*sm(icoor)*sm(jcoor)
                    Hess(LHR(icoor,jAtom,jcoor,jAtom))=
     &                        Hess(LHR(icoor,jAtom,jcoor,jAtom))
     &                        +gij*sj(icoor)*sj(jcoor)

*
                  End Do
                End Do
              Else
*
*---------------Linear case
*
                If ((abs(ymi).gt.rZero).or.
     &              (abs(xmi).gt.rZero)) Then
                   x(1)=-ymi
                   y(1)=xmi
                   z(1)=Zero
                   x(2)=-xmi*zmi
                   y(2)=-ymi*zmi
                   z(2)=xmi*xmi+ymi*ymi
                Else
                   x(1)=One
                   y(1)=Zero
                   z(1)=Zero
                   x(2)=Zero
                   y(2)=One
                   z(2)=Zero
                End If
#ifdef _DEBUG_
                nqB=nqB+2
#endif
                Do i=1,2
                   r1=sqrt(x(i)**2+y(i)**2+z(i)**2)
                   cosThetax=x(i)/r1
                   cosThetay=y(i)/r1
                   cosThetaz=z(i)/r1
                   si(1)=-cosThetax/rmi
                   si(2)=-cosThetay/rmi
                   si(3)=-cosThetaz/rmi
                   sj(1)=-cosThetax/rmj
                   sj(2)=-cosThetay/rmj
                   sj(3)=-cosThetaz/rmj
                   sm(1)=-(si(1)+sj(1))
                   sm(2)=-(si(2)+sj(2))
                   sm(3)=-(si(3)+sj(3))
*
                   Do icoor=1,3
                      Do jCoor=1,3
                         If (mAtom.gt.iAtom) Then
                          Hess(LHR(icoor,mAtom,jcoor,iAtom))=
     &                        Hess(LHR(icoor,mAtom,jcoor,iAtom))
     &                         +gij*sm(icoor)*si(jcoor)
                        else
                           Hess(LHR(icoor,iAtom,jcoor,mAtom))=
     &                        Hess(LHR(icoor,iAtom,jcoor,mAtom))
     &                         +gij*si(icoor)*sm(jcoor)
                        End If
                        If (mAtom.gt.jAtom) Then
                          Hess(LHR(icoor,mAtom,jcoor,jAtom))=
     &                        Hess(LHR(icoor,mAtom,jcoor,jAtom))
     &                         +gij*sm(icoor)*sj(jcoor)
                        else
                          Hess(LHR(icoor,jAtom,jcoor,mAtom))=
     &                        Hess(LHR(icoor,jAtom,jcoor,mAtom))
     &                         +gij*sj(icoor)*sm(jcoor)
                        End If
                        If (iAtom.gt.jAtom) Then
                           Hess(LHR(icoor,iAtom,jcoor,jAtom))=
     &                        Hess(LHR(icoor,iAtom,jcoor,jAtom))
     &                         +gij*si(icoor)*sj(jcoor)
                        else
                           Hess(LHR(icoor,jAtom,jcoor,iAtom))=
     &                        Hess(LHR(icoor,jAtom,jcoor,iAtom))
     &                         +gij*sj(icoor)*si(jcoor)
                        End If
                       End Do
                     End Do
                     Do icoor=1,3
                       Do jCoor=1,icoor
                         Hess(LHR(icoor,iAtom,jcoor,iAtom))=
     &                        Hess(LHR(icoor,iAtom,jcoor,iAtom))
     &                         +gij*si(icoor)*si(jcoor)
                         Hess(LHR(icoor,mAtom,jcoor,mAtom))=
     &                        Hess(LHR(icoor,mAtom,jcoor,mAtom))
     &                         +gij*sm(icoor)*sm(jcoor)
                         Hess(LHR(icoor,jAtom,jcoor,jAtom))=
     &                         Hess(LHR(icoor,jAtom,jcoor,jAtom))
     &                         +gij*sj(icoor)*sj(jcoor)
                       End Do
                     End Do
                 End Do
              End If
            End If
*
 40         Continue
        End Do
 30         Continue
       End Do
 20    Continue
      End Do
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian after bending','(12f12.7)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Hessian for torsion
*
      If (nBonds.lt.3) Go To 999
      Do iBond = 1, nBonds
         jAtom=iTabBonds(1,iBond)
         kAtom=iTabBonds(2,iBond)
         iBondType =iTabBonds(3,iBond)
         Fact = One
         If (iBondType.gt.Magic_Bond) Fact=Two
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) '*',jAtom,kAtom,' *'
         Write (6,*)
         Write (6,*) 'BondType=',iBondType
#endif
*
*        Allow center bond to be a "magic" bond
*
C        If (iBondType.eq.vdW_Bond) Go To 444

         jr=iTabRow(iANr(jAtom))
         kr=iTabRow(iANr(kAtom))
*
         call dcopy_(3,Cart(1,jAtom),1,xyz(1,2),1)
         call dcopy_(3,Cart(1,kAtom),1,xyz(1,3),1)
*
         nNeighbor_j = iTabAtoms(1,0,jAtom)
         If (nNeighbor_j.lt.2) Go To 444
         nNeighbor_k = iTabAtoms(1,0,kAtom)
         If (nNeighbor_k.lt.2) Go To 444
*
         Do jNeighbor = 1, nNeighbor_j
            iAtom = iTabAtoms(1,jNeighbor,jAtom)
#ifdef _DEBUG_
*           Write (6,*)
*           Write (6,*) iAtom,jAtom,kAtom,' *'
*           Write (6,*)
#endif
            jBond = iTabAtoms(2,jNeighbor,jAtom)
            If (iBond.eq.jBond) Go To 333
            jBondType =iTabBonds(3,jBond)
C           If (jBondType.eq.vdW_Bond) Go To 333
            If (jBondType.gt.Magic_Bond) Go To 333
            ir=iTabRow(iANr(iAtom))
*
            call dcopy_(3,Cart(1,iAtom),1,xyz(1,1),1)
*
            Do kNeighbor = 1, nNeighbor_k
               lAtom = iTabAtoms(1,kNeighbor,kAtom)
               kBond = iTabAtoms(2,kNeighbor,kAtom)
               If (iBond.eq.kBond) Go To 222
               If (lAtom.eq.iAtom) Go To 222
               kBondType =iTabBonds(3,kBond)
C              If (kBondType.eq.vdW_Bond) Go To 222
               If (kBondType.gt.Magic_Bond) Go To 222
               lr=iTabRow(iANr(lAtom))
               Help=kr.gt.3.or.ir.gt.3.or.jr.gt.3.or.lr.gt.3
*
               call dcopy_(3,Cart(1,lAtom),1,xyz(1,4),1)
*
               rij(1)=Cart(1,iAtom)-Cart(1,jAtom)
               rij(2)=Cart(2,iAtom)-Cart(2,jAtom)
               rij(3)=Cart(3,iAtom)-Cart(3,jAtom)
               rij2=rij(1)**2+rij(2)**2+rij(3)**2
               rij0=rAv(ir,jr)**2
               aij =aAv(ir,jr)
*
               rjk(1)=Cart(1,jAtom)-Cart(1,kAtom)
               rjk(2)=Cart(2,jAtom)-Cart(2,kAtom)
               rjk(3)=Cart(3,jAtom)-Cart(3,kAtom)
               rjk2=rjk(1)**2+rjk(2)**2+rjk(3)**2
               rjk0=rAv(jr,kr)**2
               ajk =aAv(jr,kr)
*
               rkl(1)=Cart(1,kAtom)-Cart(1,lAtom)
               rkl(2)=Cart(2,kAtom)-Cart(2,lAtom)
               rkl(3)=Cart(3,kAtom)-Cart(3,lAtom)
               rkl2=rkl(1)**2+rkl(2)**2+rkl(3)**2
               rkl0=rAv(kr,lr)**2
               akl =aAv(kr,lr)
*
*              Allow only angles in the range of 35-145
               A35 = (35.0D0/180.D0)* Pi
               CosFi_Max=Cos(A35)
               CosFi2=(rij(1)*rjk(1)+rij(2)*rjk(2)+rij(3)*rjk(3))
     &               /Sqrt(rij2*rjk2)
               If (Abs(CosFi2).gt.CosFi_Max) Go To 222
               CosFi3=(rkl(1)*rjk(1)+rkl(2)*rjk(2)+rkl(3)*rjk(3))
     &               /Sqrt(rkl2*rjk2)
               If (Abs(CosFi3).gt.CosFi_Max) Go To 222
#ifdef _DEBUG_
               Write (6,*) 'CosFi2,CosFi3=',CosFi2,CosFi3
               Write (6,*) 'rij=',rij,rij2
               Write (6,*) 'rjk=',rjk,rjk2
               Write (6,*) 'rkl=',rkl,rkl2
#endif
*
               If (Schlegel.or.Help) Then
                  Rab=Sqrt(rij2)
                  RabCov=(CovRadT(iANr(iAtom))
     &                   +CovRadT(iANr(jAtom)))/bohr
                  Rbc=Sqrt(rjk2)/Fact
                  RbcCov=(CovRadT(iANr(jAtom))
     &                   +CovRadT(iANr(kAtom)))/bohr
                  Rcd=Sqrt(rkl2)
                  RcdCov=(CovRadT(iANr(kAtom))
     &                   +CovRadT(iANr(lAtom)))/bohr
                  Diff=RbcCov-Rbc
                  If (Diff.lt.Zero) Diff=Zero
                  tij=Fact*A_Trsn(1)+A_Trsn(2)*Diff
               Else
*                 Magic bond fix
                  rjk2=rjk2/Fact**2
*
                  g_ij=exp(aij*(rij0-rij2))
                  g_jk=exp(ajk*(rjk0-rjk2))
                  g_kl=exp(akl*(rkl0-rkl2))
                  If (iAnd(iOptC,1024).eq.1024) Then
                     r0_vdW_ij= r_ref_vdW(ir,jr)
                     g_vdW_ij = Exp(-alpha_vdW*
     &                              (r0_vdW_ij-SQRT(rij2))**2
     &                             )
                     r0_vdW_jk= r_ref_vdW(jr,kr)
                     g_vdW_jk = Exp(-alpha_vdW*
     &                              (r0_vdW_jk-SQRT(rjk2))**2
     &                             )
                     r0_vdW_kl= r_ref_vdW(kr,lr)
                     g_vdW_kl = Exp(-alpha_vdW*
     &                              (r0_vdW_kl-SQRT(rkl2))**2
     &                             )
                  Else
                     g_vdW_ij=0.0D0
                     g_vdW_jk=0.0D0
                     g_vdW_kl=0.0D0
                  End If
                  g_vdW_ij = g_vdW_ij * rkr_vdW / rkr
                  g_vdW_jk = g_vdW_jk * rkr_vdW / rkr
                  g_vdW_kl = g_vdW_kl * rkr_vdW / rkr
                  g_ij=g_ij + Half * g_vdW_ij
                  g_jk=g_jk + Half * g_vdW_jk
                  g_kl=g_kl + Half * g_vdW_kl
                  tij=rkt*g_ij*g_jk*g_kl
               End If
               tij = Max(tij,f_const_Min_)
               If (Torsion_Check(iAtom,jAtom,kAtom,lAtom,
     &                           xyz,iTabAtoms,
     &                           nMax,nAtoms) )Then
                  tij = Max(tij,10.0D0*f_const_Min_)
               End If
#ifdef _DEBUG_
               nqT=nqT+1
               Write (6,*)
               Write (6,*) iAtom,jAtom,kAtom,lAtom
               Write (6,*) tij
#endif
*
               Call Trsn(xyz,4,Tau,C,.False.,.False.,'        ',
     &                   Dum,.False.)
               call dcopy_(3,C(1,1),1,si,1)
               call dcopy_(3,C(1,2),1,sj,1)
               call dcopy_(3,C(1,3),1,sk,1)
               call dcopy_(3,C(1,4),1,sl,1)
#ifdef _DEBUG_
*              Call RecPrt('C',' ',C,3,4)
#endif
*
*------------- Off diagonal block
*
              Do icoor=1,3
                Do jCoor=1,3
                 Hess(LHR(icoor,iAtom,jcoor,jAtom))=
     &           Hess(LHR(icoor,iAtom,jcoor,jAtom))
     &            +tij*si(icoor) * sj(jcoor)
                 Hess(LHR(icoor,iAtom,jcoor,kAtom))=
     &           Hess(LHR(icoor,iAtom,jcoor,kAtom))
     &            +tij*si(icoor) * sk(jcoor)
                 Hess(LHR(icoor,iAtom,jcoor,lAtom))=
     &           Hess(LHR(icoor,iAtom,jcoor,lAtom))
     &            +tij*si(icoor) * sl(jcoor)
                 Hess(LHR(icoor,jAtom,jcoor,kAtom))=
     &           Hess(LHR(icoor,jAtom,jcoor,kAtom))
     &            +tij*sj(icoor) * sk(jcoor)
                 Hess(LHR(icoor,jAtom,jcoor,lAtom))=
     &           Hess(LHR(icoor,jAtom,jcoor,lAtom))
     &            +tij*sj(icoor) * sl(jcoor)
                 Hess(LHR(icoor,kAtom,jcoor,lAtom))=
     &           Hess(LHR(icoor,kAtom,jcoor,lAtom))
     &            +tij*sk(icoor) * sl(jcoor)

                End Do
              End Do
*
*-------------Diagonal block
*
              Do icoor=1,3
                Do jCoor=1,icoor
                 Hess(LHR(icoor,iAtom,jcoor,iAtom))=
     &           Hess(LHR(icoor,iAtom,jcoor,iAtom))
     &            +tij*si(icoor) * si(jcoor)
                 Hess(LHR(icoor,jAtom,jcoor,jAtom))=
     &           Hess(LHR(icoor,jAtom,jcoor,jAtom))
     &            +tij*sj(icoor) * sj(jcoor)
                 Hess(LHR(icoor,kAtom,jcoor,kAtom))=
     &           Hess(LHR(icoor,kAtom,jcoor,kAtom))
     &            +tij*sk(icoor) * sk(jcoor)
                 Hess(LHR(icoor,lAtom,jcoor,lAtom))=
     &           Hess(LHR(icoor,lAtom,jcoor,lAtom))
     &            +tij*sl(icoor) * sl(jcoor)

*
                 End Do
               End Do
*
 222           Continue
            End Do          ! iNeigbor_k
 333        Continue
         End Do             ! iNeighbor_j
 444     Continue
      End Do                ! iBonds
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian after torsion','(12f12.7)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*                                      k
*                                     /
*     Hessian for out-of-plane   j - i
*                                     \
*                                      l
*
C     Go To 999
      If (nBonds.lt.3) Go To 999
*
      Do iAtom = 1, nAtoms
*
         nNeighbor_i = iTabAtoms(1,0,iAtom)
C        Write (*,*) 'iAtom,nNeighbor_i=',iAtom,nNeighbor_i
         If (nNeighbor_i.lt.3) Go To 446
         ir=iTabRow(iANr(iAtom))
         call dcopy_(3,Cart(1,iAtom),1,xyz(1,4),1)
*
         Do iNb0 = 1, nNeighbor_i
            jAtom = iTabAtoms(1,iNb0,iAtom)
C           Write (*,*) 'jAtom=',jAtom
            jr=iTabRow(iANr(jAtom))
            iBond = iTabAtoms(2,iNb0,iAtom)
            iBondType=iTabBonds(3,iBond)
C           Write (*,*) 'iBondType=',iBondType
            nCoBond_j=nCoBond(jAtom,nAtoms,nMax,iTabBonds,nBonds,
     &                        nBonds,iTabAtoms)
            If (nCoBond_j.gt.1) Go To 447
*           If (iBondType.eq.vdW_Bond) Go To 447
            If (iBondType.gt.Magic_Bond) Go To 447
            call dcopy_(3,Cart(1,jAtom),1,xyz(1,1),1)
*
            Do iNb1 = 1, nNeighbor_i
               kAtom = iTabAtoms(1,iNb1,iAtom)
C              Write (*,*) 'kAtom=',kAtom
               kBond = iTabAtoms(2,iNb1,iAtom)
               If (kAtom.eq.jAtom) Go To 335
               kBondType=iTabBonds(3,kBond)
C              Write (*,*) 'kBondType=',kBondType
*              If (kBondType.eq.vdW_Bond) Go To 335
               If (kBondType.gt.Magic_Bond) Go To 335
               kr=iTabRow(iANr(kAtom))
*
               call dcopy_(3,Cart(1,kAtom),1,xyz(1,2),1)
*
               Do iNb2 = 1, nNeighbor_i
                  lAtom = iTabAtoms(1,iNb2,iAtom)
C                 Write (6,*) 'lAtom=',lAtom
                  lBond = iTabAtoms(2,iNb2,iAtom)

                  If (lAtom.eq.jAtom) Go To 224
                  If (lAtom.le.kAtom) Go To 224
                  lBondType=iTabBonds(3,lBond)
C                 Write (*,*) 'lBondType=',lBondType
*                 If (lBondType.eq.vdW_Bond) Go To 224
                  If (lBondType.gt.Magic_Bond) Go To 224
                  lr=iTabRow(iANr(lAtom))
                  Help=kr.gt.3.or.ir.gt.3.or.jr.gt.3.or.lr.gt.3
*
C                 Write (*,*) 'i,j,k,l=',iAtom,jAtom,kAtom,lAtom
C                 Write (*,*) 'Help=',Help
*
                  call dcopy_(3,Cart(1,lAtom),1,xyz(1,3),1)
*
                  rij(1)=Cart(1,iAtom)-Cart(1,jAtom)
                  rij(2)=Cart(2,iAtom)-Cart(2,jAtom)
                  rij(3)=Cart(3,iAtom)-Cart(3,jAtom)
                  rij0=rAv(ir,jr)**2
                  aij =aAv(ir,jr)
*
                  rik(1)=Cart(1,iAtom)-Cart(1,kAtom)
                  rik(2)=Cart(2,iAtom)-Cart(2,kAtom)
                  rik(3)=Cart(3,iAtom)-Cart(3,kAtom)
                  rik0=rAv(ir,kr)**2
                  aik =aAv(ir,kr)
*
                  ril(1)=Cart(1,iAtom)-Cart(1,lAtom)
                  ril(2)=Cart(2,iAtom)-Cart(2,lAtom)
                  ril(3)=Cart(3,iAtom)-Cart(3,lAtom)
                  ril0=rAv(ir,lr)**2
                  ail =aAv(ir,lr)
*
                  rij2=rij(1)**2+rij(2)**2+rij(3)**2
                  rik2=rik(1)**2+rik(2)**2+rik(3)**2
                  ril2=ril(1)**2+ril(2)**2+ril(3)**2
*
                  ThrFi1=Cos(90.0D0*Pi/(180.0D0))
                  ThrFi2=Cos(150.0D0*Pi/(180.0D0))
                  CosFi2=(rij(1)*rik(1)+rij(2)*rik(2)+rij(3)*rik(3))
     &                  /Sqrt(rij2*rik2)
                  If (CosFi2.gt.ThrFi1 .or. CosFi2.lt.ThrFi2) Go To 224
*
                  CosFi3=(rij(1)*ril(1)+rij(2)*ril(2)+rij(3)*ril(3))
     &                  /Sqrt(rij2*ril2)
                  If (CosFi3.gt.ThrFi1 .or. CosFi3.lt.ThrFi2) Go To 224
*
                  CosFi4=(rik(1)*ril(1)+rik(2)*ril(2)+rik(3)*ril(3))
     &                  /Sqrt(rik2*ril2)
                  If (CosFi4.gt.ThrFi1 .or. CosFi4.lt.ThrFi2) Go To 224
#ifdef _DEBUG_
                  Write (6,*) 'CosFi2,CosFi3,CosFi4=',
     &                        CosFi2,CosFi3,CosFi4
#endif
*
                  If (Schlegel.or.Help) Then
*
*------------------- I do not have a clue to how this will really work!
*
                     tij = f_const_Min_
                 Else
                     beta=rko*
     &                        exp( (aij*rij0+aik*rik0+ail*ril0))
                     tij=beta*exp(-(aij*rij2+aik*rik2+ail*ril2))
                  End If
C                 tij=Max(tij,f_const_Min_)
*
                  Call OutofP(xyz,4,Tau,C,.False.,.False.,'        ',
     &                        Dum,.False.)
                  If (Abs(Tau).gt.25.0D0*(Pi/180.D0)) Go To 224
#ifdef _DEBUG_
                  nqO=nqO+1
#endif
*
                  call dcopy_(3,C(1,4),1,si,1)
                  call dcopy_(3,C(1,1),1,sj,1)
                  call dcopy_(3,C(1,2),1,sk,1)
                  call dcopy_(3,C(1,3),1,sl,1)
#ifdef _DEBUG_
                  Write (6,*) 'iAtoms=',
     &                         iAtom,jAtom,kAtom,lAtom
                  Write(6,*) 'tij,Tau=',tij,Tau
*                 Call RecPrt('si',' ',si,1,3)
*                 Call RecPrt('sj',' ',sj,1,3)
*                 Call RecPrt('sk',' ',sk,1,3)
*                 Call RecPrt('sl',' ',sl,1,3)
#endif
*
*---------------- Off diagonal block
*
                  Do icoor=1,3
                    Do jCoor=1,3
                     Hess(LHR(icoor,iAtom,jcoor,jAtom))=
     &               Hess(LHR(icoor,iAtom,jcoor,jAtom))
     &                +tij*si(icoor) * sj(jcoor)
                     Hess(LHR(icoor,iAtom,jcoor,kAtom))=
     &               Hess(LHR(icoor,iAtom,jcoor,kAtom))
     &                +tij*si(icoor) * sk(jcoor)
                     Hess(LHR(icoor,iAtom,jcoor,lAtom))=
     &               Hess(LHR(icoor,iAtom,jcoor,lAtom))
     &                +tij*si(icoor) * sl(jcoor)
                     Hess(LHR(icoor,jAtom,jcoor,kAtom))=
     &               Hess(LHR(icoor,jAtom,jcoor,kAtom))
     &                +tij*sj(icoor) * sk(jcoor)
                     Hess(LHR(icoor,jAtom,jcoor,lAtom))=
     &               Hess(LHR(icoor,jAtom,jcoor,lAtom))
     &                +tij*sj(icoor) * sl(jcoor)
                     Hess(LHR(icoor,kAtom,jcoor,lAtom))=
     &               Hess(LHR(icoor,kAtom,jcoor,lAtom))
     &                +tij*sk(icoor) * sl(jcoor)

                    End Do
                  End Do
*
*---------------- Diagonal block
*
                  Do icoor=1,3
                    Do jCoor=1,icoor
                     Hess(LHR(icoor,iAtom,jcoor,iAtom))=
     &               Hess(LHR(icoor,iAtom,jcoor,iAtom))
     &                +tij*si(icoor) * si(jcoor)
                     Hess(LHR(icoor,jAtom,jcoor,jAtom))=
     &               Hess(LHR(icoor,jAtom,jcoor,jAtom))
     &                +tij*sj(icoor) * sj(jcoor)
                     Hess(LHR(icoor,kAtom,jcoor,kAtom))=
     &               Hess(LHR(icoor,kAtom,jcoor,kAtom))
     &                +tij*sk(icoor) * sk(jcoor)
                     Hess(LHR(icoor,lAtom,jcoor,lAtom))=
     &               Hess(LHR(icoor,lAtom,jcoor,lAtom))
     &                +tij*sl(icoor) * sl(jcoor)

*
                     End Do
                   End Do
*
 224            Continue
                End Do        ! iNb2
 335        Continue
            End Do          ! iNb1
 447        Continue
         End Do             ! iCase
 446     Continue
      End Do               ! iBond
#ifdef _DEBUG_
      Call TriPrt(' In LNM: Hessian after out-of-plane','(12f12.7)',
     &               Hess,n3)
      Call DiagMtrx_T(Hess,n3,iNeg)
#endif
*                                                                      *
************************************************************************
*                                                                      *
 999  Continue
#ifdef _DEBUG_
      Write (6,*) 'ddV: nqR, nqB, nqT, nqO=',nqR, nqB, nqT, nqO
#endif
      Call QExit('ddV_')
      Return
      End
#ifdef _DEBUG_
      Subroutine DiagMtrx_T(H,nH,iNeg)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      character*16 filnam
      Real*8 H(*)
      Logical Exist
*
      Call QEnter('DiagMtrx')
      Lu=6
      iRout=22
      iprint=nPrint(iRout)
*
      Call GetMem('EVal','Allo','Real',ipEVal,nH*(nH+1)/2)
      Call GetMem('EVec','Allo','Real',ipEVec,nH*nH)
*
*---- Copy elements for H
*
      call dcopy_(nH*(nH+1)/2,H,1,Work(ipEVal),1)
*
*---- Set up a unit matrix
*
      call dcopy_(nH*nH,[Zero],0,Work(ipEVec),1)
      call dcopy_(nH,[One],0,Work(ipEVec),nH+1)
*
*---- Compute eigenvalues and eigenvectors
*
      Call NIDiag_new(Work(ipEVal),Work(ipEVec),nH,nH,0)
      Call Jacord(Work(ipEVal),Work(ipEVec),nH,nH)
*
*---- Print out the result
*
      iNeg=0
      Do i = 1, nH
         If (Work(i*(i+1)/2+ipEVal-1).lt.Zero) iNeg=iNeg+1
      End Do
      IF (iprint.gt.5) THEN
        Write (Lu,*)
        Write (Lu,*) 'Eigenvalues of the Hessian'
        Write (Lu,*)
        Write (Lu,'(5G20.6)') (Work(i*(i+1)/2+ipEVal-1),i=1,nH)
      END IF
*
      call f_Inquire('SPCINX',Exist)
*
      If (Exist .AND. iprint.gt.5) Then
*
*        Read linear combinations from disc
*

         LuTmp=11
         filnam='SPCINX'
c         Open(luTmp,File=filnam,Form='unformatted',Status='unknown')
         Call molcas_binaryopen_vanilla(luTmp,filnam)
         ReWind (LuTmp)
*
         Read (LuTmp) nq,nQQ
*
         If (nQQ.eq.nH) Then
*
           Call GetMem('rK','Allo','Real',iprK,nq*nQQ)
           Call GetMem('qEVec','Allo','Real',ipqEVec,nq*nH)
           Call Print_qEVec(Work(ipEVec),nH,ipEVal,nq,
     &                      Work(iprK),Work(ipqEVec),LuTmp)
*
           Call GetMem('qEVec','Free','Real',ipqEVec,nq*nH)
           Call GetMem('rK','Free','Real',iprK,nq*nQQ)
*
         Else
*
           Write (Lu,*)
           Write (Lu,*) 'Eigenvectors of the Hessian'
           Write (Lu,*)
           Do i = 1, nH
              Write (Lu,'(10F10.5)')
     &              (Work((j-1)*nH+i+ipEVec-1),j=1,nH)
           End Do
         End If
*
         Close (LuTmp)
*
      ElSE IF (iprint.gt.5) THEN
         Write (Lu,*)
         Write (Lu,*) 'Eigenvectors of the Hessian'
         Write (Lu,*)
         Do i = 1, nH
            Write (Lu,'(10F10.5)') (Work((j-1)*nH+i+ipEVec-1),j=1,nH)
         End Do
*
      End If
*
      Call GetMem('EVec','Free','Real',ipEVec,nH*nH)
      Call GetMem('EVal','Free','Real',ipEVal,nH*(nH+1)/2)
*
      Call QExit('DiagMtrx')
      Return
      End
#endif
