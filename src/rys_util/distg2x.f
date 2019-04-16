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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine Distg2x(Hess,DAO,nt,mvec,FINAL,
     &                  IndGrd,
     &                  IfHss,IndHss,iStab,kOp,iChBas,MxFnc,
     &                  nIrrep,iChTbl,Tr,IfGr,ioper)
************************************************************************
*                                                                      *
* Object: trace the gradient of the ERI's with the second order        *
*         density matrix                                               *
*                                                                      *
* Called from: Rysg2                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DGeMV   (ESSL)                                          *
*              DCopy   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Anders Bernhardsson Dept. of Theoretical Chemistry,      *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
#include "real.fh"
      Real*8 g2(78),  DAO(*),Final(nt,mvec),
     &         Prmt(0:7),Hess(*)
      Logical IfHss(3,3,4,4),Tr(4),IfGr(4)
      Integer IndGrd(3,4,0:(nIrrep-1)),kOp(4),iStab(4),iChBas(MxFnc),
     &          IndHss(4,3,4,3,0:(nIrrep-1)),
     &          iChTbl(0:7,0:7),nop(4),iOper(0:nIrrep-1)
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
      ix(icn,icar,jcn,jcar)=(((icn-1)*3+icar)*((icn-1)*3+icar-1))/2+
     &                       (jcn-1)*3+jcar
*
c     iRout = 239
      iPrint = 49
c     rc=0
*     iQ = 0
c     Call qEnter('Distg2')
*

      call dcopy_(78,[0.0d0],0,g2,1)
      Call dGeMV_('T',nT,mvec,
     &           One,Final,nT,
CBS  &           PAO,1,                PAO was not declared above .. typo??
     &           DAO,1,
     &           Zero,g2,1)
      If (iPrint.ge.49) Call TriPrt('before Local Hessian ',
     &                               ' ',g2,12)

*-----Compute some of the contributions via the translational invariance
*
      Do 200 iCn = 1, 4
         Do 210 iCar = 1, 3
           Do 220 jCn = 1,iCn
             If (iCn.eq.jCn) Then
              iStop=iCar
             Else
              iStop=3
             End If
             Do 230 jCar = 1, iStop
               if (Tr(iCn).or.Tr(jCn)) Then      ! iCn = 3  or jCn = 3
                ij = ix(iCn,iCar,jCn,jCar)
C               write(6,*) 'ij = ',ij
                g2(ij)=zero
C               write(6,*) 'g2 ',g2(ij)
*------------------------------------------------------------*
*
*               Both derivatives by translation!
*
*------------------------------------------------------------*
                if (tr(iCn).and.tr(jCn))      !   iCn = jCn = 3
     &          Then
C               write(6,*) icn,jcn, ' by translation'
                 Do 240 kCn = 1, 4
                  Do 250 lCn = 1,kCn
                     If (lCn.eq.kCn) then
                      iCa2=Min(iCar,jCar)
                      iCa1=Max(iCar,jCar)
                      If (IfHss(iCa1,iCa2,kCn,lCn)) Then
                        k1=Ix(kCn,iCa1,lCn,iCa2)
                        g2(ij)=g2(ij)+g2(k1)
C                       write(6,*) 'contribution from ',k1,g2(k1)
                      End If
                     Else
                      If (IfHss(iCar,jCar,kcn,lcn)) Then
                        k1=Ix(kCn,iCar,lCn,jCar)
                        k2=Ix(kCn,jCar,lCn,iCar)
                        g2(ij)=g2(ij)+g2(k1)+
     &                                    g2(k2)
C                       write(6,*) 'contribution from ',
C    &                  k1,k2,g2(k1),g2(k2)
                      End If
                     End If
 250              Continue
 240             Continue
*------------------------------------------------------------*
*
*               Centre jCn by translation
*
*------------------------------------------------------------*
                Else If (ifgr(iCn).and.tr(jCn))
     &          Then
C               write(6,*) jcn, ' by translation'
                  Do 260 kCn = 1, 4
                    If (kCn.gt.iCn) Then
                      iCn1=kCn
                      iCn2=iCn
                      iCa1=jCar
                      iCa2=iCar
                     Else If (kCn.lt.iCn) Then
                      iCn1=iCn
                      iCn2=kCn
                      iCa1=iCar
                      iCa2=jCar
                     Else
                      iCn1=iCn
                      iCn2=kCn
                      iCa1=Max(iCar,jCar)
                      iCa2=Min(iCar,jCar)
                     End If
                    If (IfHss(iCa1,iCa2,icn1,icn2))
     &              Then
                      kl=Ix(iCn1,iCa1,iCn2,iCa2)
                      g2(ij)=g2(ij)-g2(kl)
C                       write(6,*) 'contribution from ',kl,-g2(kl)
                    End If
 260              Continue
*------------------------------------------------------------*
*
*               Centre iCn by translation
*
*------------------------------------------------------------*
                Else If (IfGr(jCn).and.tr(iCn))
     &          Then
C               write(6,*) icn, ' by translation'
                  Do 270 kCn = 1, 4
                    If (kCn.gt.jCn) Then
                      iCn1=kCn
                      iCn2=jCn
                      iCa1=iCar
                      iCa2=jCar
                     Else If (kCn.lt.jCn) Then
                      iCn1=jCn
                      iCn2=kCn
                      iCa1=jCar
                      iCa2=iCar
                     Else
                      iCn1=jCn
                      iCn2=kCn
                      iCa1=Max(iCar,jCar)
                      iCa2=Min(iCar,jCar)
                     End If
                    If (IfHss(iCa1,iCa2,icn1,icn2))
     &              Then
                      kl=Ix(iCn1,iCa1,iCn2,iCa2)
                      g2(ij)=g2(ij)-g2(kl)
C                       write(6,*) 'contribution from ',kl,-g2(kl)
                    End If
 270              Continue
                End If
               End If
 230        Continue
 220       Continue
 210     Continue
 200  Continue
      If (iPrint.ge.49) Call TriPrt('Local Hessian ',
     &                               ' ',g2,12)
      nop(1)=nrOpr(kop(1),ioper,nirrep)
      nop(2)=nrOpr(kop(2),ioper,nirrep)
      nop(3)=nrOpr(kop(3),ioper,nirrep)
      nop(4)=nrOpr(kop(4),ioper,nirrep)
*----------------------------------------------------------------------*
*
*-----Distribute contribution to the hessian.
*
*----------------------------------------------------------------------*
      Do 90 iIrrep=0,nIrrep-1
       Do 100 iCn = 1, 4
         Do 110 iCar = 1, 3
          Do 120 jCn =1,iCn
           if (iCn.eq.jCn) Then
               iStop=iCar
           Else
               iStop=3
           End If
           Do 130 jCar=1,istop
            If (IndHss(iCn,iCar,jCn,jCar,iIrrep).ne.0) Then
*----------------------------------------------------------------------*
*
*              Get indexes
*
*----------------------------------------------------------------------*
               ij = Ix(iCn,iCar,jCn,jCar)
               iHess = Abs(IndHss(iCn,iCar,jCn,jCar,iIrrep))
CBS            write(6,*) 'contribution to iHess = ',iHess
*----------------------------------------------------------------------*
*
*            Sign due to integral direction
*
*----------------------------------------------------------------------*
               ps=DBLE(iChTbl(iIrrep,nOp(iCn))*
     &                   iChTbl(iIrrep,nOp(jCn)))
*----------------------------------------------------------------------*
*
*              If over & under triangular integrals are needed
*              multiply by two instead!
*
*----------------------------------------------------------------------*
               if ((iCn.ne.jCn).and.(iCar.eq.jCar).and.
     &               (Abs(indgrd(iCar,iCn,iIrrep)).eq.
     &                Abs(indgrd(jCar,jCn,iIrrep)))) Then
                      ps=ps*Two
               End If
*----------------------------------------------------------------------*
*
*              Sign due to which symmetry group the translation is in
*
*----------------------------------------------------------------------*
               iCh1=2**(iCar-1)
               iCh2=2**(jCar-1)
               ps = ps*xPrmt(kOp(iCn),iCh1)*xPrmt(kOp(jCn),iCh2)
*              ps1 = xPrmt(kOp(jCn),iChBas(1+jCar))
*              ps=ps*ps1*ps2
*----------------------------------------------------------------------*
*
*              Multiply by number of stabilisators.
*
*----------------------------------------------------------------------*
               Fact=ps*DBLE(iStab(iCn))/DBLE(nIrrep*nirrep)
     &                * DBLE(iStab(jCn))
*----------------------------------------------------------------------*
*
*               Add to hessian
*
*----------------------------------------------------------------------*
C              if (abs(Fact * g2(ij)).ge.1d-8)
C    &         write(6,*) 'contribution to iHess = ',iHess, ' from ',
C    &         ij,Fact * g2(ij),g2(ij),Fact
               Hess(iHess) = Hess(iHess) + Fact * g2(ij)
            End If
 130       Continue
 120      Continue
 110     Continue
 100   Continue
  90  Continue
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(iChBas)
      End
