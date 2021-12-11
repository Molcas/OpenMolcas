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
* Copyright (C) 2021, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Dec. 08, 2021, created this file.               *
* ****************************************************************
      Subroutine CalcOrbOff()
#include "nq_info.fh"

      INTEGER jOffA_,jOffB_,nTri,iIrrep

      NASHT=0
      jOffA_ = 0
      jOffB_ = 0
      nPot1=0
      nTri=0
      nOrbt=0
      DO iIrrep=0,mIrrep-1
       mOrb(iIrrep)=mBas(iIrrep)-nFro(iIrrep)
       nPot1=nPot1+mOrb(iIrrep)**2
       nOrbt=nOrbt+mOrb(iIrrep)
       NASHT=NASHT+NASH(iIrrep)
       iOff_Ash(iIrrep)=jOffA_
       iOff_Bas(iIrrep)=jOffB_
       OffBasFro(iIrrep)=jOffB_+nFro(iIrrep)
       iOff_BasAct(iIrrep)=jOffB_ + nIsh(iIrrep) + nFro(iIrrep)
       OffOrbTri(iIrrep)=nTri
       nTri=nTri+mOrb(iIrrep)*(mOrb(iIrrep)+1)/2
       jOffA_=jOffA_+nAsh(iIrrep)
       jOffB_=jOffB_+mBas(iIrrep)
      END DO

      OffOrb(0)=0
      OffBas(0)=1
      OffBas2(0)=1
      OffOrb2(0)=0
      DO IIrrep=1,mIrrep-1
       OffBas(iIrrep) =OffBas(iIrrep-1) +mBas(iIrrep-1)
       OffOrb(iIrrep) =OffOrb(iIrrep-1) +mOrb(iIrrep-1)
       OffBas2(iIrrep)=OffBas2(iIrrep-1)+mBas(iIrrep-1)**2
       OffOrb2(iIrrep)=OffOrb2(iIrrep-1)+mOrb(iIrrep-1)**2
      END DO

      RETURN
      End Subroutine

      Subroutine CalcPUVXOff()
#include "nq_info.fh"

      INTEGER IOff1,iIrrep,jIrrep,kIrrep,lIrrep,iOrb,jAct,kAct,lAct,
     &        ijIrrep,klIrrep,nklAct

      IOff1=0
      DO kIrrep=0,mIrrep-1
       kAct=nAsh(kIrrep)
       Do lIrrep=0,kIrrep
        lAct=nAsh(lIrrep)
        nklAct=kAct*lAct
        If(kIrrep.eq.lIrrep) nklAct=kAct*(kAct+1)/2
        OffVX(lIrrep,kIrrep)=IOff1
        nVX(lIrrep,kIrrep)=nklAct
        IOff1=IOff1+nklAct
       End Do
      END DO
      nVXt=iOff1

      IOff1=0
      DO jIrrep=0,mIrrep-1
       jAct=nAsh(jIrrep)
       Do kIrrep=0,mIrrep-1
        kAct=nAsh(kIrrep)
        do lIrrep=0,kIrrep
         lAct=nAsh(lIrrep)
         nklAct=kAct*lAct
         If(kIrrep.eq.lIrrep) nklAct=kAct*(kAct+1)/2
          OffUVX(lIrrep,kIrrep,jIrrep)=IOff1
          nUVX(lIrrep,kIrrep,jIrrep)=jAct*nklAct
          IOff1=iOff1+jAct*nklAct
        end do
       End Do
      END DO
      nUVXt=IOff1

      IOff1=0
      DO iIrrep=0,mIrrep-1
       OffPUVX(iIrrep)=IOff1
       iOrb=mOrb(iIrrep)
       Do jIrrep=0,mIrrep-1
        jAct=nAsh(jIrrep)
        ijIrrep=1+IEOR(iIrrep,jIrrep)
        Do kIrrep=0,mIrrep-1
         kAct=nAsh(kIrrep)
         do lIrrep=0,kIrrep
          lAct=nAsh(lIrrep)
          klIrrep=1+IEOR(kIrrep,lIrrep)
          IF(ijIrrep.eq.klIrrep) THEN
           iOff1=iOff1+iOrb*nUVX(lIrrep,kIrrep,jIrrep)
          END IF
         end do
        End Do
       End Do
      END DO
      nPot2=IOff1

C      write(6,*)'OffPUVX new method',nPot2,MaxUVX
C      write(6,'(8(I5,2X))')(OffPUVX(iIrrep),iIrrep=0,mIrrep-1)
      RETURN
      End Subroutine

      Subroutine TransActMO(MOs,TabMO,mAO,mGrid,nMOs)
#include "nq_info.fh"
******Purpose:
******Trasnferring active orbitals to the MOs array.
******It records the MO values on each grid point.
******The first and the second elements are the MO values
******of the first and the second active MO at grid point 1.
******Input
      INTEGER mAO,mGrid,nMOs
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
******Output
      Real*8,DIMENSION(mGrid*NASHT)::MOs
******Auxiliary
      INTEGER nGridPi,iIrrep,IOff1,iOff2,iOff3


      nGridPi=mAO*mGrid
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT
       Do iIrrep=0,mIrrep-1
        IOff2=IOff_Ash(iIrrep)+1
        IOff3=IOff_BasAct(iIrrep)+1
        CALL DCopy_(nAsh(iIrrep),TabMO(1,iGrid,IOff3),nGridPi,
     &                             MOs(IOff1+IOff2)  ,1)
       End Do
      END DO
      RETURN
      End Subroutine

      Subroutine TransferMO(MOas,TabMO,mAO,mGrid,nMOs)
#include "nq_info.fh"

******Purpose:
******Transferring MO information to MOas to be used in dgemm.
******It records the MO values on each grid point, too.
******But the difference from TransActMO is that the first and
******the second elements are the values of the first MO at grid
******point 1 and grid point 2.

******Input
      INTEGER mAO,mGrid,nMOs
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
******Output
      Real*8,DIMENSION(mGrid*nOrbt)::MOas

******Auxiliary
      INTEGER iIrrep,IOff1,iOff2,iOff3,nCP
      IOff3=0
      DO iIrrep=0,mIrrep-1
       IOff1=OffBasFro(iIrrep)+1
       IOff2=IOff3*mGrid+1
       nCP=mOrb(iIrrep)*mGrid
       CALL DCopy_(nCP,TabMO(1,1,IOff1),mAO,MOas(IOff2),1)
       IOff3=IOff3+mOrb(iIrrep)
      END DO
      RETURN
      End Subroutine


      Subroutine PackPot1(Packed,Full,nPack,Factor)
#include "nq_info.fh"

******Input
      Real*8 Factor
      Real*8,DIMENSION(NPot1)::Full
******Output
      Real*8,DIMENSION(nPack)::Packed
******Auxiliary
      INTEGER iIrrep,p,q,iOff1,IOff2,nOrbs
      DO iIrrep=0,mIrrep-1
       nOrbs=mOrb(iIrrep)
       IOff1=OffOrbTri(iIrrep)
       IOff2=OffOrb2(iIrrep)
       Do P=1,nOrbs
        do Q=1,P
      Packed(IOff1+(P-1)*P/2+Q)=
     &Full(IOff2+(P-1)*nOrbs+Q)+Full(IOff2+(Q-1)*nOrbs+P)
        end do
       End Do
C       write(6,*)'Pot1 in irrep',iirrep
C       CALL RecPrt(' ','(10(F9.5,1X))',Full(iOff2+1),nOrbs,nOrbs)
C       CALL TriPrt(' ','(10(F9.5,1X))',Packed(iOff1+1),nOrbs)
      END DO
      CALL DScal_(nPack,Factor,Packed,1)
      RETURN
      End Subroutine

      Subroutine UnzipD1(D1Unzip,D1MO,nD1MO)
#include "nq_info.fh"

******Input
      INTEGER nD1MO
      Real*8,DIMENSION(nD1MO)::D1MO
******Output
      Real*8,DIMENSION(NASHT**2)::D1Unzip
******Intermediate
      INTEGER iv,ix,iIrrep,IOff1,iOff2,IOffD1,nAct

      CALL FZero(D1Unzip,NASHT**2)
      IOffD1=0
      DO iIrrep=0,mIrrep-1
       nAct=nAsh(iIrrep)
       Do iv=1,nAct
        IOff1=(iv-1)*NASHT
        IOff2=(iv-1)*iv/2
C        do ix=1,iv-1
        do ix=1,iv
         D1Unzip(iOff1+ix)=D1MO(iOff2+ix+IOffD1)*0.5d0
         D1Unzip((ix-1)*nASHT+iv)=D1Unzip(iOff1+ix)
        end do
C        ix=iv
C        D1Unzip(iOff1+ix)=0.5d0*D1MO(iOff2+ix+IOffD1)
       End Do
C       write(6,*)'triangular D1 in irrep',iIrrep
C       CALL TriPrt(' ','(10(F9.5,1X))',D1MO(iOffD1+1),nAct)
       IOffD1=IOffD1+nAct*(nAct+1)/2
      END DO
C      write(6,*)'Rectangular D1'
C      CALL RecPrt(' ','(10(F9.5,1X))',D1Unzip,nAsht,nAsht)

      RETURN
      End Subroutine



      Subroutine UnzipP2(P2Unzip,P2MO,nP2Act)
#include "nq_info.fh"

******Input
      INTEGER nP2Act
      Real*8,DIMENSION(nP2Act)::P2MO
******Output
      Real*8,DIMENSION(NASHT4)::P2Unzip
******AUXILIARY
      INTEGER NASHT2,NASHT3,IOFF1,IOff2,IOff3,
     &I,J,K,L,IAct,JAct,kAct,LAct,iIrrep,jIrrep,kIrrep,lIrrep,
     &IJ,KL,IJKL
      Real*8 Fact

************************************************************************
*                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************

      IF(NASHT4.eq.0) RETURN

      NASHT2=NASHT**2
      NASHT3=NASHT2*NASHT

      DO IIrrep = 0, mIrrep-1
      DO I=1,NASH(iIrrep)
       IAct=iOff_Ash(iIrrep)+I
       IOff1=(IAct-1)*NASHT3
       Do jIrrep = 0, mIrrep-1
       Do J=1,NASH(JIrrep)
        JAct=iOff_Ash(JIrrep)+J
        IOff2=IOff1+(JAct-1)*NASHT2
        IJ=iTri(IAct,JAct)
        do kIrrep = 0, mIrrep-1
        do K=1,NASH(KIrrep)
         KAct=IOff_Ash(KIrrep)+K
         IOff3=IOff2+(KAct-1)*NASHT
         do lIrrep = 0, mIrrep-1
         do L=1,NASH(lIrrep)
          LAct=IOff_Ash(LIrrep)+L
          KL=iTri(KAct,LAct)
          IJKL=iTri(ij,kl)
          Fact=0.5d0
         if((ij.ge.kl).and.(kAct.eq.lAct)) Fact=1.0d0
         if((kl.ge.ij).and.(iAct.eq.jAct)) Fact=1.0d0
          P2Unzip(IOff3+LAct)=P2MO(ijkl)*Fact
         end do
         end do
        end do
        end do
       End Do
       End Do
      END DO
      END DO

      RETURN
      End Subroutine






      Subroutine TranslateDens(Pi,dRhodR,dPi,
     & Thrsrho,ThrsZ2,iTrans,nRho,mGrid,nPi,ndRhodR,nEGrad,
     & DoGrad)
      use nq_Grid, only: Rho, GradRho, nGradRho
******Input
      INTEGER iTrans,nRho,mGrid,nPi,ndRhodR,nEGrad
      REAL*8 ThrsRho,ThrsZ2
      Real*8,DIMENSION(nPi,mGrid)::Pi
      Real*8,DIMENSION(nPi,nEGrad,mGrid)::dPi
      Logical DoGrad
******Input & Output
      Real*8,DIMENSION(ndRhodR,mGrid,nEGrad)::dRhodR
******In-subroutine
      INTEGER iGrid,iEGrad,ngragri,iOff1,nGRho
      Real*8 TempR,RRatio,ScaleFact
      Real*8,DIMENSION(mGrid)::OnePZeta,OneMZeta,Zeta,Ratio,
     &                         RhoAB,dRhodx,dRhody,dRhodz,
     &                         tanhrx,tanhry,tanhrz
      Real*8,DIMENSION(mGrid*nEGrad)::dRhoABdR,dRhoxdR,dRhoydR,dRhozdR,
     &                                dRatio,dZeta

******iTrans
*     1. translated functionals
*     2. fully translated functionals (not implemented yet)
*     3. tanh
*********************************************************************

      nGRho=nGradRho
*********************************************************************
*     calculating total density at each grid
*********************************************************************
      CALL DCopy_(mGrid,Rho(1,1),nRho,RhoAB,1)
      CALL DAXPY_(mGrid,1.0d0,Rho(2,1),nRho,RhoAB,1)


*********************************************************************
*     calculating x, y, z components of density gradient
*********************************************************************
      IF(nGRho.eq.6) THEN
       CALL DCopy_(mGrid,GradRho(1,1),nGRho,dRhodx,1)
       CALL DAXPY_(mGrid,1.0d0,GradRho(4,1),nGRho,dRhodx,1)
       CALL DCopy_(mGrid,GradRho(2,1),nGRho,dRhody,1)
       CALL DAXPY_(mGrid,1.0d0,GradRho(5,1),nGRho,dRhody,1)
       CALL DCopy_(mGrid,GradRho(3,1),nGRho,dRhodz,1)
       CALL DAXPY_(mGrid,1.0d0,GradRho(6,1),nGRho,dRhodz,1)
      END IF


*********************************************************************
*    Ratio and Zeta at each grid point
*********************************************************************
      CALL FZero( Zeta,mGrid)
      CALL FZero(Ratio,mGrid)
      DO iGrid=1,mGrid
       IF(RhoAB(iGrid).ge.ThrsRho) THEN
        RRatio=4.0d0*Pi(1,iGrid)/(RhoAB(iGrid)**2)
        If(iTrans.eq.3) RRatio=tanh(RRatio)
        If((iTrans.eq.1).or.(iTrans.eq.3)) Then
         If((1.0d0-Rratio).gt.ThrsZ2) Zeta(iGrid)=sqrt(1.0d0-Rratio)
        End If
        Ratio(iGrid)=Rratio
       END IF
      END DO

*********************************************************************
*    (1 + zeta)/2 and (1 - zeta)/2
*********************************************************************
      CALL DCopy_(mGrid,[0.5d0],0,OnePZeta,1)
      CALL DCopy_(mGrid,[0.5d0],0,OneMZeta,1)
      CALL DAXPY_(mGrid, 0.5d0,Zeta,1,OnePZeta,1)
      CALL DAXPY_(mGrid,-0.5d0,Zeta,1,OneMZeta,1)


*********************************************************************
*     translating rho_a and rho_b
*********************************************************************
      DO iGrid=1,mGrid
       IF(RhoAB(iGrid).ge.ThrsRho) THEN
        Rho(1,iGrid)=OnePZeta(iGrid)*RhoAB(iGrid)
        Rho(2,iGrid)=OneMZeta(iGrid)*RhoAB(iGrid)
       END IF
      END DO

*********************************************************************
*     translating gradient component of rho_a and rho_b
*********************************************************************
      IF(nGRho.eq.6) THEN
       DO iGrid=1,mGrid
        If(RhoAB(iGrid).ge.ThrsRho) Then
         GradRho(1,iGrid)=OnePZeta(iGrid)*dRhodX(iGrid)
         GradRho(2,iGrid)=OnePZeta(iGrid)*dRhodY(iGrid)
         GradRho(3,iGrid)=OnePZeta(iGrid)*dRhodZ(iGrid)
         GradRho(4,iGrid)=OneMZeta(iGrid)*dRhodX(iGrid)
         GradRho(5,iGrid)=OneMZeta(iGrid)*dRhodY(iGrid)
         GradRho(6,iGrid)=OneMZeta(iGrid)*dRhodZ(iGrid)
        End If
       END DO
      END IF

*********************************************************************
*     Additional terms in the tanh translation
*********************************************************************
      IF(iTrans.eq.3) THEN
       CALL FZero(tanhrx,mGrid)
       CALL FZero(tanhry,mGrid)
       CALL FZero(tanhrz,mGrid)
       DO iGrid=1,mGrid
        If(RhoAB(iGrid).gt.thrsrho) Then
         RRatio=Ratio(iGrid)
         TempR=4.0d0*Pi(1,iGrid)/RhoAB(iGrid)
         TanhrX(iGrid)=(RRatio**2-1.0d0)*(Pi(2,iGrid)-
     &(dRhodX(iGrid)*TempR))/(RhoAB(iGrid)*Zeta(iGrid))
         TanhrY(iGrid)=(RRatio**2-1.0d0)*(Pi(3,iGrid)-
     &(dRhodY(iGrid)*TempR))/(RhoAB(iGrid)*Zeta(iGrid))
         TanhrZ(iGrid)=(RRatio**2-1.0d0)*(Pi(4,iGrid)-
     &(dRhodZ(iGrid)*TempR))/(RhoAB(iGrid)*Zeta(iGrid))
        End If
       END DO
       CALL DAXPY_(mGrid, 1.0d0,TanhrX,1,GradRho(1,1),nRho)
       CALL DAXPY_(mGrid,-1.0d0,TanhrX,1,GradRho(4,1),nRho)
       CALL DAXPY_(mGrid, 1.0d0,TanhrY,1,GradRho(2,1),nRho)
       CALL DAXPY_(mGrid,-1.0d0,TanhrY,1,GradRho(5,1),nRho)
       CALL DAXPY_(mGrid, 1.0d0,TanhrZ,1,GradRho(3,1),nRho)
       CALL DAXPY_(mGrid,-1.0d0,TanhrZ,1,GradRho(6,1),nRho)
      END IF



*********************************************************************
*     calculating terms needed in gradient calculation
*********************************************************************
*     if not doing gradient, code ends here
      IF(.not.DoGrad) RETURN
*********************************************************************
*     calculating density gradient wrt geometrical changes
*********************************************************************
      ngragri=mGrid*nEGrad
      CALL DCopy_(ngragri,dRhodr(1,1,1),ndRhodR,dRhoABdR,1)
      CALL DAXPY_(ngragri,1.0d0,dRhodr(2,1,1),ndRhodR,dRhoABdR,1)

      IF(ndRhodR.eq.8) Then
       CALL DCopy_(ngragri,dRhodr(3,1,1),ndRhodR,dRhoxdR,1)
       CALL DAXPY_(ngragri,1.0d0,dRhodr(6,1,1),ndRhodR,dRhoxdR,1)
       CALL DCopy_(ngragri,dRhodr(4,1,1),ndRhodR,dRhoydR,1)
       CALL DAXPY_(ngragri,1.0d0,dRhodr(7,1,1),ndRhodR,dRhoydR,1)
       CALL DCopy_(ngragri,dRhodr(5,1,1),ndRhodR,dRhozdR,1)
       CALL DAXPY_(ngragri,1.0d0,dRhodr(8,1,1),ndRhodR,dRhozdR,1)
      END IF

*********************************************************************
*    dRatio and dZeta at each grid point
*********************************************************************
*     Calculate dRatio
      CALL Fzero(dRatio,nGraGri)
      DO iGrid=1,mGrid
       IF(RhoAB(iGrid).ge.ThrsRho) THEN
        Do iEGrad=1,nEGrad
         IOff1=(iEGrad-1)*mGrid
         dRatio(IOff1+iGrid)=4.0d0*dPi(1,iEGrad,iGrid)/(RhoAB(iGrid)**2)
     &        -8.0d0*Pi(1,iGrid)*dRhoABdR(IOff1+iGrid)/(RhoAB(iGrid)**3)
        End Do
       END IF
      END DO
*     alculate dZeta
      CALL Fzero(dZeta,nGraGri)
      DO iGrid=1,mGrid
       IF((1.0d0-Ratio(iGrid)).gt.ThrsZ2) THEN
       ScaleFact=-0.5d0/Zeta(iGrid)
       CALL DAxpy_(nEGrad,ScaleFact,dRatio(iGrid),mGrid,
     &                               dZeta(iGrid),mGrid)
       END IF
      END DO

      DO iEGrad=1,nEGrad
       IOff1=(iEGrad-1)*mGrid
       Do iGrid=1,mGrid
        If(RhoAB(iGrid).ge.ThrsRho) Then
         dRhodR(1,iGrid,iEGrad)=OnePZeta(iGrid)*dRhoABdR(IOff1+iGrid)+
     &                          0.50d0*dZeta(IOFf1+iGrid)*RhoAB(iGrid)
         dRhodR(2,iGrid,iEGrad)=OneMZeta(iGrid)*dRhoABdR(IOff1+iGrid)-
     &                          0.50d0*dZeta(IOFf1+iGrid)*RhoAB(iGrid)
        End If
       End Do
      END DO

      IF(ndRhodR.eq.8) THEN
       DO iEGrad=1,nEGrad
        IOff1=(iEGrad-1)*mGrid
        Do iGrid=1,mGrid
         If(RhoAB(iGrid).ge.ThrsRho) Then
          dRhodR(3,iGrid,iEGrad)=OnePZeta(iGrid)*dRhoxdR(IOff1+iGrid)+
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodx(iGrid)
          dRhodR(6,iGrid,iEGrad)=OneMZeta(iGrid)*dRhoxdR(IOff1+iGrid)-
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodx(iGrid)
          dRhodR(4,iGrid,iEGrad)=OnePZeta(iGrid)*dRhoydR(IOff1+iGrid)+
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhody(iGrid)
          dRhodR(7,iGrid,iEGrad)=OneMZeta(iGrid)*dRhoydR(IOff1+iGrid)-
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhody(iGrid)
          dRhodR(5,iGrid,iEGrad)=OnePZeta(iGrid)*dRhozdR(IOff1+iGrid)+
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodz(iGrid)
          dRhodR(8,iGrid,iEGrad)=OneMZeta(iGrid)*dRhozdR(IOff1+iGrid)-
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodz(iGrid)
         End If
        End Do
       END DO
      END IF
      RETURN
      END SUBROUTINE

***********************************************************************


***********************************************************************
      Subroutine CalcP2MOCube(P2MOCube,MOs,MOx,MOy,MOz,TabMO,P2Unzip,
     &                       mAO,mGrid,nMOs)
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "stdalloc.fh"

******Input
      INTEGER mAO,mGrid,nMOs
      REAL*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
      Real*8,DIMENSION(NASHT4)::P2Unzip

******Output
      REAL*8,DIMENSION(mGrid*NASHT)::P2MOCube,MOs,MOx,MOy,MOz

******Auxiliary
      INTEGER iOff1,IOff2,IOff3,IIrrep,nGridPi,NASHT2,NASHT3,icount
      Real*8,DIMENSION(NASHT**3)::P2MO1
      Real*8,DIMENSION(NASHT**2)::P2MOSquare

      nGridPi=mAO*mGrid
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT
       Do iIrrep=0,mIrrep-1
        IOff2=IOff_Ash(iIrrep)+1
        IOff3=IOff_BasAct(iIrrep)+1
        CALL DCopy_(nAsh(iIrrep),TabMO(1,iGrid,IOff3),nGridPi,
     &                             MOs(IOff1+IOff2)  ,1)
        do icount=1,nAsh(iIrrep)
        end do
       End Do
      END DO


      IF (mAO.eq.4) THEN
       DO iGrid=1,mGrid
        IOff1=(iGrid-1)*NASHT
        Do iIrrep=1,mIrrep-1
         IOff2=IOff_Ash(iIrrep)+1
         IOff3=IOff_BasAct(iIrrep)+1
         CALL DCopy_(nAsh(iIrrep),TabMO(2,iGrid,IOff3),nGridPi,
     &                              MOx(IOff1+IOff2)  ,1)
         CALL DCopy_(nAsh(iIrrep),TabMO(3,iGrid,IOff3),nGridPi,
     &                              MOy(IOff1+IOff2)  ,1)
         CALL DCopy_(nAsh(iIrrep),TabMO(4,iGrid,IOff3),nGridPi,
     &                              MOz(IOff1+IOff2)  ,1)
        End Do
       END DO
      END IF

      NASHT2=NASHT**2
      NASHT3=NASHT2*NASHT
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT+1

C       write(6,*) 'MOs array'
C       CALL RecPrt(' ','(10(F9.5,1X))',MOs(IOff1),1,NASHT)
C
C       write(6,*) '2RDM array'
C       CALL RecPrt(' ','(10(F9.5,1X))',P2Unzip,NASHT3,NASHT)

       CALL DGEMM_('T','N',NASHT3,1,NASHT,1.0d0,
     & P2UnZip,NASHT,MOs(IOff1),NASHT,
     & 0.0d0,P2MO1,NASHT3)

C       write(6,*) 'P2MO1 array'
C       CALL RecPrt(' ','(10(F9.5,1X))',P2MO1,NASHT2,NASHT)

       CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,
     & P2MO1,NASHT,MOs(IOff1),NASHT,
     & 0.0d0,P2MOSquare,NASHT2)

C       write(6,*) 'P2MOSquare array'
C       CALL RecPrt(' ','(10(F9.5,1X))',P2MOSquare,NASHT,NASHT)

       CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,
     & P2MOSquare,NASHT,MOs(IOff1),NASHT,
     & 0.0d0,P2MOCube(iOff1),NASHT)

C       write(6,*) 'P2MOCube array'
C       CALL RecPrt(' ','(10(F9.5,1X))',P2MOCube(IOff1),1,NASHT)
      END DO


      RETURN
      END SUBROUTINE

      Subroutine ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)

      INTEGER mAO,mGrid,nMOs,iGrid,nAOGrid,iGridOff,iCoordOff
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabSO
      Real*8,DIMENSION(mAO*mGrid*nMOs)::TabSO2

      INTEGER iCoord

      nAOGrid=mAO*mGrid


      DO iGrid=1,mGrid
       IGridOff=(iGrid-1)*mAO*nMOs
       Do iCoord=1,3
        ICoordOff=IGridOff+(iCoord-1)*nMOs+1
        CALL DCopy_(nMOs,TabSO((iCoord+1),iGrid,1),nAOGrid,
     &                   TabSO2(iCoordOff),1)
       End Do
      END DO
      RETURN
      End Subroutine
