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


      Subroutine TransActMO2(MOs,MOas,mGrid)
#include "nq_info.fh"
******Purpose:
******obtaining an active MO array with a structure of MOs in
******TransActMO from an MO array with a structure of that in
******TransferMO
******Input
      INTEGER mGrid
      Real*8,DIMENSION(mGrid*nOrbt)::MOas
******Output
      Real*8,DIMENSION(mGrid*NASHT)::MOs
******Auxiliary
      INTEGER iIrrep,IOff1,iOff2,iOff3

      DO iGrid=1,mGrid
       IOff3=(iGrid-1)*nAsht
       Do iIrrep=0,mIrrep-1
        IOff2=IOff3+iOff_Ash(iIrrep)+1
        IOff1=(OffOrb(iIrrep)+nIsh(iIrrep))*mGrid+iGrid
        CALL DCopy_(nAsh(iIrrep),MOas(iOff1),mGrid,
     &                           MOs(IOff2) ,1    )
       End Do
      END DO

      RETURN
      End Subroutine


      Subroutine TransferMO(MOas,TabMO,mAO,mGrid,nMOs,iAO)
#include "nq_info.fh"

******Purpose:
******Transferring MO information to MOas to be used in dgemm.
******It records the MO values on each grid point, too.
******But the difference from TransActMO is that the first and
******the second elements are the values of the first MO at grid
******point 1 and grid point 2.

******Input
      INTEGER mAO,mGrid,nMOs,iAO
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
       CALL DCopy_(nCP,TabMO(iAO,1,IOff1),mAO,MOas(IOff2),1)
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
      INTEGER iv,ix,iLoc1,iLoc2,iLoc3

      CALL FZero(D1Unzip,NASHT**2)
      DO iv=1,NASHT
       Do ix=1,iv-1
        iLoc1=(iv-1)*NASHT+ix
        iLoc2=(ix-1)*NASHT+iv
        iLoc3=(iv-1)*iv/2+ix
        D1Unzip(iLoc1)=0.5d0*D1MO(iLoc3)
        D1Unzip(iLoc2)=D1Unzip(iLoc1)
       End Do
       ix=iv
       iLoc1=(iv-1)*NASHT+ix
       iLoc3=(iv+1)*iv/2
       D1Unzip(iLoc1)=0.5d0*D1MO(iLoc3)
      END DO

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

***********************************************************************


***********************************************************************
      Subroutine CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,
     &                        nPMO3p,MOs,MOx,MOy,MOz,TabMO,P2Unzip,
     &                        mAO,mGrid,nMOs,do_grad)
      use nq_pdft, only: lft, lGGA
      Implicit Real*8 (A-H,O-Z)
#include "nq_info.fh"
#include "stdalloc.fh"

******Input
      INTEGER mAO,mGrid,nMOs,nPMO3p
      REAL*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
      Real*8,DIMENSION(NASHT4)::P2Unzip
      Logical do_grad
******Output
      REAL*8,DIMENSION(mGrid*NASHT)::P2MOCube,MOs,MOx,MOy,MOz
      REAL*8,DIMENSION(nPMO3p)::P2MOCubex,P2MOCubey,P2MOCubez

******Auxiliary
      INTEGER iOff1,IOff2,IOff3,IIrrep,nGridPi,NASHT2,NASHT3,icount
      Real*8,DIMENSION(NASHT**3)::P2MO1
      Real*8,DIMENSION(NASHT**2)::P2MOSquare
      Logical lftGGA

      lftGGA=.false.
      IF(lft.and.lGGA) lftGGA=.true.
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


      IF (lGGA) THEN
       DO iGrid=1,mGrid
        IOff1=(iGrid-1)*NASHT
        Do iIrrep=0,mIrrep-1
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

       IF(lftGGA.and.Do_Grad) THEN
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,
     &  P2MOSquare,NASHT,MOx(IOff1),NASHT,
     &  0.0d0,P2MOCubex(iOff1),NASHT)
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,
     &  P2MOSquare,NASHT,MOy(IOff1),NASHT,
     &  0.0d0,P2MOCubey(iOff1),NASHT)
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,
     &  P2MOSquare,NASHT,MOz(IOff1),NASHT,
     &  0.0d0,P2MOCubez(iOff1),NASHT)

        CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,
     &  P2MO1,NASHT,MOx(IOff1),NASHT,
     &  0.0d0,P2MOSquare,NASHT2)
        CALL DGEMM_('T','N',NASHT,1,NASHT,2.0d0,
     &  P2MOSquare,NASHT,MOs(IOff1),NASHT,
     &  1.0d0,P2MOCubex(iOff1),NASHT)

        CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,
     &  P2MO1,NASHT,MOy(IOff1),NASHT,
     &  0.0d0,P2MOSquare,NASHT2)
        CALL DGEMM_('T','N',NASHT,1,NASHT,2.0d0,
     &  P2MOSquare,NASHT,MOs(IOff1),NASHT,
     &  1.0d0,P2MOCubey(iOff1),NASHT)

        CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,
     &  P2MO1,NASHT,MOz(IOff1),NASHT,
     &  0.0d0,P2MOSquare,NASHT2)
        CALL DGEMM_('T','N',NASHT,1,NASHT,2.0d0,
     &  P2MOSquare,NASHT,MOs(IOff1),NASHT,
     &  1.0d0,P2MOCubez(iOff1),NASHT)
       END IF

C       write(6,*) 'P2MOCube array'
C       CALL RecPrt(' ','(10(F9.5,1X))',P2MOCube(IOff1),1,NASHT)
      END DO


      RETURN
      END SUBROUTINE

      Subroutine ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)
      use nq_pdft, only: lft, lGGA

      INTEGER mAO,mGrid,nMOs,iGrid,nAOGrid,iGridOff,iCoordOff
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabSO
      Real*8,DIMENSION(mAO*mGrid*nMOs)::TabSO2

      INTEGER iCoord

      nAOGrid=mAO*mGrid   ! TabSO : mAO*mGrid x nMOs
                          ! TabSO2:

!     Pull out the derivatives and transpose

      Do iGrid=1,mGrid
         IGridOff=(iGrid-1)*mAO*nMOs
         Do iCoord=1,3
            ICoordOff=IGridOff+(iCoord-1)*nMOs+1
            CALL DCopy_(nMOs,TabSO(iCoord+1,iGrid,1),nAOGrid,
     &                       TabSO2(iCoordOff),1)
         End Do
      End Do

      IF(lft.and.lGGA) THEN
       DO iGrid=1,mGrid
        IGridOff=(iGrid-1)*mAO*nMOs
        Do iCoord=4,9
         ICoordOff=IGridOff+(iCoord-1)*nMOs+1
         CALL DCopy_(nMOs,TabSO((iCoord+1),iGrid,1),nAOGrid,
     &                    TabSO2(iCoordOff),1)
        End Do
       END DO
      END IF
      RETURN
      End Subroutine ConvertTabSO
