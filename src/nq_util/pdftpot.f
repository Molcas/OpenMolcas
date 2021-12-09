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
      Subroutine Calc_Pot2(Pot2,mGrid,Pi,nPi,dF_dRho,ndF_dRho,
     &                     Ratio,Zeta,dRhodX,dRhodY,dRhodZ,MOs,Weights,
     &                     Rhos,dEdRhoax,dEdRhoay,dEdRhoaz,dEdRhobx,
     &                     dEdRhoby,dEdRhobz,Fact1,lGGa)
#include "nq_info.fh"
#include "stdalloc.fh"
******Input
      INTEGER mGrid,nPi,ndF_dRho
      Real*8,DIMENSION(nPi,mGrid)::Pi
      Real*8,DIMENSION(mGrid*nOrbt)::MOs
      Real*8,DIMENSION(ndF_dRho,mGrid)::dF_dRho
      Real*8,DIMENSION(mGrid)::Ratio,Zeta,dRhodX,dRhodY,dRhodZ,Weights,
     &                         Rhos,dEdRhoax,dEdRhoay,dEdRhoaz,dEdRhobx,
     &                         dEdRhoby,dEdRhobz
      Logical lGGA
******Output
      Real*8,DIMENSION(nPot2)::Pot2
******Internal
      Real*8,DIMENSION(mGrid)::Fact1
      Real*8,DIMENSION(:),Allocatable::MOUVX,MOVX,PreMO
      INTEGER iGrid,nGOrb,iOff0,iOff1,iOff2,iOff3,iStack,
     &       nnUVX,iVX,
     &       pIrrep,uIrrep,vIrrep,xIrrep,xMax,puIrrep,
     &       u,v,x,vorb,xorb,ioffu,nporb
      Real*8 ThrsRho,ThrsZ2,ThrsPi,ggaterm

      ThrsRho=1.0d-15
      ThrsZ2=1.0d-15
      ThrsPi=1.0d-30


      DO iGrid=1,mGrid
       IF((Rhos(iGrid).ge.ThrsRho).and.(Pi(1,iGrid).gt.ThrsPi)) THEN
        If((1.0d0-Ratio(iGrid)).gt.ThrsZ2) Then
         if(lGGA) then
          ggaterm=
     &    dRhodX(iGrid)*(dEdRhoBX(iGrid)-dEdRhoAX(iGrid))+
     &    dRhodY(iGrid)*(dEdRhoBY(iGrid)-dEdRhoAY(iGrid))+
     &    dRhodZ(iGrid)*(dEdRhoBZ(iGrid)-dEdRhoAZ(iGrid))
         else
          ggaterm=0.0d0
         end if
          Fact1(iGrid)=Weights(iGrid)/(Zeta(iGrid)*Rhos(iGrid)**2)*
     &   (Rhos(iGrid)*(dF_dRho(2,iGrid)-dF_dRho(1,iGrid))+ggaterm)
        Else
         Fact1(iGrid)=0.0d0
        End If
       ELSE
        Fact1(iGrid)=0.0d0
       END IF
      END DO
      nGOrb=mGrid*nOrbt

C      write(6,*)'Zeta,Rho,Fact1'
C      CALL RecPrt(' ','(10(F4.1,1X))',Zeta ,1,mGrid)
C      CALL RecPrt(' ','(10(F4.1,1X))',Rhos ,1,mGrid)
C      CALL RecPrt(' ','(10(F4.1,1X))',Fact1,1,mGrid)

      CALL mma_allocate(PreMO,nGOrb)
      CALL DCopy_(nGOrb,MOs,1,PreMO,1)

      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,Fact1(iGrid),PreMO(iGrid),mGrid)
      END DO

      CALL mma_allocate(MOVX,nVXt*mGrid)

C      write(6,*)'MOs'
C      CALL RecPrt(' ','(10(F9.4,1X))',MOs,mGrid,nOrbt)
C
C      write(6,*)'PreMO'
C      CALL RecPrt(' ','(10(F9.4,1X))',PreMO,mGrid,nOrbt)
******Calculate Orbv*Orbx at each grid point
      DO vIrrep=0,mIrrep-1
       DO xIrrep=0,vIrrep
        Do v=1,nAsh(vIrrep)
         vorb=v+OffOrb(vIrrep)+nIsh(vIrrep)
         IOff1=(vorb-1)*mGrid
         IF(xIrrep.eq.vIrrep) THEN
          xMax=v
          iStack=v*(v-1)/2
         ELSE
          xMax=nAsh(xIrrep)
          iStack=(v-1)*xMax
         END IF
         Do x=1,xMax
          xorb=x+OffOrb(xIrrep)+nIsh(xIrrep)
          IOff2=(xorb-1)*mGrid
          IOff3=(OffVX(xIrrep,vIrrep)+iStack+x-1)*mGrid
          do iGrid=1,mGrid
           MOVX(iGrid+IOff3)=MOs(iGrid+IOff1)*MOs(iGrid+IOff2)
          end do
         End Do
        End Do
       END DO
      END DO

C      write(6,*)'MOVX'
C      CALL RecPrt(' ','(10(F9.4,1X))',MOVX,mGrid,nVXt)

      CALL mma_allocate(MOUVX,nUVXt*mGrid)
******Calculate Orbu'*Orbv*Orbx at each grid point
******Orbu' is the PreMO array, namely Orbu times the "fact1"
      DO uIrrep=0,mIrrep-1
       IOffU=OffOrb(uIrrep)+nIsh(uIrrep)
       DO vIrrep=0,mIrrep-1
        DO xIrrep=0,vIrrep
         Do iVX=1,nVX(xIrrep,vIrrep)
          IOff1=(Offvx(xIrrep,vIrrep)+iVX-1)*mGrid
          IOff0=OffUVX(xIrrep,vIrrep,uIrrep)+(iVX-1)*nAsh(uIrrep)
          Do u=1,nAsh(uIrrep)
           IOff2=(iOffU+u-1)*mGrid
           IOff3=(IOff0+u-1)*mGrid
      do iGrid=1,mGrid
       MOUVX(iGrid+IOff3)=PreMO(iGrid+IOff2)*MOVX(iGrid+IOff1)
C       write(6,*) MOUVX(iGrid+IOff3),PreMO(iGrid+IOff2),
C     & MOVX(iGrid+IOff1)
      end do
          End Do
         End Do
        End Do
       END DO
      END DO
      CALL mma_deallocate(MOVX )

******Use dgemm to calculate PUVX at this grid point
      DO pIrrep=0,mIrrep-1
       nporb=mOrb(pIrrep)
       IOff1=OffOrb(pIrrep)*mGrid+1
       IOff2=OffPUVX(pIrrep)+1
       DO uIrrep=0,mIrrep-1
        puIrrep=IEOR(pIrrep,uIrrep)
        DO vIrrep=0,mIrrep-1
         xIrrep=IEOR(puIrrep,vIrrep)
         nnUVX=nUVX(xIrrep,vIrrep,uIrrep)
         IOff3=OffUVX(xIrrep,vIrrep,uIrrep)*mGrid+1
       CALL DGEMM_('T','N',npOrb,nnUVX,mGrid,
     &             1.0d0,MOs(iOff1),mGrid,MOUVX(IOff3),mGrid,
     &             1.0d0,Pot2(iOff2),npOrb)
C       write(6,*)'irreps',pIrrep,uIrrep,vIrrep,xIrrep
C       write(6,*)iOff1,iOff3,IOff2
C
C       write(6,*)
C       CALL RecPrt(' ','(10(F9.4,1X))',MOs(iOff1),mGrid,npOrb)
C       write(6,*)
C       CALL RecPrt(' ','(10(F9.4,1X))',MOUVX(iOff3),mGrid,nnUVX)
C       write(6,*)
C       CALL RecPrt(' ','(10(F9.4,1X))',Pot2(IOff2),npOrb,nnUVX)
         IOff2=IOff2+nnUVX*npOrb
        END DO
       END DO
      END DO

C      write(6,*)'PUVX new code'
C      CALL RecPrt(' ','(10(F9.4,1X))',Pot2,1,nPot2)

      CALL mma_deallocate(PreMO)
      CALL mma_deallocate(MOUVX)
      RETURN
      End Subroutine




       Subroutine Calc_Pot1(Pot1,TabMO,mAO,mGrid,nMOs,Pi,nPi,
     &            Rho,nRho,dF_dRho,ndF_dRho,Weights,
     &            OnePz,OneMz,Rhos,Ratio,Zeta,dEdRhoax,dEdRhoay,
     &            dEdRhoaz,dEdRhobx,dEdRhoby,dEdRhobz,
     &            dRhodX,dRhodY,dRhodZ,
     &            MOs,lGGA)
      use nq_Grid, only: GradRho
#include "nq_info.fh"
#include "ksdft.fh"
#include "stdalloc.fh"
******Input
      INTEGER mAO,mGrid,nMOs,nPi,nRho,ndF_dRho
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
      Real*8,DIMENSION(mGrid*nOrbt)::MOs
      Real*8,DIMENSION(mGrid)::Weights
      Real*8,DIMENSION(nPi,mGrid)::Pi
      Real*8,DIMENSION(ndF_dRho,mGrid)::dF_dRho
      Real*8,DIMENSION(nRho,mGrid)::Rho
      Logical lGGA
******Output
      Real*8,DIMENSION(nPot1)::Pot1
      Real*8,DIMENSION(mGrid)::OnePZ,OneMZ,Rhos,Ratio,
     &dEdRhoax,dEdRhoay,dEdRhoaz,dEdRhobx,dEdRhoby,dEdRhobz,Zeta,
     &dRhodX,dRhodY,dRhodZ
******Internal
      Real*8,DIMENSION(:),Allocatable::PreMO
C     PreMO is MO multiplied with things needed for potential
C     calculation
      Real*8,DIMENSION(mGrid)::Fact1,Term2
      INTEGER iGrid,iOff1,iMO,iOrb
      Real*8 ThrsRho,ThrsZ2

      thrsrho=1.0d-15
      thrsZ2=1.0d-15
      CALL mma_allocate(PreMO,mGrid*nOrbt)
      CALL dcopy_(mGrid*nOrbt,MOs,1,PreMO,1)

      DO iGrid=1,mGrid
       Rhos(iGrid)=Rho(1,iGrid)+Rho(2,iGrid)
       IF(Rhos(iGrid).ge.ThrsRho) THEN
        Ratio(iGrid)=4.0d0*Pi(1,iGrid)/(Rhos(iGrid)**2)
        If((1.0d0-Ratio(iGrid)).gt.ThrsZ2) Then
         Zeta(iGrid)=sqrt(1.0d0-Ratio(iGrid))
         OnePz(iGrid)=0.5d0*(1.0d0+Zeta(iGrid))
         OneMz(iGrid)=0.5d0*(1.0d0-Zeta(iGrid))
         Term2(iGrid)=0.5d0*Ratio(iGrid)/Zeta(iGrid)
        Else
         OnePz(iGrid)=0.5d0
         OneMz(iGrid)=0.5d0
         Term2(iGrid)=0.0d0
          Zeta(iGrid)=0.0d0
        End If
        Fact1(iGrid)=dF_dRho(1,iGrid)*(OnePz(iGrid)+Term2(iGrid))+
     &               dF_dRho(2,iGrid)*(OneMz(iGrid)-Term2(iGrid))
       ELSE
        OnePz(iGrid)=0.0d0
        OneMz(iGrid)=0.0d0
        Term2(iGrid)=0.0d0
        Fact1(iGrid)=0.0d0
         Zeta(iGrid)=0.0d0
       END IF
      END DO

      IF(lGGA) THEN
       DO iGrid=1,mGrid
        If(Rhos(iGrid).ge.ThrsRho) Then
         dEdRhoax(iGrid)=2.0D0*dF_dRho(3,iGrid)*GradRho(1,iGrid)+
     &                         dF_dRho(4,iGrid)*GradRho(4,iGrid)
         dEdRhobx(iGrid)=2.0D0*dF_dRho(5,iGrid)*GradRho(4,iGrid)+
     &                         dF_dRho(4,iGrid)*GradRho(1,iGrid)
         dEdRhoay(iGrid)=2.0D0*dF_dRho(3,iGrid)*GradRho(2,iGrid)+
     &                         dF_dRho(4,iGrid)*GradRho(5,iGrid)
         dEdRhoby(iGrid)=2.0D0*dF_dRho(5,iGrid)*GradRho(5,iGrid)+
     &                         dF_dRho(4,iGrid)*GradRho(2,iGrid)
         dEdRhoaz(iGrid)=2.0D0*dF_dRho(3,iGrid)*GradRho(3,iGrid)+
     &                         dF_dRho(4,iGrid)*GradRho(6,iGrid)
         dEdRhobz(iGrid)=2.0D0*dF_dRho(5,iGrid)*GradRho(6,iGrid)+
     &                         dF_dRho(4,iGrid)*GradRho(3,iGrid)

         dRhodx(iGrid)=Rho(3,iGrid)+Rho(6,iGrid)
         dRhody(iGrid)=Rho(4,iGrid)+Rho(7,iGrid)
         dRhodz(iGrid)=Rho(5,iGrid)+Rho(8,iGrid)
         Fact1(iGrid)=Fact1(iGrid)+
     &   ((dEdRhoax(iGrid)-dEdRhobx(iGrid))*dRhodx(iGrid)+
     &    (dEdRhoay(iGrid)-dEdRhoby(iGrid))*dRhody(iGrid)+
     &    (dEdRhoaz(iGrid)-dEdRhobz(iGrid))*dRhodz(iGrid))
     &   *Term2(iGrid)/Rhos(iGrid)
        Else
         dEdRhoax(iGrid)=0.0d0
         dEdRhobx(iGrid)=0.0d0
         dEdRhoay(iGrid)=0.0d0
         dEdRhoby(iGrid)=0.0d0
         dEdRhoaz(iGrid)=0.0d0
         dEdRhobz(iGrid)=0.0d0
        End If
       END DO
      END IF

      CALL DScal_(mGrid,0.5d0,Fact1,1)

C      write(6,*) 'Fact1 for GGA'
C      CALL RecPRt(' ','(10(F9.6,1X)))',Fact1,1,mGrid)

      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,Fact1(iGrid),PreMO(iGrid),mGrid)
      END DO

C      write(6,*)'PreMO after Dscal_'
C      CALL RecPrt(' ','(10(F9.6,1X)))',PreMO,mGrid,nOrbt)

      IF(lGGA) THEN
       DO iIrrep=0,mIrrep-1
        Do iOrb=1,mOrb(iIrrep)
         IOff1=(iOrb+OffOrb(iIrrep)-1)*mGrid
         iMO=iOrb+OffBasFro(iIrrep)
         do iGrid=1,mGrid
          PreMO(IOff1+iGrid)=PreMO(IOff1+iGrid)+
     &     TabMO(2,iGrid,iMO)*
     &(OnePz(iGrid)*dEdRhoax(iGrid)+OneMz(iGrid)*dEdRhobx(iGrid))+
     &     TabMO(3,iGrid,iMO)*
     &(OnePz(iGrid)*dEdRhoay(iGrid)+OneMz(iGrid)*dEdRhoby(iGrid))+
     &     TabMO(4,iGrid,iMO)*
     &(OnePz(iGrid)*dEdRhoaz(iGrid)+OneMz(iGrid)*dEdRhobz(iGrid))
         end do
        End Do
       END DO
C      write(6,*)'PreMO after lGGA'
C      CALL RecPrt(' ','(10(F9.6,1X)))',PreMO,mGrid,nOrbt)
      END IF



      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,Weights(iGrid),PreMO(iGrid),mGrid)
      END DO

C      write(6,*)'PreMO after weights'
C      CALL RecPrt(' ','(10(F9.6,1X)))',PreMO,mGrid,nOrbt)

      DO iIrrep=0,mIrrep-1
       IOff1=OffOrb(iIrrep)*mGrid+1
       IOff2=OffOrb2(iIrrep)+1

C       write(6,*) 'PreMO in irrep',iIrrep+1
C       CALL RecPrt(' ','(10(F9.6,1X)))',
C     & PreMO(IOff1),mGrid,mOrb(iIrrep))
C       write(6,*)
C       CALL RecPrt(' ','(10(F9.6,1X)))',
C     & MOs(IOff1),mGrid,mOrb(iIrrep))

       CALL DGEMM_('T','N',mOrb(iIrrep),mOrb(iIrrep),mGrid,1.0d0,
     & PreMO(IOff1),mGrid,MOs(IOff1),mGrid,
     & 1.0d0,Pot1(iOff2),mOrb(iIrrep))

C       write(6,*)
C       CALL RecPrt(' ','(10(F9.6,1X)))',
C     & Pot1(iOff2),mOrb(iIrrep),mOrb(iIrrep))
      END DO

      CALL mma_deallocate(PreMO)


      RETURN
      End Subroutine


      Subroutine PDFTFock(FI,FA,nFock,Fact1,D1,mGrid,MOas,MOs)
#include "nq_info.fh"

******Input
      INTEGER mGrid,nFock
      Real*8,DIMENSION(mGrid*nOrbt)::MOas
      Real*8,DIMENSION(mGrid*NASHT)::MOs
      Real*8,DIMENSION(mGrid)::Fact1
      Real*8,DIMENSION(NASHT)::D1

******Output
      Real*8,DIMENSION(nFock)::FI,FA
******Intermediate
      Real*8,DIMENSION(mGrid)::Fact2
      Real*8,DIMENSION(mGrid*NASHT)::SumDX
      Real*8 TempD1
      INTEGER iGrid,iIrrep,ik,k,iOff1
      Real*8 ddot_
      External DDot_


******calculate FI. FI=2*pq*sum_k{kk*Fact1}

******TempD1: sum_k{kk}
******Fact2 : 2sum_k{kk}*Fact1
      DO iGrid=1,mGrid
       TempD1=0.0d0
       Do iIrrep=0,mIrrep-1
        do ik=1,nIsh(iIrrep)
         k=ik+OffOrb(iIrrep)
         IOff1=(k-1)*mGrid+iGrid
         TempD1=TempD1+MOas(IOff1)**2
        end do
       End Do
       Fact2(iGrid)=2.0d0*TempD1*Fact1(iGrid)
      END DO
      CALL PDFTFock_Inner(FI,nFock,Fact2,MOas,mGrid)

******calculate FA. FA=pq*sum_vx{vx*Dvx*Fact1}
******First calcualte sum_x{Dvx*x}
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT+1
       CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,
     & D1,NASHT,MOs(IOff1),NASHT,0.0d0,SumDX(iOff1),NASHT)

       Fact2(iGrid)=Fact1(iGrid)*
     & ddot_(NASHT,MOs(iOff1),1,SumDX(iOff1),1)

C       write(6,*)'iGrid',iGrid,Fact1(iGrid),Fact2(iGrid)
C       write(6,*)'D1'
C       CALL RecPrt(' ','(10(F9.5,1X))',D1,NASHT,NASHT)
C       CALL RecPrt(' ','(10(F9.5,1X))',  MOs(iOff1),1,NASHT)
C       CALL RecPrt(' ','(10(F9.5,1X))',SumDX(iOff1),1,NASHT)

      END DO

      CALL PDFTFock_Inner(FA,nFock,Fact2,MOas,mGrid)

      RETURN
      End Subroutine


      Subroutine PDFTFock_Inner(Fock,nFock,Kern,MOas,mGrid)
#include "nq_info.fh"
******Input
      INTEGER nFock,mGrid
      Real*8,DIMENSION(mGrid*nOrbt)::MOas
      Real*8,DIMENSION(mGrid)::Kern
******Output
      Real*8,DIMENSION(nFock)::Fock
******Intermediate
      Real*8,DIMENSION(mGrid*nOrbt)::PreMO
      INTEGER iGrid,iIrrep,iOff1,iOff2

      CALL dcopy_(mGrid*nOrbt,MOas,1,PreMO,1)

      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,Kern(iGrid),PreMO(iGrid),mGrid)
      END DO


      DO iIrrep=0,mIrrep-1
       IOff1=OffOrb(iIrrep)*mGrid+1
       IOff2=OffOrb2(iIrrep)+1

       CALL DGEMM_('T','N',mOrb(iIrrep),mOrb(iIrrep),mGrid,1.0d0,
     & PreMO(IOff1),mGrid,MOas(IOff1),mGrid,
     & 1.0d0,Fock(iOff2),mOrb(iIrrep))

      END DO


      RETURN
      End Subroutine
