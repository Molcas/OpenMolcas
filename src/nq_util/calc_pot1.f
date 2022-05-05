!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 22, 2021, created this file.               *
! ****************************************************************
       Subroutine Calc_Pot1(Pot1,TabMO,mAO,mGrid,nMOs,P2_ontop,         &
     &                      nP2_ontop,                                  &
     &                      MOs)
      use nq_Grid, only: GradRho,Weights
      use nq_Grid, only: vRho, vSigma
      use nq_pdft
      use nq_Info
#include "stdalloc.fh"
!*****Input
      INTEGER mAO,mGrid,nMOs,nP2_ontop
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
      Real*8,DIMENSION(mGrid*nOrbt)::MOs
      Real*8,DIMENSION(nP2_ontop,mGrid)::P2_ontop
!*****Output
      Real*8,DIMENSION(nPot1)::Pot1
!*****Internal
      Real*8,DIMENSION(:),Allocatable::PreMO
      Real*8 dF_dRhoax,dF_dRhoay,dF_dRhoaz,                             &
     &       dF_dRhobx,dF_dRhoby,dF_dRhobz,                             &
     &       dRdx,dRdy,dRdz,Diff1,dEdRhop2
!     PreMO is MO multiplied with things needed for potential
!     calculation
      INTEGER iGrid,iOff1,iMO,iOrb


      CALL mma_allocate(PreMO,mGrid*nOrbt)
      CALL dcopy_(mGrid*nOrbt,MOs,1,PreMO,1)


!      black terms in the notes
      DO iGrid=1,mGrid
       IF(Pass1(iGrid)) THEN
        dF_dRhoapb(iGrid)=vRho(1,iGrid)+vRho(2,iGrid)
        dF_dRhoamb(iGrid)=vRho(1,iGrid)-vRho(2,iGrid)
        dRdRho(iGrid)=RatioA(iGrid)*(-2.0d0/RhoAB(iGrid))
        dZdRho(iGrid)=dZdR(iGrid)*dRdRho(iGrid)
        dEdRho(iGrid)=dF_dRhoapb(iGrid)+dF_dRhoamb(iGrid)*              &
     &                (ZetaA(iGrid)+RhoAB(iGrid)*dZdRho(iGrid))
         dRdPi(iGrid)=4.0d0/RhoAB(iGrid)**2
       ELSE
        dRdPi(iGrid)=0.0d0
        dF_dRhoapb(iGrid)=0.0d0
        dF_dRhoamb(iGrid)=0.0d0
        dEdRho(iGrid)=0.0d0
       END IF
      END DO

!     red terms in the notes
      IF(lGGA) THEN
       DO iGrid=1,mGrid
        If(Pass1(iGrid)) Then
         dF_dRhoax=2.0D0*vSigma(1,iGrid)*GradRho(1,iGrid)+              &
     &                  vSigma(2,iGrid)*GradRho(4,iGrid)
         dF_dRhobx=2.0D0*vSigma(3,iGrid)*GradRho(4,iGrid)+              &
     &                  vSigma(2,iGrid)*GradRho(1,iGrid)
         dF_dRhoay=2.0D0*vSigma(1,iGrid)*GradRho(2,iGrid)+              &
     &                  vSigma(2,iGrid)*GradRho(5,iGrid)
         dF_dRhoby=2.0D0*vSigma(3,iGrid)*GradRho(5,iGrid)+              &
     &                  vSigma(2,iGrid)*GradRho(2,iGrid)
         dF_dRhoaz=2.0D0*vSigma(1,iGrid)*GradRho(3,iGrid)+              &
     &                  vSigma(2,iGrid)*GradRho(6,iGrid)
         dF_dRhobz=2.0D0*vSigma(3,iGrid)*GradRho(6,iGrid)+              &
     &                  vSigma(2,iGrid)*GradRho(3,iGrid)
         dF_dRhoxapb(iGrid)=dF_dRhoax+dF_dRhobx
         dF_dRhoxamb(iGrid)=dF_dRhoax-dF_dRhobx
         dF_dRhoyapb(iGrid)=dF_dRhoay+dF_dRhoby
         dF_dRhoyamb(iGrid)=dF_dRhoay-dF_dRhoby
         dF_dRhozapb(iGrid)=dF_dRhoaz+dF_dRhobz
         dF_dRhozamb(iGrid)=dF_dRhoaz-dF_dRhobz

         dRhodx(iGrid)=GradRho(1,iGrid)+GradRho(4,iGrid)
         dRhody(iGrid)=GradRho(2,iGrid)+GradRho(5,iGrid)
         dRhodz(iGrid)=GradRho(3,iGrid)+GradRho(6,iGrid)

         GradRhodFdRho(iGrid)=                                          &
     &   (dF_dRhoxamb(iGrid)*dRhodx(iGrid)+                             &
     &    dF_dRhoyamb(iGrid)*dRhody(iGrid)+                             &
     &    dF_dRhozamb(iGrid)*dRhodz(iGrid))
         dEdRho(iGrid)=dEdRho(iGrid)+dZdRho(iGrid)*                     &
     &   GradRhodFdRho(iGrid)

         dEdRhox(iGrid)=dF_dRhoxapb(iGrid)+                             &
     &                  ZetaA(iGrid)*dF_dRhoxamb(iGrid)
         dEdRhoy(iGrid)=dF_dRhoyapb(iGrid)+                             &
     &                  ZetaA(iGrid)*dF_dRhoyamb(iGrid)
         dEdRhoz(iGrid)=dF_dRhozapb(iGrid)+                             &
     &                  ZetaA(iGrid)*dF_dRhozamb(iGrid)

        Else
         dF_dRhoxapb(iGrid)=0.0d0
         dF_dRhoxamb(iGrid)=0.0d0
         dF_dRhoyapb(iGrid)=0.0d0
         dF_dRhoyamb(iGrid)=0.0d0
         dF_dRhozapb(iGrid)=0.0d0
         dF_dRhozamb(iGrid)=0.0d0
         GradRhodFdRho(iGrid)=0.0d0
         dEdRhox(iGrid)=0.0d0
         dEdRhoy(iGrid)=0.0d0
         dEdRhoz(iGrid)=0.0d0
        End If
       END DO
!      green and blue terms in the notes
       If(lft) Then
        DO iGrid=1,mGrid
         if(Pass1(iGrid)) then
          dRdX=dRdRho(iGrid)*dRhodX(iGrid)+                             &
     &         dRdPi(iGrid) *P2_ontop(2,iGrid)
          dRdY=dRdRho(iGrid)*dRhodY(iGrid)+                             &
     &         dRdPi(iGrid) *P2_ontop(3,iGrid)
          dRdZ=dRdRho(iGrid)*dRhodZ(iGrid)+                             &
     &         dRdPi(iGrid) *P2_ontop(4,iGrid)
          GradRdFdRho(iGrid)=dRdX*dF_dRhoxamb(iGrid)+                   &
     &                       dRdY*dF_dRhoyamb(iGrid)+                   &
     &                       dRdZ*dF_dRhozamb(iGrid)
          GradPidFdRho(iGrid)=                                          &
     &    P2_ontop(2,iGrid)*dF_dRhoxamb(iGrid)+                         &
     &    P2_ontop(3,iGrid)*dF_dRhoyamb(iGrid)+                         &
     &    P2_ontop(4,iGrid)*dF_dRhozamb(iGrid)
          d2RdRho2(iGrid)=6.0d0*RatioA(iGrid)/RhoAB(iGrid)**2
          d2RdRhodPi(iGrid)=-2.0d0*dRdPi(iGrid)/RhoAB(iGrid)
          if(Pass2(iGrid)) then
           d2ZdR2(iGrid)=2.0d0*dZdR(iGrid)**3
          else if(Pass3(iGrid)) then
           Diff1=RatioA(iGrid)-ThrsNT
           d2ZdR2(iGrid)=                                               &
     &    (2.0d1*fta*Diff1**2+1.2d1*ftb*Diff1+6.0d0*ftc)*Diff1
          else
           d2ZdR2(iGrid)=0.0d0
          end if
          dEdRho(iGrid)=dEdRho(iGrid)+                                  &
     &    (dZdR(iGrid)+RhoAB(iGrid)*d2ZdR2(iGrid)*dRdRho(iGrid))*       &
     &    GradRdFdRho(iGrid)+                                           &
     &    RhoAB(iGrid)*dZdR(iGrid)*d2RdRho2(iGrid)*                     &
     &    GradRhodFdRho(iGrid)+                                         &
     &    RhoAB(iGrid)*dZdR(iGrid)*d2RdRhodPi(iGrid)*GradPidFdRho(iGrid)
          dEdRhop2=RhoAB(iGrid)*dZdRho(iGrid)
!         end of dEdRho term
!         now dEdRhoprime terms
          dEdRhox(iGrid)=dEdRhox(iGrid)+dEdRhop2*dF_dRhoxamb(iGrid)
          dEdRhoy(iGrid)=dEdRhoy(iGrid)+dEdRhop2*dF_dRhoyamb(iGrid)
          dEdRhoz(iGrid)=dEdRhoz(iGrid)+dEdRhop2*dF_dRhozamb(iGrid)
         else
          GradRdFdRho(iGrid)=0.0d0
          GradPidFdRho(iGrid)=0.0d0
          d2RdRho2(iGrid)=0.0d0
          d2RdRhodPi(iGrid)=0.0d0
          d2ZdR2(iGrid)=0.0d0
         end if
        END DO
       End If
      END IF

      CALL DScal_(mGrid,0.5d0,dEdRho,1)

      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,dEdRho(iGrid),PreMO(iGrid),mGrid)
      END DO

      IF(lGGA) THEN
       DO iIrrep=0,mIrrep-1
        Do iOrb=1,mOrb(iIrrep)
         IOff1=(iOrb+OffOrb(iIrrep)-1)*mGrid
         iMO=iOrb+OffBasFro(iIrrep)
         do iGrid=1,mGrid
          PreMO(IOff1+iGrid)=PreMO(IOff1+iGrid)+                        &
     &     TabMO(2,iGrid,iMO)*dEdRhox(iGrid)+                           &
     &     TabMO(3,iGrid,iMO)*dEdRhoy(iGrid)+                           &
     &     TabMO(4,iGrid,iMO)*dEdRhoz(iGrid)
         end do
        End Do
       END DO
      END IF

      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,Weights(iGrid),PreMO(iGrid),mGrid)
      END DO

      DO iIrrep=0,mIrrep-1
       IOff1=OffOrb(iIrrep)*mGrid+1
       IOff2=OffOrb2(iIrrep)+1
       CALL DGEMM_('T','N',mOrb(iIrrep),mOrb(iIrrep),mGrid,1.0d0,       &
     & PreMO(IOff1),mGrid,MOs(IOff1),mGrid,                             &
     & 1.0d0,Pot1(iOff2),mOrb(iIrrep))
      END DO

      CALL mma_deallocate(PreMO)


      RETURN
      End Subroutine
