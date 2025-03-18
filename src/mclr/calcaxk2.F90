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
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      subroutine CalcAXk2(AXk,D1,D2,PUVX,NPUVX,IndPUVX,Off_Act,Off_Orb)
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: nDens2, nNA, ipMat
      use input_mclr, only: ntBas,ntAsh,nSym,nAsh,nOrb,nIsh
      Implicit None
      Integer NPUVX
      Real*8,DIMENSION(NPUVX)::PUVX
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh)::IndPUVX
      Real*8,DIMENSION(ntAsh**2)::D1,D2
      INTEGER,DIMENSION(nSym)::Off_Act,Off_Orb

      Real*8,DIMENSION(nDens2)::AXk
      Real*8,DIMENSION(:),Allocatable::Opu
      INTEGER p,q,iSym,t,u,v,x,ip,iq,it,loc1,loc2
      Real*8 tempa

      CALL mma_allocate(Opu,ntBas*ntAsh)
      DO p=1,ntBas
       Do u=1,ntAsh
        tempa=0.0d0
        do v=1,ntAsh
         do x=1,ntAsh
          if(IndPUVX(p,u,v,x).ne.0)                                     &
     &tempa=tempa+PUVX(IndPUVX(p,u,v,x))*D2((v-1)*nnA+x)
         end do
        end do
        Opu((p-1)*ntAsh+u)=tempa
       End Do
      END DO
      DO iSym=1,nSym
       Do iq=1,nAsh(iSym)
       q=iq+Off_Act(ISym)
        Do ip=1,nIsh(iSym)  ! p is inactive
         p=ip+Off_Orb(ISym)
         tempa=0.0d0
         do it=1,nAsh(ISym)
          t=it+Off_Act(ISym)
          tempa=tempa+(D1((t-1)*nnA+q)+D1((q-1)*nnA+t))                 &
     &         *Opu((p-1)*nnA+t)
         end do
!         write(6,*)'tempa after sum over t',tempa
         loc1=ipMat(iSym,iSym)+(iq-1)*nOrb(iSym)+ip-1                   &
     &   +nOrb(iSym)*NIsh(iSym)
         loc2=ipMat(iSym,iSym)+(ip-1)*nOrb(iSym)+iq-1                   &
     &   +NIsh(iSym)
         AXK(loc1)=AXK(loc1)+tempa
         AXK(loc2)=AXK(loc2)-tempa
        End Do
        Do ip=nIsh(iSym)+nAsh(iSym)+1,nOrb(iSym)  ! p is virtual
         p=ip+Off_Orb(ISym)
         tempa=0.0d0
         do it=1,nAsh(ISym)
          t=it+Off_Act(ISym)
          tempa=tempa+(D1((t-1)*nnA+q)+D1((q-1)*nnA+t))                 &
     &         *Opu((p-1)*nnA+t)
         end do
         loc1=ipMat(iSym,iSym)+(iq-1)*nOrb(iSym)+ip-1                   &
     &   +nOrb(iSym)*NIsh(iSym)
         loc2=ipMat(iSym,iSym)+(ip-1)*nOrb(iSym)+iq-1                   &
     &   +NIsh(iSym)
         AXK(loc1)=AXK(loc1)+tempa
         AXK(loc2)=AXK(loc2)-tempa
        End Do
       END Do
      END DO
      CALL mma_deallocate(Opu)
      END SUBROUTINE CalcAXk2
