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
      Subroutine CMSRHSGDMat(GDMat)
      use ipPage, only: W
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: nNA,n2Dens,ipCI,n1Dens
      use MCLR_Data, only: XISPSM
      use input_mclr, only: State_Sym,nRoots,nCSF
      Implicit None
!      Input
!      Output
       Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
!      Auxiliary quantities
       Real*8,DIMENSION(:),Allocatable::GDArray
       Real*8,DIMENSION(n2dens)::rdum
       INTEGER I,J,IOrb,JOrb,NIJ
!      I:state index, "excited state" of state J when I .ne. J
!      IOrb:row index,    orbital index for state I
!      JOrb:column index, orbital index for state J
       INTEGER nConfL,nConfR,iL,iR
       Real*8,Allocatable::CIL(:),CIR(:)
!      (D^IJ)_pq = <I|E_pq|J>, setting I>=J
!       <I|E_pq|J>=<J|E_qp|I>
       Call mma_allocate(GDArray,n1dens)
       iL=state_sym
       iR=state_sym
       nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))
       nConfL=Max(ncsf(iL),nint(xispsm(iL,1)))
       Call mma_allocate(CIR,nConfR)
       Call mma_allocate(CIL,nConfL)
       DO I=1,nRoots
        Call CSF2SD(W(ipCI)%Vec(1+(I-1)*ncsf(iL)),CIL,iL)
        Do J=1,I
        Call CSF2SD(W(ipCI)%Vec(1+(J-1)*ncsf(iR)),CIR,iR)
         Call Densi2_mclr(1,GDArray,rdum,CIL,CIR,0,0,0,n1dens,n2dens)
         NIJ=I*(I-1)/2+J
         do IOrb=1,nnA
          do JOrb=1,nnA
           GDMat(NIJ,IOrb,JOrb)=GDArray((JOrb-1)*nnA+IOrb)
          end do
         end do
        End Do
       END DO
       Call mma_deallocate(GDArray)
       Call mma_deallocate(CIL)
       Call mma_deallocate(CIR)
       END Subroutine CMSRHSGDMat
