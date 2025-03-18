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
      Subroutine QaaP2MO(G2q,ng2,GDMat,IKL,IKK,ILL)
      use stdalloc, only : mma_allocate, mma_deallocate
      use MCLR_Data, only: nNA
      use input_mclr, only: nRoots
      implicit none
!*****  Input
      INTEGER nG2,IKL,IKK,ILL
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
!*****  Output
      Real*8,DIMENSION(nG2)::G2q
!*****  Auxiliaries
      INTEGER i,j,k,l,ij,kl,ijkl,nD, lMax, itri
      Real*8 Fact
      Real*8,DIMENSION(:),Allocatable::Dsum,Ddif
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)

!****** Calculating dQ_aa/dX_KL, original purpose of this subroutine
      nD=nnA*(nnA+1)/2
      CALL mma_allocate(Dsum,nD)
      CALL mma_allocate(Ddif,nD)
      DO i=1,nnA
       Do j=1,i
        ij=iTri(i,j)
        Dsum(ij)=GDMat(IKL,i,j)+GDMat(IKL,j,i)
        Ddif(ij)=GDMat(IKK,i,j)-GDMat(ILL,i,j)
       End Do
      END DO
       ijkl=0
       do i=1,nna
         do j=1,i
           ij = iTri(i,j)
           do k=1,i
             if(i.eq.k) then
               lmax = j
             else
               lmax = k
             end if
             do l=1,lmax
               kl = iTri(k,l)
               ijkl = ijkl + 1
               fact=0.5d0
               if(k.eq.l) fact=0.25d0
               G2q(ijkl)=                                               &
     & fact*(Dsum(ij)*Ddif(kl)+Dsum(kl)*Ddif(ij))*2.0d0
             end do
           end do
         end do
       end do
      CALL mma_deallocate(Dsum)
      CALL mma_deallocate(Ddif)
      End Subroutine QaaP2MO
