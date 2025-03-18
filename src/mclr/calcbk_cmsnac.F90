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
! Copyright (C) 2021, Paul B Calio                                     *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Based on cmsbk.f from Jie J. Bao &&                            *
! Additional work from  rhs_nac.f                                *
! ****************************************************************
      Subroutine Calcbk_CMSNAC(bk,R,nTri,GDMat,zX)
      use stdalloc, only : mma_allocate, mma_deallocate
      use MCLR_Data, only: nDens2, nNA
      use input_mclr, only: nRoots,ntAsh
      Implicit None

!*****Output
      Real*8,DIMENSION(nDens2)::bk
!*****Input
      Real*8,DIMENSION(nRoots**2)::R
      INTEGER nTri
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
!*****Auxiliaries
      Real*8,DIMENSION(:),Allocatable::FOccMO,P2MOt
      INTEGER nP2,nG1
      Integer i, j, iTri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      ng1=itri(ntash,ntash)
      nP2=itri(ng1,ng1)
      CALL mma_allocate(FOccMO,nDens2)
      CALL mma_allocate(P2MOt,nP2)
      CALL FZero(bk,nDens2)
      CALL GetWFFock_NAC(FOccMO,bk,R,nTri,P2MOt,nP2)
      CALL GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)
      CALL GetPDFTFock_NAC(bk)
      CALL PutCMSFockOcc(FOccMO,nTri)

      CALL mma_deallocate(FOccMO)
      CALL mma_deallocate(P2MOt)

      end subroutine Calcbk_CMSNAC
