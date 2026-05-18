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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      Subroutine Recompute_DPT2AO(NDIM,DPT2,DPT2C,DPT2AO,DPT2CAO)

      use definitions, only: wp, iwp
      use stdalloc, only: mma_allocate,mma_deallocate
      use caspt2_global, only: LUONEM
      use caspt2_module, only: NSYM, NBAS, IAD1M
      use Constants, only: Zero

      implicit none

      integer(kind=iwp), intent(in) :: NDIM ! = NBSQT
      real(kind=wp), intent(inout) :: DPT2(NDIM), DPT2C(NDIM),          &
     &  DPT2AO(NDIM), DPT2CAO(NDIM)

      real(kind=wp), allocatable :: WRK1(:), WRK2(:), WRK3(:), WRK4(:)

      real(kind=wp) :: val
      integer(kind=iwp) :: IDISK, iBasTr, iBasSq, iSym, nBasI, liBasTr, &
     &                     liBasSq, ljBasSq, iBasI, jBasI, NBSQT

      NBSQT = NDIM
      call mma_allocate(WRK1,NDIM,Label='WRK1')
      call mma_allocate(WRK2,NDIM,Label='WRK2')
      call mma_allocate(WRK3,NDIM,Label='WRK3')
      call mma_allocate(WRK4,NDIM,Label='WRK4')
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,WRK1,NDIM,IDISK)

      DPT2AO(1:NDIM) = Zero
      DPT2CAO(1:NDIM) = Zero

      iBasTr = 1
      iBasSq = 1
      Do iSym = 1, nSym
        call OLagTrf(1,iSym,NBSQT,WRK1,DPT2 ,WRK3,WRK2)
        call OLagTrf(1,iSym,NBSQT,WRK1,DPT2C,WRK4,WRK2)
        nBasI = nBas(iSym)
        liBasTr = iBasTr
        liBasSq = iBasSq
        ljBasSq = iBasSq
        Do iBasI = 1, nBasI
          Do jBasI = 1, iBasI
            liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
            ljBasSq = iBasSq + jBasI-1 + nBasI*(iBasI-1)
            If (iBasI == jBasI) Then
              DPT2AO (liBasTr) = WRK3(liBasSq)
              DPT2CAO(liBasTr) = WRK4(liBasSq)
            Else
              val = WRK3(liBasSq)+WRK3(ljBasSq)
              DPT2AO (liBasTr) = val
              val = WRK4(liBasSq)+WRK4(ljBasSq)
              DPT2CAO(liBasTr) = val
            End If
            liBasTr = liBasTr + 1
          End Do
        End Do
        iBasTr = iBasTr + nBasI*(nBasI+1)/2
      End Do

      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
      call mma_deallocate(WRK3)
      call mma_deallocate(WRK4)

      End Subroutine Recompute_DPT2AO
