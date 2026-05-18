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

      SUBROUTINE CLagEigT(CLag,RDMEIG,SLag,EINACT)

      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NCONF, ISCF, NSTATE
      use Constants, only: One, Two
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: nProcs, Is_Real_Par
#endif

      implicit none

      real(kind=wp), intent(inout) :: CLag(nConf,nState)
      real(kind=wp), intent(in) :: RDMEIG(*), SLag(*), EINACT

      real(kind=wp),allocatable :: CI1(:),CI2(:)
      real(kind=wp) :: Scal
      integer(kind=iwp) :: iStat, jStat

      !! RDMEIG
      call mma_allocate(CI1,nConf,Label='CI1')
      call mma_allocate(CI2,nConf,Label='CI2')
      Do iStat = 1, nState
        If (ISCF == 0) Then
          Call LoadCI(CI1,iStat)
        Else
          CI1(1) = One
        End If
        !! Skip iStat = jStat because SLag is zero
        Do jStat = 1, nState !! iStat-1
          If (ISCF == 0) Then
            Call LoadCI(CI2,jStat)
          Else
            CI2(1) = One
          End If
          !! One of doubling is due to the scaling factor
          !! One of doubling is due to the symmetry of iStat and jStat
!         Scal = SLag(iStat+nState*(jStat-1))*Four
          Scal = SLag(iStat+nState*(jStat-1))*Two
          If (ABS(Scal) <= 1.0e-09_wp) Cycle

          Call Poly1_CLagT(CI1,CI2,                                     &
     &                     CLag(1,iStat),CLag(1,jStat),RDMEIG,Scal)
          !! Inactive terms
#ifdef _MOLCAS_MPP_
          !! The inactive contributions are computed in all processes,
          !! whereas GADGOP will be done later, so divide
          if (is_real_par()) Scal = Scal/real(nProcs,kind=wp)
#endif
          CLag(1:nconf,jStat) = CLag(1:nconf,jStat)                     &
     &      + Scal*EINACT*CI1(1:nconf)
          CLag(1:nconf,iStat) = CLag(1:nconf,iStat)                     &
     &      + Scal*EINACT*CI2(1:nconf)
        End Do
      End Do

      call mma_deallocate(CI1)
      call mma_deallocate(CI2)

      Return

      End Subroutine CLagEigT
