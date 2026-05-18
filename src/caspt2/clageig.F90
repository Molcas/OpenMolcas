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

      SUBROUTINE CLagEig(if_SSDMloc,force_equal,nConf,nRoots,nState,    &
     &                   nLev,CLag,RDMEIG)

      use caspt2_global, only: DREF, DWGT
      use caspt2_global, only: OMGDER, Weight
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: IFSADREF, IFDW, NASHT, ISCF,             &
     &                         JSTATE, ZETA
      use Constants, only: Zero, One, Half

      implicit none

      logical(kind=iwp), intent(in) :: if_SSDMloc, force_equal
      integer(kind=iwp), intent(in) :: nConf, nRoots, nState, nLev
      real(kind=wp), intent(inout) :: CLag(nConf,nRoots)
      real(kind=wp), intent(in) :: RDMEIG(nLev**2)

      real(kind=wp),allocatable :: CI1(:),WRK(:)

      integer(kind=iwp) :: iState
      real(kind=wp) :: WGT, Scal
      real(kind=wp), external :: DDOT_
!
!     MODE=0: Either state-averaged or DWGT matrix
!     MODE=1: XMS-specific term, always state-averaged DM
!
      !! RDMEIG
      call mma_allocate(CI1,nConf,Label='LCI')
      call mma_allocate(WRK,MAX(NLEV,NASHT)**2,Label='WRK')

      Do iState = 1, nState
        If (.not.if_SSDMloc) Then
          if (force_equal .or. .not.IFSADREF) then
            WGT = One/nState ! force equal-weight for XMS
          else if (IFSADREF) then
            WGT = Weight(iState) ! can be unequal weight
          else
            WGT = One/nState ! this should not happen...
          end if
          if (abs(wgt) <= 1.0e-09_wp) cycle
          If (ISCF == 0) Then
            Call LoadCI(CI1,iState)
          Else
            CI1(1) = One
          End If
          WRK(1:NLEV**2) = RDMEIG(1:NLEV**2)*WGT
          Call Poly1_CLag(NCONF,NLEV,CI1,CLag(1,iState),WRK)
        Else
          Wgt = DWgt(iState,jState)
          If (abs(wgt) > 1.0e-09_wp) Then
            If (ISCF == 0) Then
              Call LoadCI(CI1,iState)
            Else
              CI1(1) = One
            End If
            WRK(1:NLEV**2) = RDMEIG(1:NLEV**2)*WGT
            Call Poly1_CLag(NCONF,NLEV,CI1,CLag(1,iState),WRK)
          End If

          !! Derivative of omega for dynamically weighted density
          If (IFDW .and. zeta >= Zero) Then
            If (ISCF == 0) Then
              Call LoadCI(CI1,iState)
            Else
              CI1(1) = One
            End If
            call POLY1(CI1,nConf)
            call GETDREF(DREF,SIZE(DREF))
            Call SQUARE(DREF,WRK,1,nAshT,nAshT)
            !! probably it is doubled somewhere, so should half
            Scal = DDOT_(nAshT**2,RDMEIG,1,WRK,1)*Half
!           write (*,*) 'scal = ', scal
            OMGDER(iState,jState) = OMGDER(iState,jState) + Scal
          End If

        End If
      End Do

      call mma_deallocate(WRK)
!     write(u6,*) 'clag before projection'
!     do istate = 1, nstate
!       write(u6,*) 'state = ', istate
!       do i = 1, nconf
!         write(u6,'(i3,f20.10)') i,clag(i,istate)
!       end do
!     end do
!     write(u6,*) 'debug'
!     IF(ORBIN == 'TRANSFOR') Call CLagX_TrfCI(NCONF,CLAG)
!     if (proj) then
!     ovl = ddot_(nconf*nstate,ci1,1,clag,1)
!     write(u6,*) 'projection coeff = ',ovl
!     call daxpy_(nconf*nstate,-ovl,ci1,1,clag,1)
!     write(u6,*) 'clag after projection'
!     do istate = 1, nstate
!       write(u6,*) 'state = ', istate
!       do i = 1, nconf
!         write(u6,'(i3,f20.10)') i,clag(i,istate)
!       end do
!     end do
!     end if
!
      call mma_deallocate(CI1)

      Return

      End Subroutine CLagEig
