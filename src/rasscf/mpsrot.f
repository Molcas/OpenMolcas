!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!                                                                      *
! Copyright (C) 2016-2017, Stefan Knecht                               *
!                                                                      *
!***********************************************************************
#ifdef _DMRG_
!-------------------------------------------------------------------------------
      subroutine mpsrot(smat,nasht,nash,nsym)


      use qcmaquis_interface_cfg
      use qcmaquis_interface_wrapper

      implicit none

#include "stdalloc.fh"
!-------------------------------------------------------------------------------
!
!  Purpose:  driver routine for the MPS rotation wrt the orbital
!            transformation matrix given in t.
!    --> rotate MPS state ISTATE
!    NOTE: t contains square matrices, one per symmetry
!
!-------------------------------------------------------------------------------

      real*8,  intent(inout) :: smat(nasht*nasht)
      integer, intent(in)    :: nasht
      integer, intent(in)    :: nash(nsym)
      integer, intent(in)    :: nsym

      real*8, allocatable    :: x(:)
      real*8, allocatable    :: t(:)
      real*8, allocatable    :: scr(:)

      integer                :: ioff, iadr, icol, irow, i, j, isym, jorb
      real*8                 :: fac(1,1)

      fac(1,1) = 1.0d0
      !> create file with info for inactive orbital rotation
      call dmrg_interface_ctl(
     &                        task    = 'tra dump',
     &                        x1      = fac,
     &                        ndim    =  1,
     &                        mdim    = -1,
     &                        state   =  0,
     &                        stateL  =  0
     &                       )

      call mma_allocate(x, nasht**2)

      x    = 0
      ioff = 0
      iadr = 1
      icol = 1

#ifdef BLUBB
      !> debug
      do i = 1, nasht
      do j = 1, nasht
        write(6,*) 'row i, column j, value ',i,j,smat(i+(j-1)*nasht)
      end do
      end do
#endif

      do isym = 1, nsym
        if(nash(isym) > 0)then
          irow = icol
          do i = 1, nash(isym)
            iadr = (icol-1)*nasht+irow
            do j = 1, nash(isym)
              x(ioff+nash(isym)*(j-1)+i) = smat(iadr+j-1)
            end do
            icol = icol + 1
          end do
          ioff = ioff + nash(isym)**2
        end if
      end do


      call mma_allocate(scr, nasht**2 +nasht*(nasht+1)/2)
      call mma_allocate(t, nasht**2)
      t = 0; scr = 0;

      jorb = 0
      ioff = 1
      do isym = 1, nsym
        if(isym > 1) then
          ioff = ioff + nash(isym-1)**2
        end if
        if(nash(isym) > 0)then
         call pamtmt(x(ioff),t(ioff),scr,nash(isym))
         !> create file(s) with info for active orbital rotation
#ifdef BLUBB
         !> debug
         write(6,*) 'T ... '
         do i = 1, nash(isym)
         do j = 1, nash(isym)
           write(6,*) 'row, column, value ',i,j,t(ioff-1+i+(j-1)*nasht)
         end do
         end do
#endif
         do i = 1, nash(isym)
            jorb = jorb + 1
            call dmrg_interface_ctl(
     &                              task    = 'tra dump',
     &                              x1      = t(ioff),
     &                              ndim    = nash(isym),
     &                              mdim    = 0,
     &                              state   = jorb,
     &                              stateL  = isym
     &                             )
          end do
        end if
      end do
      call mma_deallocate(x)
      call mma_deallocate(scr)
      call mma_deallocate(t)

      do i = 1, dmrg_state%maxroot
        !> make a backup of the actual MPS
        call dmrg_interface_ctl(
     &                          task    = 'MPS back',
     &                          state   = i,
     &                          stateL  = -1
     &                          )
        !> counterrotate MPS
        write(6,*) ' counterrotate MPS #',i
        call dmrg_interface_ctl(
     &                          task    = 'MPS crot',
     &                          state   = i
     &                          )
      end do

      end subroutine mpsrot
#else
      subroutine mpsrot()
      implicit none
      end subroutine mpsrot
#endif
