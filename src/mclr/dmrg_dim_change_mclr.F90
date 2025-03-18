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
      Subroutine dmrg_dim_change_mclr(orbspc,ndim,iflag)
        Implicit None

        integer::orbspc(8)
        integer::iflag
        integer::ndim
        integer i,n1,n2

!        write(6,*)"==================================================="
!        write(6,*)" Currently, only valid for no symmetry calculation"
!        write(6,*)"==================================================="

!        I remember it was already valid for all symmetry.
!                         -- yma, need check it again 2015.5.14

!        write(6,*)"orbspc",orbspc(1:8)
        ndim=0

        n1 = 0

        if(iflag.eq.0)then
          do i=1,8
            ndim=ndim+orbspc(i)
          end do
        else if(iflag.eq.1)then
          do i=1,1
            n1=n1+orbspc(i)
            ndim=n1*n1
          end do
        else if(iflag.eq.2)then
          do i=1,1
            n1=n1+orbspc(i)
            ndim=n1**4
          end do
        else if(iflag.eq.3)then
          do i=1,1
            n1=n1+orbspc(i)
            ndim=(n1+1)*n1/2
          end do
        else if(iflag.eq.4)then
          do i=1,1
            n1=n1+orbspc(i)
            n2=n1*n1
            ndim=(n2+1)*n2/2
          end do
        else
          write(6,*)"unknow iflag"
          call Quit_OnUserError()
        end if

      End Subroutine dmrg_dim_change_mclr
