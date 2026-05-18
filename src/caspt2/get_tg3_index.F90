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
      !> @brief Calculate linear index for 3-body transition density
      !> matrix elements
      !>
      !> Given 6 active orbital indices (t,u,v,x,y,z), this subroutine
      !> calculates the linear index for the corresponding transition
      !> density matrix element used in MKTG3. The indices are first
      !> converted to three orbital pair indices, which are then sorted
      !> in descending order to determine the final index.
      !>
      !> @param t First active orbital index
      !> @param u Second active orbital index
      !> @param v Third active orbital index
      !> @param x Fourth active orbital index
      !> @param y Fifth active orbital index
      !> @param z Sixth active orbital index
      !> @param nasht Number of active orbitals
      !> @param linear_index Output - the calculated linear index
      subroutine get_tg3_index(t, u, v, x, y, z, nasht, ituvxyz)

        use Definitions, only: iwp
        implicit none

        integer(kind=iwp), intent(in) :: t, u, v, x, y, z, nasht
        integer(kind=iwp), intent(out) :: ituvxyz

        integer(kind=iwp) :: itu, ivx, iyz, jtu, jvx, jyz

        ! Convert individual orbital indices to pair indices
        itu = t + nasht*(u-1)
        ivx = v + nasht*(x-1)
        iyz = y + nasht*(z-1)

        ! Sort the pair indices in descending order (jtu >= jvx >= jyz)
        if (itu < ivx) then
          if (itu >= iyz) then
            jtu = ivx
            jvx = itu
            jyz = iyz
          else if (ivx < iyz) then
            jtu = iyz
            jvx = ivx
            jyz = itu
          else
            jtu = ivx
            jvx = iyz
            jyz = itu
          end if
        else
          if (itu < iyz) then
            jtu = iyz
            jvx = itu
            jyz = ivx
          else if (ivx >= iyz) then
            jtu = itu
            jvx = ivx
            jyz = iyz
          else
            jtu = itu
            jvx = iyz
            jyz = ivx
          end if
        end if

        ! Calculate the linear index using the sorted pair indices
        ituvxyz = ((jtu+1)*jtu*(jtu-1))/6 + (jvx*(jvx-1))/2 + jyz

      end subroutine get_tg3_index
