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
! Copyright (C) 2022, Andy Kaiser                                      *
!***********************************************************************
      module frenkel_global_vars

        use Definitions, only: iwp, wp

        implicit none
        private

        integer(kind=iwp) :: ityp, jtyp, valst, corest, nesta, nestb, &
                             nh1
        integer(kind=iwp), allocatable :: nestla(:), nestlb(:)
        logical :: docoul, doexcitonics, doexch, labb, &
        rixs, laba, excl, aux2
        real(kind=wp), allocatable :: vnucb(:), enucb(:)

        public :: ityp, jtyp, valst, corest, nesta, nestb, nh1, &
                  nestla, nestlb, docoul, doexcitonics, doexch, labb, &
                  rixs, laba, excl, aux2, vnucb, enucb

      end module frenkel_global_vars
