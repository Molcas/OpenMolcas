************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Oskar Weser                                      *
************************************************************************
      module orthonormalization
        use fortran_strings, only : to_upper
        implicit none
        save
        private
        public ::
     &    t_ON_scheme, ON_scheme,
     &    t_procrust_metric, procrust_metric

        type :: t_ON_scheme
          character(:), allocatable :: val
        end type

        type(t_ON_scheme) :: ON_scheme

        type :: t_procrust_metric
          character(:), allocatable :: val
        end type

        type(t_procrust_metric) :: procrust_metric

      contains

      function orthonormalize(A, S, scheme) result(ONB)
        real*8, intent(in) :: A(:, :)
        real*8, intent(in), optional :: S(:, :)
        type(t_ON_scheme), intent(in), optional :: scheme
        type(t_ON_scheme) :: scheme_
        real*8 :: ONB

        scheme_%val = merge(scheme%val, ON_scheme%val, present(scheme))

        select case (to_upper(trim(scheme_%val)))
          case("LOWDIN")
          case("GRAHM-SCHMITT")
          case default
        end select
        ONB = 1.d0
      end function

!>  Return an orthogonal transformation to make A match B as closely as possible.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  The orthogonal transformation (\f$ T \f$) is given by
!>  the minimization of the distance between (\f$ RA \f$) and
!>  (\f$ B \f$).
!>  The distance is measured by the metric (\f$ d \f$) which
!>  leads to
!>  \f[ T = \text{argmin}\limits_{R \in OG(n)} d(RA, B) \f]
!>  If the metric is induced by the Frobenius-Norm
!>  \f[ T = \text{argmin}\limits_{R \in OG(n)} |RA -  B|_F \f]
!>  it becomes the classical orthogonal procrust's problem.
!>
!>  @paramin[in] A Matrix that should be rotated/mirrored etc.
!>  @paramin[in] B Target matrix.
!>  @paramin[in] metric (Optional parameter) Can be "FROBENIUS", "MAX-4EL-TRACE".
      function procrust(A, B, metric) result(R)
        real*8, intent(in) :: A(:, :), B(:, :)
        type(t_procrust_metric), intent(in), optional :: metric
        type(t_procrust_metric) :: metric_
        real*8 :: R

        metric_%val =
     &    merge(metric%val, procrust_metric%val, present(metric))

        select case (to_upper(trim(metric_%val)))
          case ("FROBENIUS")
          case ("MAX-4EL-TRACE")
          case default
        end select

        R = 1.d0
      end function procrust

      end module orthonormalization
