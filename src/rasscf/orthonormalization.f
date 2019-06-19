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
        use stdalloc, only : mma_allocate, mma_deallocate
        use fortran_strings, only : to_upper
        implicit none
        save
        private
        public ::
     &    t_ON_scheme, ON_scheme, ON_scheme_values,
     &    t_procrust_metric, procrust_metric, metric_values,
     &    orthonormalize, procrust, v_orthonormalize

        type :: t_ON_scheme_values
          integer ::
     &      no_ON = 1,
     &      Grahm_Schmidt = 2,
     &      Lowdin = 3
        end type
        type(t_ON_scheme_values), parameter ::
     &    ON_scheme_values = t_ON_scheme_values()

        type :: t_ON_scheme
          integer :: val = ON_scheme_values%Lowdin
        end type
        type(t_ON_scheme) :: ON_scheme


        type :: t_metric_values
          integer ::
     &      Frobenius = 1,
     &      Max_4el_trace = 2
        end type
        type(t_metric_values), parameter ::
     &    metric_values = t_metric_values()

        type :: t_procrust_metric
          integer :: val = metric_values%Frobenius
        end type
        type(t_procrust_metric) :: procrust_metric


      contains

      function orthonormalize(A, scheme) result(ONB)
        real*8, intent(in) :: A(:, :)
        type(t_ON_scheme), intent(in) :: scheme
        real*8, allocatable :: ONB(:, :), ONB_v(:)

        call mma_allocate(ONB, size(A, 1), size(A, 2))
        call mma_allocate(ONB_v, size(A, 1)**2)

        select case (scheme%val)
          case(ON_scheme_values%Lowdin)
          case(ON_scheme_values%Grahm_Schmidt)
            call ONCMO(pack(A, .true.), ONB_v)
            ONB = reshape(ONB_v, shape(A))
        end select

        call mma_deallocate(ONB_v)
      end function

      function v_orthonormalize(CMO, scheme) result(ONB)
        real*8, intent(in) :: CMO(:)
        type(t_ON_scheme), intent(in) :: scheme
        real*8 :: ONB(size(CMO))

        select case (scheme%val)
          case(ON_scheme_values%Lowdin)
          case(ON_scheme_values%Grahm_Schmidt)
            call ONCMO(CMO, ONB)
        end select
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
        type(t_procrust_metric), intent(in) :: metric
        real*8 :: R

        select case (metric%val)
          case (metric_values%Frobenius)
          case (metric_values%Max_4el_trace)
          case default
!            abort_
        end select

        R = 1.d0
      end function procrust

      end module orthonormalization
