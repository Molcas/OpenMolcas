************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICEnSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Oskar Weser                                      *
************************************************************************
#include "intent.h"
      module orthonormalization
        use stdalloc, only : mma_allocate, mma_deallocate
        use fortran_strings, only : to_upper
        use blockdiagonal_matrices, only : t_blockdiagonal, new, delete,
     &    from_raw, to_raw, from_symm_raw, blocksizes

        implicit none
        save
        private
        public ::
     &    t_ON_scheme, ON_scheme, ON_scheme_values,
     &    t_procrust_metric, procrust_metric, metric_values,
     &    orthonormalize, procrust, v_orthonormalize, Grahm_Schmidt

        type :: t_ON_scheme_values
          integer ::
     &      no_ON = 1,
     &      Grahm_Schmidt = 2,
     &      Lowdin = 3
        end type
        type(t_ON_scheme_values), parameter ::
     &    ON_scheme_values = t_ON_scheme_values()

        type :: t_ON_scheme
          integer :: val = ON_scheme_values%Grahm_Schmidt
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

        interface
          real*8 function ddot_(n_,dx,incx_,dy,incy_)
            implicit none
            integer n_, incx_, incy_
            real*8 dx(*), dy(*)
            real*8 ddot
          end function
        end interface

        interface Grahm_Schmidt
          module procedure Grahm_Schmidt_Array, Grahm_Schmidt_Blocks
        end interface


      contains

      function orthonormalize(basis, scheme) result(ONB)
        use general_data, only : nSym, nBAs, nDel, nDelt, nSSH, nOrb
        use rasscf_data, only : nSec, nOrbt, nTot3, nTot4
        implicit none
        type(t_blockdiagonal), intent(in) :: basis(:)
        type(t_ON_scheme), intent(in) :: scheme
        type(t_blockdiagonal) :: ONB(size(basis)), S(size(basis))

        integer :: i, n_to_ON(nSym), n_new(nSym)


        call new(S, blocksizes=blocksizes(basis))
        call read_S(S)

        call new(ONB, blocksizes=blocksizes(basis))

        select case (scheme%val)
          case(ON_scheme_values%Lowdin)
          case(ON_scheme_values%Grahm_Schmidt)
            n_to_ON(:) = nBas(:nSym) - nDel(:nSym)
            call Grahm_Schmidt(basis, S, n_to_ON, ONB, n_new)
            call update_orb_numbers(n_to_ON, n_new, nBas,
     &          nDel, nSSH, nOrb, nDelt, nSec, nOrbt, nTot3, nTot4)
        end select
        call delete(S)

      end function

      function v_orthonormalize(CMO, scheme) result(ONB_v)
        use general_data, only : nBas, nSym
        implicit none
        real*8, intent(in) :: CMO(:)
        type(t_ON_scheme), intent(in) :: scheme
        real*8 :: ONB_v(size(CMO))

        type(t_blockdiagonal) :: basis(nSym)
        type(t_blockdiagonal) :: ONB(nSym)

        call new(basis, blocksizes=nBAS(:nSym))
        call new(ONB, blocksizes=nBAS(:nSym))

        call from_raw(CMO, basis)

        ONB = orthonormalize(basis, scheme)

        call to_raw(ONB, ONB_v)
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
        implicit none
        real*8, intent(in) :: A(:, :), B(:, :)
        type(t_procrust_metric), intent(in) :: metric
        real*8 :: R(size(B, 1), size(B, 2))

        select case (metric%val)
          case (metric_values%Frobenius)
          case (metric_values%Max_4el_trace)
          case default
!            abort_
        end select

        R = matmul(A, B)
      end function procrust

! The elemental keyword automatically loops over the blocks
! and even allows the compiler to parallelize.
      impure elemental subroutine Grahm_Schmidt_Blocks(
     &    basis, S, n_to_ON, ONB, n_new)
        implicit none
        type(t_blockdiagonal), intent(in) :: basis, S
        integer, intent(in) :: n_to_ON
        type(t_blockdiagonal), intent(_OUT_) :: ONB
        integer, intent(out) :: n_new

        call Grahm_Schmidt(
     &        basis%block, S%block,  n_to_ON, ONB%block, n_new)
      end subroutine Grahm_Schmidt_Blocks

      subroutine Grahm_Schmidt_Array(basis, S, n_to_ON, ONB, n_new)
        implicit none
        real*8, intent(in) :: basis(:, :), S(:, :)
        integer, intent(in) :: n_to_ON
        real*8, intent(out) :: ONB(:, :)
        integer, intent(out) :: n_new

        real*8, allocatable :: SCTMP(:), OVL(:)
        real*8 :: L
        integer :: i, nB
        logical :: lin_dep_detected, improve_solution

        nB = size(basis, 1)

        call mma_allocate(SCTMP, nB)
        call mma_allocate(OVL, nB)

        n_new = 0
        ONB(:, n_to_ON + 1 :) = basis(:, n_to_ON + 1 :)
        do i = 1, n_to_ON
          ONB(:, n_new + 1) = basis(:, i)

          improve_solution = .true.
          lin_dep_detected = .false.
          do while (improve_solution .and. .not. lin_dep_detected)
            SCTMP(:nB) = matmul(S, ONB(:, n_new + 1))
! NOTE: One could use DGEMM_ routines,
!       But the matmul routines seem a lot more readable and
!       performance is not really a problem here.
            if (n_new > 0) then
              ovl(:n_new) =
     &          matmul(transpose(ONB(:, :n_new)), sctmp(:nB))
              ONB(:, n_new + 1) =
     &          ONB(:, n_new + 1) - matmul(ONB(:, :n_new) , ovl(:n_new))
            end if
            L = ddot_(nB, SCTMP, 1, ONB(:, n_new + 1), 1)

            lin_dep_detected = L < 1.0d-10
            improve_solution = L < 0.2d0
            if (.not. lin_dep_detected) then
              ONB(:, n_new + 1) = ONB(:, n_new + 1) / sqrt(L)
            end if
            if (.not. (improve_solution .or. lin_dep_detected)) then
              n_new = n_new + 1
            end if
          end do
        end do
        ONB(:, n_new + 1 : n_to_ON) = basis(:, n_new + 1 : n_to_ON)
        call mma_deallocate(SCTMP)
        call mma_deallocate(OVL)
      end subroutine Grahm_Schmidt_Array

      subroutine update_orb_numbers(
     &    n_to_ON, nNew, nBas,
     &    nDel, nSSH, nOrb, nDelt, nSec, nOrbt, nTot3, nTot4)
      use general_data, only : nSym
      implicit none
#include "warnings.fh"
#include "output_ras.fh"
      integer, intent(in) :: n_to_ON(:), nNew(:), nBas(:)
      integer, intent(inout) :: nDel(:), nSSH(:), nOrb(:),
     &  nDelt, nSec, nOrbt, nTot3, nTot4

      integer :: iSym, nSnew(nSym), remove(nSym), total_remove

      remove = n_to_ON(:nSym) - nNew(:nSym)
      total_remove = sum(remove(:nSym))
      if (total_remove > 0) then
        do iSym = 1, nSym
          if (nSSH(iSym) < remove(iSym)) then
            call WarningMessage(2,'Orthonormalization Error')
            Write(LF,*) 'Exact or very near linear dependence '
            Write(LF,*) 'forces RASSCF to stop execution.'
            Write(LF,*) 'Symmetry block:', iSym
            Write(LF,*) 'Effective NR of orthonormal orbs:', nNew(iSym)
            Write(LF,*) 'Earlier number of deleted orbs:', nDel(iSym)
            Write(LF,*) 'Earlier number of secondary orbs:', nSSH(iSym)
            Write(LF,*) 'New number of deleted orbs:',
     &                  nDel(iSym) + remove(iSym)
            Write(LF,*) 'New number of secondary orbs:',
     &                  nSSH(iSym) - remove(iSym)
            call quit(_RC_GENERAL_ERROR_)
          else if (iPRLoc(1) >= usual) then
            call WarningMessage(1,'Orthonormalization Warning')
            Write(LF,*) 'Exact or very near linear dependence'
            Write(LF,*) 'forces RASSCF to delete additional orbitals.'
            Write(LF,*) 'Symmetry block:', iSym
            Write(LF,*) 'Earlier number of deleted orbs =', nDel(iSym)
            Write(LF,*) 'New number of deleted orbs =',
     &                  nDel(iSym) + remove(iSym)
          end if
        end do
        nDel(:nSym) = nDel(:nSym) + remove(:nSym)
        nSSH(:nSym) = nSSH(:nSym) - remove(:nSym)
        nOrb(:nSym) = nOrb(:nSym) - remove(:nSym)
        nDelt = nDelt + total_remove
        nSec = nSec - total_remove
        nOrbt = nOrbt - total_remove
        nTot3 = sum((nOrb(:nSym) + nOrb(:nSym)**2) / 2)
        nTot4 = sum(nOrb(:nSym)**2)
      end if
      end subroutine update_orb_numbers


      subroutine read_raw_S(S_buffer)
        implicit none
        real*8, intent(inout) :: S_buffer(:)
        integer :: i_Rc, i_Opt, i_Component, i_SymLbl
#include "warnings.fh"

        i_Rc = 0
        i_Opt = 2
        i_Component = 1
        i_SymLbl = 1
        Call RdOne(i_Rc, i_Opt, 'Mltpl  0', i_Component,
     &             S_buffer, i_SymLbl)
        if ( i_rc /= 0 ) then
          write(6,*)' RASSCF is trying to orthonormalize orbitals but'
          write(6,*)' could not read overlaps from ONEINT. Something'
          write(6,*)' is wrong with the file, or possibly with the'
          write(6,*)' program. Please check.'
          call quit(_RC_IO_ERROR_READ_)
        end if
      end subroutine


      subroutine read_S(S)
        use general_data, only : nBas, nSym, nActEl
        use rasscf_data, only : nFr, nIn, Tot_Nuc_Charge
        implicit none
#include "warnings.fh"
#include "output_ras.fh"
        type(t_blockdiagonal) :: S(nSym)

        real*8 :: Mol_Charge
        real*8, allocatable :: S_buffer(:)
        integer :: size_S_buffer

        size_S_buffer = sum(nBas(:nSym) * (nBas(:nSym) + 1) / 2)
        call mma_allocate(S_buffer, size_S_buffer + 4)
        call read_raw_S(S_buffer)
        Tot_Nuc_Charge = S_buffer(size_S_buffer + 4)
        call from_symm_raw(S_buffer, S)
        call mma_deallocate(S_buffer)

        Mol_Charge = Tot_Nuc_Charge - dble(2 * (nFr + nIn) + nActEl)
        call put_dscalar('Total Charge    ', Mol_Charge)
        if (IPRLOC(1) >= usual) then
          write(6,*)
          write(6,'(6x,A,f8.2)') 'Total molecular charge',Mol_Charge
        end If
      end subroutine
      end module orthonormalization
