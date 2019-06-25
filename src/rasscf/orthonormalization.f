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
        use blockdiagonal_matrices, only : t_blockdiagonal, new, delete,
     &    fill_from_buffer, fill_from_symm_buffer, fill_to_buffer

        implicit none
        save
        private
        public ::
     &    t_ON_scheme, ON_scheme, ON_scheme_values,
     &    t_procrust_metric, procrust_metric, metric_values,
     &    orthonormalize, procrust, v_orthonormalize, ONCMO

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


      contains

      function orthonormalize(A, scheme) result(ONB)
        implicit none
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
        implicit none
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

      subroutine fill_overlap_matrix(S_buffer, S)
        implicit none
        real*8, intent(in) :: S_buffer(:)
        type(t_blockdiagonal), intent(inout) :: S(:)
        integer :: i_block, block_size, idx_block

        idx_block = 1
        do i_block = 1, size(S)
          block_size = size(S(i_block)%block, 1)
          if (block_size > 0) then
            call square(S_buffer(idx_block), S(i_block)%block,1,
     &                  block_size, block_size)
          end if
          idx_block = idx_block + (block_size**2 + block_size) / 2
        end do
      end subroutine

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


      subroutine ONCMO(oCMO1, oCMO2)
      use general_data, only :
     &    nSym, nBAs, nDel, nActEl, nDelt, nSSH, nDel, nOrb, nTot
      use rasscf_data, only :
     &    nSec, nTOT3, nTOT4, Tot_Nuc_Charge, nFr, nIn, nOrbT
      implicit none
      real*8, intent(in) :: oCMO1(:)
      real*8, intent(out) :: oCMO2(:)
#include "warnings.fh"
#include "output_ras.fh"
      type(t_blockdiagonal) :: S(nSym), CMO1(nSym), CMO2(nSym)
      integer :: iPRLEV, nBM, n_to_ON, NOM, size_S_buffer, iSYM, nB,
     &    nDSAVe, NNEGSS, CMO_block, NNEW, i,
     &    IPOLD, IPNEW, NREMOV, ND, NS, NDNEW, NSNEW
      real*8 :: xMol_Charge
      real*8, allocatable :: S_buffer(:), SCTMP(:), OVL(:)

      logical :: improve_solution, lin_dep_detected

      real*8 :: L, XSCL
      Parameter (ROUTINE='ONCMO   ')

      Call qEnter(ROUTINE)

      if (maxval(nBas(:nSym) - nDel(:nSym)) == 0) call qExit(routine)

      oCMO2(:) = oCMO1(:)
      call new(CMO1, blocksizes=nBas(:nSym))
      call new(CMO2, blocksizes=nBas(:nSym))
      call fill_from_buffer(oCMO1, CMO1)
      call fill_from_buffer(oCMO1, CMO2)

      size_S_buffer = sum(nBas(:nSym) * (nBas(:nSym) + 1) / 2)
      call mma_allocate(S_buffer, size_S_buffer + 4)
      call read_raw_S(S_buffer)

      Tot_Nuc_Charge = S_buffer(size_S_buffer + 4)
      xMol_Charge = Tot_Nuc_Charge - dble(2 * (nFr + nIn) + nActEl)
      Call put_dscalar('Total Charge    ',xMol_Charge)
      if (iprlev >= usual) then
        write(LF,*)
        write(LF,'(6x,A,f8.2)') 'Total molecular charge',xMol_Charge
      end If

      call new(S, blocksizes=nBas(:nSym))
      call fill_from_symm_buffer(S_buffer, S)

      call mma_deallocate(S_buffer)

      call mma_allocate(SCTMP, maxval(nBas(:nSym)))
      call mma_allocate(OVL, maxval(nBas(:nSym)))
* Orthonormalize symmetry blocks:
      do iSym = 1, nSym
        nB = size(CMO2(iSym)%block, 1)
        n_to_ON = nB - nDel(iSym)
        if (nB > 0 .and. n_to_ON > 0) then
* nNew: Nr of already orthonormal new CMOs
          nNew = 0
          do i = 1, n_to_ON
            if (nNew + 1 < i) then
              CMO2(iSym)%block(:, nNew + 1) = CMO1(iSym)%block(:, i)
            end if

            improve_solution = .true.
            lin_dep_detected = .false.
            do while (improve_solution .and. .not. lin_dep_detected)
              SCTMP(:nB) =
     &            matmul(S(iSym)%block, CMO2(iSym)%block(:, nNew + 1))
              if (nnew > 0) then
                ovl(:nNew) = matmul(
     &              transpose(CMO2(iSym)%block(:, :nNew)), sctmp(:nB))
                CMO2(iSym)%block(:, nNew + 1) =
     &              CMO2(iSym)%block(:, nNew + 1)
     &              - matmul(CMO2(iSym)%block(:, :nNew) , ovl(:nNew))
              end if
              L = ddot_(nB, SCTMP, 1,
     &                      CMO2(iSym)%block(:, nNew + 1), 1)

              lin_dep_detected = L < 1.0d-10
              improve_solution = L < 0.2d0
              if (.not. lin_dep_detected) then
                call dscal_(nB, 1.0d0 / sqrt(L),
     &                      CMO2(isym)%block(:, nNew + 1), 1)
              end if
              if (.not. (improve_solution .or. lin_dep_detected)) then
                nNew = nNew + 1
              end if
            end do
          end do

          NDSAVE = NDELT
          NNEGSS = 0
          NREMOV = n_to_ON-NNEW
          IF (NREMOV.GT.0) THEN
            ND=NDEL(ISYM)
            NS=NSSH(ISYM)
            NDNEW=NB-NNEW
            NSNEW=NS+ND-NDNEW
            IF(NSNEW.GE.0) THEN
              IF(IPRLEV.GE.TERSE) THEN
                Call WarningMessage(1,'ONCMO Warning')
                Write(LF,*)' * Exact or very near linear dependence '
                Write(LF,*)' * forces RASSCF to delete additional '//
     &                      'orbitals.'
                Write(LF,*)' *                  Symmetry block:',ISYM
                Write(LF,*)' * Earlier number of deleted orbs =',ND
                Write(LF,*)' *     New number of deleted orbs =',NDNEW
              END IF
            ELSE
              Write(LF,*)' **** ONCMO Error *************************'
              Write(LF,*)' * Exact or very near linear dependence '
              Write(LF,*)' * forces RASSCF to stop execution.'
              Write(LF,*)' *                  Symmetry block:',ISYM
              Write(LF,*)' * Effective nr of orthonormal orbs =',NNEW
              Write(LF,*)' *   Earlier number of deleted orbs =',ND
              Write(LF,*)' * Earlier number of secondary orbs =',NS
              Write(LF,*)' *       New number of deleted orbs =',NDNEW
              Write(LF,*)' *     New number of secondary orbs =',NSNEW
              NNEGSS=NNEGSS+1
            END IF
            NDEL(ISYM)=NDNEW
            NSSH(ISYM)=NSNEW
            NORB(ISYM)=NORB(ISYM)-NREMOV
            NDELT=NDELT+NREMOV
            NSEC =NSEC -NREMOV
            NORBT=NORBT-NREMOV
          END IF
        End If
      End Do
      IF(NNEGSS.GT.0) CALL QUIT(_RC_GENERAL_ERROR_)

      IF (NDSAVE /= NDELT) THEN
        nTot3 = sum((nOrb(:nsym) + nOrb(:nSym)**2) / 2)
        nTot3 = sum(nOrb(:nSym)**2)
      END IF

      call mma_deallocate(SCTMP)
      call mma_deallocate(OVL)

      call fill_to_buffer(CMO2, oCMO2)
      call delete(CMO1)
      call delete(CMO2)
      call delete(S)

      Call qExit(routine)
      end subroutine ONCMO

      subroutine Grahm_Schmidt(basis, S, n_to_ON, ONB, n_new)
        implicit none
        real*8, intent(in) :: basis(:, :), S(:, :)
        integer, intent(in) :: n_to_ON
        real*8, intent(out) :: ONB(:, :)
        integer, intent(out) :: n_new

        real*8, allocatable :: SCTMP(:), OVL(:)

        integer :: i, nB
        real*8 :: L

        logical :: lin_dep_detected, improve_solution

        nB = size(basis, 1)

        call mma_allocate(SCTMP, nB)
        call mma_allocate(OVL, nB)

        n_new = 0
        do i = 1, n_to_ON
          if (n_new + 1 < i) then
            ONB(:, n_new + 1) = basis(:, i)
          end if

          improve_solution = .true.
          lin_dep_detected = .false.
          do while (improve_solution .and. .not. lin_dep_detected)
            SCTMP(:nB) = matmul(S, ONB(:, n_new + 1))
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
              call dscal_(nB, 1.0d0 / sqrt(L), ONB(:, n_new + 1), 1)
            end if
            if (.not. (improve_solution .or. lin_dep_detected)) then
              n_new = n_new + 1
            end if
          end do
        end do

        call mma_deallocate(SCTMP)
        call mma_deallocate(OVL)


      end subroutine

      end module orthonormalization
