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

        type :: t_blockdiagonal_matrix
          real*8, allocatable :: block(:, :)
        end type

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
        real*8 :: R

        select case (metric%val)
          case (metric_values%Frobenius)
          case (metric_values%Max_4el_trace)
          case default
!            abort_
        end select

        R = 1.d0
      end function procrust

      subroutine fill_overlap_matrix(S)
        use general_data, only :
     &    nSym, nBAs, nDel, nActEl, nDelt, nSSH, nDel, nOrb
        implicit none
#include "warnings.fh"
        type(t_blockdiagonal_matrix), intent(inout) :: S(nSym)
        integer :: iSym, i_Rc, i_Opt, i_Component, i_SymLbl,
     &    nB, size_S_buffer, idx_block
        real*8, allocatable :: S_buffer(:)

        size_S_buffer = sum(nBas(:nSym) * (nBas(:nSym) + 1) / 2) + 4
        call mma_allocate(S_buffer, size_S_buffer)

        i_Rc = 0; i_Opt = 2; i_Component = 1; i_SymLbl = 1
        Call RdOne(i_Rc, i_Opt, 'Mltpl  0', i_Component,
     &             S_buffer, i_SymLbl)

        if ( i_rc /= 0 ) then
          Write(6, *)' RASSCF is trying to orthonormalize orbitals but'
          Write(6, *)' could not read overlaps from ONEINT. Something'
          Write(6, *)' is wrong with the file, or possibly with the'
          Write(6, *)' program. Please check.'
          Call Quit(_RC_IO_ERROR_READ_)
        end if

        idx_block = 1
        do iSym = 1, nSym
          nB = nBas(iSym)
          if ( nB > 0 .and. nB - nDel(iSym) > 0) then
            call square(S_buffer(idx_block), S(iSym)%block,1, nB, nB)
          end if
          idx_block = idx_block + (nB**2 + nB) / 2
        end do
        call mma_deallocate(S_buffer)
      end subroutine


      subroutine ONCMO(CMO1,CMO2)
      use general_data, only :
     &    nSym, nBAs, nDel, nActEl, nDelt, nSSH, nDel, nOrb
      use rasscf_data, only :
     &    nSec, nTOT3, nTOT4, Tot_Nuc_Charge, nFr, nIn, nOrbT
      implicit none
      real*8, intent(in) :: CMO1(:)
      real*8, intent(out) :: CMO2(:)
#include "warnings.fh"
#include "output_ras.fh"
      type(t_blockdiagonal_matrix) :: S(nSym)
      integer :: iPRLEV, nBM, nO, NOM, NSBUF, iSYM, nB,
     &    i_Rc, i_Opt, i_Component, i_SymLbl,
     &    xMol_Charge, nDSAVe, NNEGSS, ip_SBUF, ip_CMO, NNEW,
     &    IOLD, IPOLD, IPNEW, NREMOV, ND, NS, NDNEW, NSNEW
      real*8, allocatable :: SBUF(:), SMAT(:), SCTMP(:), OVL(:)

      real*8 :: XNRM2, XSCL
      Parameter (ROUTINE='ONCMO   ')
      Call qEnter('ONCMO')

C Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      CMO2(:) = CMO1(:)

      do iSym = 1, size(S)
        call mma_allocate(S(iSym)%block, nBas(iSym), nBas(iSym))
      end do

      call fill_overlap_matrix(S)

      nBM = maxval(nBas(:nSym))
      nOM = maxval(nBas(:nSym) - nDel(:nSym))
      NSBUF = sum(nBas(:nSym) * (nBas(:nSym) + 1) / 2)

      if (NOM == 0) call qExit(routine)

      call mma_allocate(SBUF, nSBUF + 4)
      call mma_allocate(SMAT, nBM**2)
      call mma_allocate(SCTMP, nBM)
      call mma_allocate(OVL, nBM)
* Read overlap matrix SMAT:
      i_Rc=0
      i_Opt=2
      i_Component=1
      i_SymLbl=1
      Call RdOne(i_Rc,i_Opt,'Mltpl  0',i_Component,SBUF,i_SymLbl)
      If ( i_Rc.ne.0 ) Then
        Write(LF,*)' RASSCF is trying to orthonormalize orbitals but'
        Write(LF,*)' could not read overlaps from ONEINT. Something'
        Write(LF,*)' is wrong with the file, or possibly with the'
        Write(LF,*)' program. Please check.'
        Call Quit(_RC_IO_ERROR_READ_)
      End If
*
*---- Print out nuclear charge
*
      Tot_Nuc_Charge = SBUF(NSBUF + 4)
      xMol_Charge=Tot_Nuc_Charge-DBLE(2*(NFR+NIN)+NACTEL)
      Call put_dscalar('Total Charge    ',xMol_Charge)
      IF(IPRLEV.GE.USUAL) THEN
        Write(LF,*)
        Write(LF,'(6x,A,f8.2)') 'Total molecular charge',xMol_Charge
      End If
*
* Orthonormalize symmetry blocks:
*
      NDSAVE=NDELT
      NNEGSS=0
      ip_SBUF=1
      ip_CMO=1
      Do  iSym=1,nSym
        NB = nBas(iSym)
        NO = NB - nDel(iSym)
        if ( nB > 0 ) then
          IF (nO > 0) THEN
            Call SQUARE(SBUF(ip_SBUF),SMAT,1,NB,NB)
* NNEW=Nr of already orthonormal new CMO''s
            NNEW=0
            DO IOLD=1,NO
              IPOLD=ip_CMO+NB*(IOLD-1)
              IPNEW=ip_CMO+NB*NNEW
              if (ipnew < ipold) then
                 CALL DCOPY_(NB,CMO1(IPOLD),1,CMO2(IPNEW),1)
              end if
  10          CONTINUE
              CALL DGEMM_('N','N',NB,1,NB,1.0D0,S(iSym)%block,NB,
     &                   CMO2(IPNEW),NB,0.0D0,SCTMP,NB)
              IF(NNEW.GT.0) THEN
                CALL DGEMM_('T','N',NNEW,1,NB,1.0D0,CMO2(ip_CMO),NB,
     &                    SCTMP,NB,0.0D0,OVL,NNEW)
                CALL DGEMM_('N','N',NB,1,NNEW,-1.0D0,CMO2(ip_CMO),NB,
     &                    OVL,NNEW,1.0D0,CMO2(IPNEW),NB)
              END IF
              XNRM2=DDOT_(NB,SCTMP,1,CMO2(IPNEW:),1)
              XSCL=1.0D0/SQRT(XNRM2)
              IF (XNRM2.GT.1.0D-10) THEN
                CALL DSCAL_(NB,XSCL,CMO2(IPNEW),1)
                IF(XNRM2.LT.0.2D0) GOTO 10
                NNEW=NNEW+1
              END IF
            END DO
            NREMOV=NO-NNEW
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
     &                        'orbitals.'
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
          END IF
          ip_SBUF=ip_SBUF+(NB*NB+NB)/2
          ip_CMO=ip_CMO+NB*NB
        End If
      End Do
      IF(NNEGSS.GT.0) CALL QUIT(_RC_GENERAL_ERROR_)

      IF (NDSAVE /= NDELT) THEN
        nTot3 = sum((nOrb(:nsym) + nOrb(:nSym)**2) / 2)
        nTot3 = sum(nOrb(:nSym)**2)
      END IF

      call mma_deallocate(SBUF)
      call mma_deallocate(SMAT)
      call mma_deallocate(SCTMP)
      call mma_deallocate(OVL)

      do iSym = 1, size(S)
        call mma_deallocate(S(iSym)%block)
      end do
*
      Call qExit(routine)
      end subroutine ONCMO

      end module orthonormalization
