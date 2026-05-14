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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

module mcpdft_input

use ontop_functional, only: OTFNAL_t
use Definitions, only: wp, iwp, u6

implicit none
private

type :: McpdftInputOptions_t
  logical(kind=iwp) :: wjob = .false.
  logical(kind=iwp) :: mspdft = .false.
  logical(kind=iwp) :: grad = .false.
  logical(kind=iwp) :: meci = .false.
  logical(kind=iwp) :: nac = .false.
  logical(kind=iwp) :: is_hdf5_wfn = .false.
  logical(kind=iwp) :: extparam = .false.
  character(len=256) :: wfn_file = ''
  character(len=256) :: extparamfile = ''
  integer(kind=iwp) :: rlxroot = 0
  integer(kind=iwp) :: nac_states(2) = 0
  type(OTFNAL_t) :: otfnal
end type

type(McpdftInputOptions_t) :: mcpdft_options

public :: mcpdft_options, parse_input

contains

! Reads in the input from the user and sets the appropriate flags in mcpdft_options.
subroutine parse_input()

  use spool, only: Close_LuSpool, Spoolinp
  use text_file, only: next_non_comment
  use unixinfo, only: supername
  use KSDFT_Info, only: CoefR, CoefX
# ifdef _HDF5_
  use mh5, only: mh5_is_hdf5
# endif
  use stdalloc, only: mma_deallocate
  use Constants, only: Zero

  integer(kind=iwp) :: ierror, lu_input
  real(kind=wp) :: lambda
  character(len=80) :: otxc
  character(len=4) :: command
  character(len=:), allocatable :: buffer
  ! lu_input: Logical unit of an ASCII file with a copy of the presently used input.

  ! Resets the input options on subsequent calls
  ! such as in numerical gradients...
  mcpdft_options = mcpdftinputoptions_t(otfnal=OTFNAL_t())

  ierror = 0
  lambda = Zero
  otxc = ''

  call StatusLine('MCPDFT: ','Reading in input')

  call spoolinp(lu_input)
  rewind(lu_input)
  call rdnlst(lu_input,'MCPDFT')

  do

    if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)

    command = buffer(1:min(4,len(buffer)))
    call upcase(command)

    select case (command)
      case ('FILE')
        ! Comsume the following line
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        if (supername(1:18) == 'numerical_gradient') then
          call WarningMessage(1,'Ignoring FILE keyword during numerical gradients')
        else
          ! This will abort if the file does not exist.
          call fileorb(buffer,mcpdft_options%wfn_file)
#         ifdef _HDF5_
          mcpdft_options%is_hdf5_wfn = mh5_is_hdf5(mcpdft_options%wfn_file)
#         endif
        end if
      case ('FUNC','KSDF')
        if (command == 'KSDF') call WarningMessage(1,'Deprecation warning: KSDF should be replaced with FUNC keyword')
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        read(buffer,*) otxc

      case ('LAMB')
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        read(buffer,*,iostat=ierror) lambda
        if (ierror /= 0) call IOError(buffer)

      case ('DFCF')
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        read(buffer,*,iostat=ierror) coefx,coefr
        if (ierror /= 0) call IOError(buffer)

      case ('MSPD')
        mcpdft_options%mspdft = .true.

      case ('WJOB')
        mcpdft_options%wjob = .true.

      case ('GRAD')
        mcpdft_options%grad = .true.

      case ('NAC ')
        mcpdft_options%nac = .true.
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        read(buffer,*,iostat=ierror) mcpdft_options%nac_states(1),mcpdft_options%nac_states(2)
        if (ierror /= 0) call IOError(buffer)
      case ('MECI')
        mcpdft_options%meci = .true.

      case ('EXPM')
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        read(buffer,*,iostat=ierror) mcpdft_options%extparamfile
        call f_inquire(mcpdft_options%extparamfile,mcpdft_options%extparam)
        if (.not. mcpdft_options%extparam) call FileLocatingError(buffer,mcpdft_options%extparamfile)

      case ('RLXR')
        if (.not. next_non_comment(lu_input,buffer)) call EOFError(buffer)
        read(buffer,*,iostat=ierror) mcpdft_options%rlxroot
        if (ierror /= 0) call IOError(buffer)

        ! Done with reading input
      case ('END ')
        exit

      case default
        call warningmessage(2,'Unrecognized keyword: '//command)
        call Quit_OnUserError()
    end select

  end do
  call close_luspool(lu_input)
  call mma_deallocate(buffer)
  mcpdft_options%otfnal = OTFNAL_t(otxc,lambda)

  if (mcpdft_options%rlxroot == 0) call get_iScalar('Relax CASSCF root',mcpdft_options%rlxroot)
  call put_iScalar('NumGradRoot',mcpdft_options%rlxroot)

  mcpdft_options%grad = decide_on_grad(mcpdft_options%grad)

  call verify_input()

end subroutine parse_input

! Validates mcpdft_options object. Ensures that the options provided from the user are valid.
subroutine verify_input()

  use Functionals, only: Get_Func_Type
  use ontop_functional, only: get_base
  use nq_Info, only: meta_GGA_Type1
  use Fock_util_global, only: DoCholesky

  logical(kind=iwp) :: file_exists

  if (mcpdft_options%mspdft) then
    call f_inquire('ROT_HAM',file_exists)
    if (.not. file_exists) then
      call WarningMessage(2,'No H0_Rotate.txt found')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MSPDFT specified but rotated      '
      write(u6,*) ' Hamiltonian file not found!       '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
  end if

  if (mcpdft_options%rlxroot < 1) then
    call WarningMessage(2,'Invalid RLXRoot specified')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' RLXRoot keyword must specify a    '
    write(u6,*) ' valid RASSCF root                 '
    write(u6,*) ' **********************************'
    call Quit_OnUserError()
  end if

  if (mcpdft_options%mspdft .and. mcpdft_options%grad) then
    if (mcpdft_options%otfnal%is_hybrid()) then
      call WarningMessage(2,'Hybrid MS-PDFT gradients not implemented')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MS-PDFT gradients are not         '
      write(u6,*) ' implemented LAMBDA keyword        '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    if (DoCholesky) then
      call WarningMessage(2,'MS-PDFT gradients with density fitting not implemented')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MS-PDFT gradients are not         '
      write(u6,*) ' implemented with density fitting  '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
  end if

  if (mcpdft_options%grad) then
    if (Get_Func_Type(get_base(mcpdft_options%otfnal%otxc)) == meta_GGA_Type1) then
      call WarningMessage(2,'MC-PDFT gradients with translated meta-GGA not supported')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MC-PDFT gradients are not         '
      write(u6,*) ' implemented with meta-GGAs        '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
  end if

  if (mcpdft_options%nac) then
    if (.not. mcpdft_options%mspdft) then
      call WarningMessage(2,'NACs implemented only for MS-PDFT')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' NACs are only implemented         '
      write(u6,*) ' for Multistate PDFT               '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    if (.not. mcpdft_options%grad) then
      call WarningMessage(2,'NACs implemented with GRAD keyword')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' NACs require the GRAD Keywords    '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
  end if
  if (mcpdft_options%meci) then
    if (.not. mcpdft_options%mspdft) then
      call WarningMessage(2,'NACs implemented only for MS-PDFT')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MECI are only implemented         '
      write(u6,*) ' for Multistate PDFT               '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    if (.not. mcpdft_options%grad) then
      call WarningMessage(2,'NACs implemented with GRAD keyword')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MECI require the GRAD Keywords    '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    if (.not. mcpdft_options%nac) then
      call WarningMessage(2,'MECI requires NAC keyword')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MECI require the NAC Keywords     '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
  end if

  if (len_trim(mcpdft_options%wfn_file) /= 0) then
    if (.not. mcpdft_options%is_hdf5_wfn) then
      call WarningMessage(2,'FILE requires hdf5 file')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' FILE keyword only works with hdf5 '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
  end if

end subroutine verify_input

!> @brief decide on calculating potential terms for analytical gradients
!>
!> @details
!>   If grad is true, determines if we still need gradients terms given
!>   request for numerical gradients in GATEWAY, or if this is calculation
!>   is being called from programs LAST_ENERGY or NUMERICAL_GRADIENT
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in] grad whether analytical potential terms were requested to be computed in input
function decide_on_grad(grad)

  use UnixInfo, only: SuperName

  logical(kind=iwp) :: decide_on_grad
  logical(kind=iwp), intent(in) :: grad
  integer(kind=iwp) :: dng
  logical(kind=iwp) :: do_numgrad

  if (.not. grad) then
    decide_on_grad = .false.
    return
  end if
  ! numerical gradients requested in GATEWAY
  call qpg_iscalar('DNG',do_numgrad)
  if (do_numgrad) then
    call get_iscalar('DNG',DNG)
    do_numgrad = (dng == 1)
  end if
  decide_on_grad = .not. (do_numgrad .or. (supername(:11) == 'last_energy') .or. (supername(:18) == 'numerical_gradient'))

end function decide_on_grad

subroutine EOFError(buffer)

  character(len=*), intent(in) :: buffer

  call warningmessage(2,'EOF error when reading line.')
  write(u6,*) 'Last line read from input: ',buffer
  call Quit_OnUserError()

end subroutine EOFError

subroutine IOError(buffer)

  character(len=*), intent(in) :: buffer

  call WarningMessage(2,'I/O error when reading line.')
  write(u6,*) 'Last line read from input: ',buffer
  call Quit_OnUserError()

end subroutine IOError

subroutine FileLocatingError(buffer,filename)

  implicit none
  character(len=*), intent(in) :: buffer, filename

  call WarningMessage(2,'Error in locating file '//filename)
  write(u6,*) 'Last line read from input: ',buffer
  call Quit_OnUserError()

end subroutine FileLocatingError

end module mcpdft_input
