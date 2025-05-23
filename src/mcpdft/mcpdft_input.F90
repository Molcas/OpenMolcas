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
  use definitions,only:iwp,wp
  use ontop_functional,only:OTFNAL_t,func_type,get_base
  use spool,only:Spoolinp,Close_LuSpool

  implicit none
  private

  type :: McpdftInputOptions_t
    logical :: wjob = .false.
    logical :: mspdft = .false.
    logical :: grad = .false.
    logical :: meci = .false.
    logical :: nac = .false.
    logical :: is_hdf5_wfn = .false.
    logical :: extparam = .false.
    character(len=256) :: wfn_file = ''
    character(len=256) :: extparamfile = ''

    integer(kind=iwp) :: rlxroot = 0
    integer(kind=iwp),dimension(2) :: nac_states = 0
    type(OTFNAL_t) :: otfnal

  endtype

  type(McpdftInputOptions_t) :: mcpdft_options

  public :: mcpdft_options,parse_input

contains

  subroutine parse_input()
    ! Reads in the input from the user and sets the appropriate flags in mcpdft_options.
    use constants,only:zero
    use stdalloc,only:mma_deallocate
    use text_file,only:next_non_comment
    use unixinfo,only:supername
    use KSDFT_Info,only:CoefR,CoefX
#ifdef _HDF5_
    use mh5,only:mh5_is_hdf5
#endif

    real(kind=wp) :: lambda
    character(len=80) :: otxc

    ! Logical unit of an ASCII file with a copy of the presently used input.
    integer(kind=iwp) :: lu_input,ierror
    character(len=:),allocatable :: buffer
    character(len=4) :: command

    ! Resets the input options on subsequent calls
    ! such as in numerical gradients...
    mcpdft_options = mcpdftinputoptions_t(otfnal=OTFNAL_t())

    ierror = 0
    lambda = zero
    otxc = ""

    call StatusLine("MCPDFT: ","Reading in input")

    call spoolinp(lu_input)
    rewind(lu_input)
    call rdnlst(lu_input,"MCPDFT")

    do

      if(.not. next_non_comment(lu_input,buffer)) then
        call EOFError(buffer)
      endif

      command = buffer(1:min(4,len(buffer)))
      call upcase(command)

      select case(command)
      case("FILE")
        ! Comsume the following line
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        if(supername(1:18) == 'numerical_gradient') then
          call WarningMessage(1,'Ignoring FILE keyword during numerical gradients')
        else
          ! This will abort if the file does not exist.
          call fileorb(buffer,mcpdft_options%wfn_file)
#ifdef _HDF5_
          mcpdft_options%is_hdf5_wfn = mh5_is_hdf5(mcpdft_options%wfn_file)
#endif
        endif
      case('FUNC','KSDF')
        if(command == 'KSDF') then
          call WarningMessage(1,'Deprecation warning: KSDF should be replaced with FUNC keyword')
        endif
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*) otxc

      case("LAMB")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=ierror) lambda
        if(ierror /= 0) then
          call IOError(buffer)
        endif

      case("DFCF")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=ierror) coefx,coefr
        if(ierror /= 0) then
          call IOError(buffer)
        endif

      case("MSPD")
        mcpdft_options%mspdft = .true.

      case("WJOB")
        mcpdft_options%wjob = .true.

      case("GRAD")
        mcpdft_options%grad = .true.

      case("NAC ")
        mcpdft_options%nac = .true.
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=ierror) mcpdft_options%nac_states(1),mcpdft_options%nac_states(2)
        if(ierror /= 0) then
          call IOError(buffer)
        endif
      case("MECI")
        mcpdft_options%meci = .true.

      case("EXPM")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=iError) mcpdft_options%extparamfile
        call f_inquire(mcpdft_options%extparamfile,mcpdft_options%extparam)
        if(.not. mcpdft_options%extparam) then
          call FileLocatingError(buffer,mcpdft_options%extparamfile)
        endif

      case("RLXR")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=ierror) mcpdft_options%rlxroot
        if(ierror /= 0) then
          call IOError(buffer)
        endif

        ! Done with reading input
      case("END ")
        exit

      case("Default")
        call warningmessage(2,"Unrecognized keyword: "//command)
        call Quit_OnUserError()
      endselect

    enddo
    call close_luspool(lu_input)
    call mma_deallocate(buffer)
    mcpdft_options%otfnal = OTFNAL_t(otxc,lambda)

    if(mcpdft_options%rlxroot == 0) then
      call get_iScalar('Relax CASSCF root',mcpdft_options%rlxroot)
    endif
    call put_iScalar('NumGradRoot',mcpdft_options%rlxroot)

    mcpdft_options%grad = decide_on_grad(mcpdft_options%grad)

    call verify_input()

  endsubroutine

  subroutine verify_input()
    ! Validates mcpdft_options object. Ensures that the options provided from the user are valid.
    use nq_Info,only:meta_GGA_Type1
    use definitions,only:u6
    use Fock_util_global,only:DoCholesky
    implicit none

    logical :: file_exists

    if(mcpdft_options%mspdft) then
      call f_inquire('ROT_HAM',file_exists)
      if(.not. file_exists) then
        call WarningMessage(2,"No H0_Rotate.txt found")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MSPDFT specified but rotated      '
        write(u6,*) ' Hamiltonian file not found!       '
        write(u6,*) ' **********************************'
        call Quit_OnUserError
      endif
    endif

    if(mcpdft_options%rlxroot < 1) then
      call WarningMessage(2,"Invalid RLXRoot specified")
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' RLXRoot keyword must specify a    '
      write(u6,*) ' valid RASSCF root                 '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    endif

    if(mcpdft_options%mspdft .and. mcpdft_options%grad) then
      if(mcpdft_options%otfnal%is_hybrid()) then
        call WarningMessage(2,"Hybrid MS-PDFT gradients not implemented")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MS-PDFT gradients are not         '
        write(u6,*) ' implemented LAMBDA keyword        '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
      if(DoCholesky) then
        call WarningMessage(2,"MS-PDFT gradients with density fitting not implemented")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MS-PDFT gradients are not         '
        write(u6,*) ' implemented with density fitting  '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
    endif

    if(mcpdft_options%grad) then
      if(func_type(get_base(mcpdft_options%otfnal%otxc)) == meta_GGA_Type1) then
        call WarningMessage(2,"MC-PDFT gradients with translated meta-GGA not supported")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MC-PDFT gradients are not         '
        write(u6,*) ' implemented with meta-GGAs        '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
    endif

    if(mcpdft_options%nac) then
      if(.not. mcpdft_options%mspdft) then
        call WarningMessage(2,"NACs implemented only for MS-PDFT")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' NACs are only implemented         '
        write(u6,*) ' for Multistate PDFT               '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
      if(.not. mcpdft_options%grad) then
        call WarningMessage(2,"NACs implemented with GRAD keyword")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' NACs require the GRAD Keywords    '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
    endif
    if(mcpdft_options%meci) then
      if(.not. mcpdft_options%mspdft) then
        call WarningMessage(2,"NACs implemented only for MS-PDFT")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MECI are only implemented         '
        write(u6,*) ' for Multistate PDFT               '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
      if(.not. mcpdft_options%grad) then
        call WarningMessage(2,"NACs implemented with GRAD keyword")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MECI require the GRAD Keywords    '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
      if(.not. mcpdft_options%nac) then
        call WarningMessage(2,"MECI requires NAC keyword")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MECI require the NAC Keywords     '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
    endif

    if(len_trim(mcpdft_options%wfn_file) /= 0) then
      if(.not. mcpdft_options%is_hdf5_wfn) then
        call WarningMessage(2,"FILE requires hdf5 file")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' FILE keyword only works with hdf5 '
        write(u6,*) ' **********************************'
        call Quit_OnUserError()
      endif
    endif

  endsubroutine

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
    use UnixInfo,only:SuperName
    logical(kind=iwp) :: decide_on_grad
    logical(kind=iwp),intent(in) :: grad

    logical(kind=iwp) :: do_numgrad
    integer(kind=iwp) :: dng

    if(.not. grad) then
      decide_on_grad = .false.
      return
    endif
    ! numerical gradients requested in GATEWAY
    call qpg_iscalar('DNG',do_numgrad)
    if(do_numgrad) then
      call get_iscalar('DNG',DNG)
      do_numgrad = (dng == 1)
    endif
    decide_on_grad = .not.(do_numgrad .or. supername(:11) == 'last_energy' .or. supername(:18) == 'numerical_gradient')
  endfunction decide_on_grad

  subroutine EOFError(buffer)
    use definitions,only:u6
    character(len=*),intent(in) :: buffer

    call warningmessage(2,"EOF error when reading line.")
    write(u6,*) "Last line read from input: ",buffer
    call Quit_OnUserError()
  endsubroutine

  subroutine IOError(buffer)
    use definitions,only:u6
    character(len=*),intent(in) :: buffer

    call WarningMessage(2,"I/O error when reading line.")
    write(u6,*) "Last line read from input: ",buffer
    call Quit_OnUserError()
  endsubroutine

  subroutine FileLocatingError(buffer,filename)
    use definitions,only:u6
    implicit none
    character(len=*),intent(in) :: buffer,filename

    call WarningMessage(2,"Error in locating file "//filename)
    write(u6,*) "Last line read from input: ",buffer
    call Quit_OnUserError()
  endsubroutine

endmodule
