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
  use ontop_functional,only:OTFNAL_t

  implicit none
  private

  type :: McpdftInputOptions_t
    logical :: wjob = .false.
    logical :: mspdft = .false.
    logical :: grad = .false.
    logical :: meci = .false.
    logical :: nac = .false.
    logical :: is_hdf5_wfn = .false.
    character(len=256) :: wfn_file = ""

    integer(kind=iwp),dimension(2) :: nac_states = 0
    type(OTFNAL_t) :: otfnal

  endtype

  type(McpdftInputOptions_t) :: mcpdft_options

  public :: mcpdft_options,parse_input

contains

  subroutine parse_input()
    ! Reads in the input from the user and sets the appropriate flags in mcpdft_options.

    use stdalloc,only:mma_deallocate
    use text_file,only:next_non_comment
    use KSDFT_Info,only:CoefR,CoefX
#ifdef _HDF5_
    use mh5,only:mh5_is_hdf5
#endif
    implicit none

    real(kind=wp) :: lambda = 0.0d0
    character(len=80) :: otxc = ""

    ! Logical unit of an ASCII file with a copy of the presently used input.
    integer(kind=iwp) :: lu_input,ierror = 0
    character(len=:),allocatable :: buffer
    character(len=4) :: command

    call StatusLine("MCPDFT:","Reading in input")

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
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        ! This will abort if the file does not exist.
        call fileorb(buffer,mcpdft_options%wfn_file)
#ifdef _HDF5_
        mcpdft_options%is_hdf5_wfn = mh5_is_hdf5(mcpdft_options%wfn_file)
#endif

      case("KSDF")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*) otxc

      case("LAMB")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=iError) lambda
        if(iError /= 0) then
          call IOError(buffer)
        endif

      case("DFCF")
        if(.not. next_non_comment(lu_input,buffer)) then
          call EOFError(buffer)
        endif
        read(buffer,*,IOStat=iError) coefx,coefr
        if(iError /= 0) then
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
        read(buffer,*,IOStat=iError) mcpdft_options%nac_states(1),mcpdft_options%nac_states(2)
        if(iError /= 0) then
          call IOError(buffer)
        endif

      case("MECI")
        mcpdft_options%meci = .true.

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

    call verify_input()

  endsubroutine

  subroutine verify_input()
    ! Validates mcpdft_options object. Ensures that the options provided from the user are valid.
    use definitions,only:u6
    use Fock_util_global,only:DoCholesky
    implicit none

    if(mcpdft_options%mspdft .and. mcpdft_options%grad) then
      if(mcpdft_options%otfnal%is_hybrid()) then
        call WarningMessage(2,"Hybrid MS-PDFT gradients not implemented")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MS-PDFT gradients are not         '
        write(u6,*) ' implemented LAMBDA keyword        '
        write(u6,*) ' **********************************'
      endif
      if(DoCholesky) then
        call WarningMessage(2,"MS-PDFT gradients with density fitting not implemented")
        write(u6,*) ' ************* ERROR **************'
        write(u6,*) ' MS-PDFT gradients are not         '
        write(u6,*) ' implemented density fitting       '
        write(u6,*) ' **********************************'
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

  subroutine EOFError(buffer)
    use definitions,only:u6
    implicit none
    character(len=*),intent(in) :: buffer

    call warningmessage(2,"EOF error when reading line.")
    write(u6,*) "Last line read from input: ",buffer
    call Quit_OnUserError()
  endsubroutine

  subroutine IOError(buffer)
    use definitions,only:u6
    implicit none
    character(len=*),intent(in) :: buffer

    call WarningMessage(2,"I/O error when reading line.")
    write(u6,*) "Last line read from input: ",buffer
    call Quit_OnUserError()
  endsubroutine

endmodule
