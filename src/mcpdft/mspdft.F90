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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!***********************************************************************

module mspdft
  use definitions,only:iwp,wp
  implicit none
  private

  character(len=8) :: mspdftmethod = " Unknown"
  integer(kind=iwp) :: iIntS

  real(kind=wp),allocatable :: heff(:,:)

  public :: mspdftmethod,heff
  public :: iIntS

  public :: mspdft_finalize,mspdft_init

contains

  !> @brief Initializes mspdft calculation
  !>
  !>   Loads effective Hamiltonian from disk and initializes gradient variables as needed
  !>
  !> @author Matthew R. Hennefarth
  subroutine mspdft_init()
    use rasscf_global,only:lroots
    use stdalloc,only:mma_allocate
    use mcpdft_input,only:mcpdft_options
    use mspdftgrad,only:mspdftgrad_init

    character(len=18) :: matrix_info

    call mma_allocate(heff,lroots,lroots,label="heff")

    call readmat2('ROT_HAM',matrix_info,heff,size(heff,1),size(heff,2),7,len(matrix_info),'T')
    call determine_method(matrix_info)

    if(mcpdft_options%grad) then
      call mspdftgrad_init()
    endif

  endsubroutine

  !> @brief determines method by which effective Hamiltonian is constructed
  !>
  !>   Sets mspdftmethodand and writes to stdout
  !>
  !> @author Matthew R. Hennefarth
  !>
  !> @param[in] matrix_info string
  subroutine determine_method(matrix_info)
    use definitions,only:u6
    use PrintLevel,only:verbose
    use mcpdft_output,only:iprloc

    character(len=*),intent(IN) :: matrix_info
    character(len=len(mspdftmethod)) :: buffer
    integer(kind=iwp) :: print_level

    print_level = iPrLoc(1)
    buffer = trim(adjustl(matrix_info))

    select case(buffer)
    case("XMS-PDFT","CMS-PDFT","FMS-PDFT","VMS-PDFT")
      mspdftmethod = buffer
      if(print_level > verbose) then
        write(u6,'(6X,A,A)') 'The MS-PDFT method is ',mspdftmethod
      endif
    case default
      if(print_level > verbose) then
        write(u6,'(6X,A)') 'The MS-PDFT calculation is based on a user-supplied rotation matrix'
      endif
    endselect

  endsubroutine

  !> @brief Performs the final parts of a MS-PDFT calculation
  !>
  !>   1. diagonalize heff
  !>   2. print out the final results
  !>   3. save to jobiph if requested
  !>
  !> @author Matthew R. Hennefarth
  !>
  !> @param[in]   nroots number of roots, or dimension of heff
  subroutine mspdft_finalize(nroots)
    use definitions,only:iwp,wp,u6
    use stdalloc,only:mma_allocate,mma_deallocate
    use mcpdft_output,only:iprloc
    use mspdft_util,only:print_effective_ham,print_final_energies,print_mspdft_vectors
    use mcpdft_input,only:mcpdft_options
    use write_pdft_job,only:writejob
    use printlevel,only:usual,terse
    use mspdftgrad,only:mspdftgrad_free

    integer(kind=iwp),intent(in) :: nroots

    real(kind=wp) :: e_mspdft(nroots),si_pdft(nroots,nroots),wgronk(1)
    real(kind=wp),allocatable :: scratch(:)
    integer(kind=iwp) :: info,dim_scratch,iprlev
    character(len=120) :: Line

    iprlev = iprloc(1)

    if(iprlev >= usual) then
      write(u6,*)
      write(Line,'(6X,2A)') MSPDFTMethod,' FINAL RESULTS'
      call CollapseOutput(1,Line)
      write(u6,'(6X,A)') repeat('-',len_trim(Line)-3)
      write(u6,*)

      if(mcpdft_options%otfnal%is_hybrid()) then
        write(u6,'(6X,3A)') 'Hybrid ',MSPDFTMethod,' Effective Hamiltonian'
      else
        write(u6,'(6X,2A)') mspdftmethod,' Effective Hamiltonian'
      endif
      call print_effective_ham(heff,nroots,10)
    endif

    ! Now we diagonalize the final MS-PDFT H-matrix
    ! Eigenvectors will be stored in heff.
    ! Eigenvalues will be stored in e_mspdft

    ! Since the dsyev_ call will override the heff with the orthonormal
    ! eigenvectors, I am going to set the new variable now.
    si_pdft = heff

    ! This first call is to get how big of a scratch we need
    call dsyev_('V','U',nroots,si_pdft,nroots,e_mspdft,wgronk,-1,info)

    dim_scratch = int(wgronk(1))
    call mma_allocate(scratch,dim_scratch,label="XScratch")
    ! Now we actually do the diagonalization
    call dsyev_('V','U',nroots,si_pdft,nroots,e_mspdft,scratch,dim_scratch,info)
    call mma_deallocate(scratch)

    if(iprlev >= terse) call print_final_energies(e_mspdft,nroots)

    ! Update information on the runfile for possible gradient calculations.
    call put_dArray('Last energies',e_mspdft,nroots)
    if(mcpdft_options%grad) then
      call Put_dScalar('Last energy',e_mspdft(mcpdft_options%rlxroot))
    endif

    ! Add info the checkfile for testing!
    call Add_Info("MSPDFTE",e_mspdft,nroots,8)

    if(iprlev >= usual) then
      if(mcpdft_options%otfnal%is_hybrid()) then
        write(u6,'(6X,3A)') 'Hybrid ',MSPDFTMethod,' Eigenvectors:'
      else
        write(u6,'(6X,2A)') MSPDFTMethod,' Eigenvectors:'
      endif
      call print_mspdft_vectors(si_pdft,nroots)
    endif

    ! Added by Chen to write energies and states of MS-PDFT into JOBIPH
    if(mcpdft_options%wjob) call writejob(e_mspdft,nroots,si_pdft=si_pdft)

    if(iprlev >= usual) then
      call CollapseOutput(0,Line)
      write(u6,*)
    endif

    if(mcpdft_options%grad) then
      call put_lscalar('CalcNAC_Opt     ',.false.)
      call put_lscalar('MECI_via_SLAPAF',.false.)
      if(mcpdft_options%nac) then
        call mspdftgrad_misc(si_pdft,mcpdft_options%nac_states)
      else
        call put_iscalar('Relax CASSCF root',mcpdft_options%rlxroot)
        call mspdftgrad_misc(si_pdft,(/mcpdft_options%rlxroot,mcpdft_options%rlxroot/))
      endif
    endif

    ! Deallocate here!
    call mma_deallocate(heff)
    if(mcpdft_options%grad) then
      call mspdftgrad_free()
    endif
  endsubroutine mspdft_finalize

endmodule mspdft
