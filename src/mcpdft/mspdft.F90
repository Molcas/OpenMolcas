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
  real(kind=wp),allocatable :: DIDA(:,:)
  real(kind=wp),allocatable :: D1SAOMS(:,:)
  real(kind=wp),allocatable,dimension(:,:) :: heff

  public :: mspdftmethod,heff
  public :: iIntS,DIDA,D1SAOMS

  public :: mspdft_finalize,mspdft_init

contains
  subroutine mspdft_init()
    use rasscf_global,only:lroots
    use stdalloc,only:mma_allocate
    use mcpdft_input,only:mcpdft_options
    use mspdftgrad,only:mspdftgrad_init
    implicit none

    call mma_allocate(heff,lroots,lroots,label="heff")
    call load_rotham()

    if(mcpdft_options%grad) then
      call mspdftgrad_init()
    endif

  endsubroutine
  subroutine load_rotham()
    implicit none

    character(len=18) :: matrix_info

    call readmat2('ROT_HAM',matrix_info,heff,size(heff,1),size(heff,2),7,len(matrix_info),'T')
    call determine_method(matrix_info)

  endsubroutine

  subroutine determine_method(matrix_info)
    use definitions,only:u6
    use PrintLevel,only:verbose
    use mcpdft_output,only:iprloc
    implicit none

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

  subroutine mspdft_finalize(nroots,iadr19)
    ! Performs the final parts of MS-PDFT, namely:
    !   1. diagonalize heff
    !   2. print out the final results
    !   3. save to jobiph if requested
    !
    ! Args:
    !   heff: ndarray of length nroots*nroots
    !     MS-PDFT effective hamiltonian matrix. Diagonal elements should
    !     already be replaced with the correct energies.
    !
    !   nroots: integer
    !     number of roots, or dimension of heff
    !
    !   iadr19: integer list of length 15
    !     Holds some information when writing out

    use definitions,only:iwp,wp,u6
    use stdalloc,only:mma_allocate,mma_deallocate
    use mcpdft_output,only:iprloc
    use mspdft_util,only:print_effective_ham,print_final_energies, &
                          print_mspdft_vectors
    use mcpdft_input,only:mcpdft_options
    use write_pdft_job,only:writejob
    use printlevel,only:usual
    use mspdftgrad,only:mspdftgrad_free
    implicit none

    integer(kind=iwp),intent(in) :: nroots
    integer(kind=iwp),dimension(15),intent(in) :: IADR19

    real(kind=wp),dimension(nroots) :: e_mspdft
    real(kind=wp),dimension(nroots,nroots) :: si_pdft
    integer(kind=iwp) :: info,dim_scratch,iprlev
    real(kind=wp),dimension(1) :: wgronk
    real(kind=wp),dimension(:),allocatable :: scratch
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

    if(iprlev >= usual) call print_final_energies(e_mspdft,nroots,mspdftmethod)

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
    if(mcpdft_options%wjob) call writejob(iadr19,e_mspdft=e_mspdft,si_pdft=si_pdft)

    if(iprlev >= usual) then
      call CollapseOutput(0,Line)
      write(u6,*)
    endif

    if(mcpdft_options%grad) then
      if(mcpdft_options%nac) then
        call MSPDFTNAC_Misc(si_pdft)
      else
        call mspdftgrad_misc(si_pdft)
      endif
    endif

    ! Deallocate here!
    call mma_deallocate(heff)
    if(mcpdft_options%grad) then
      call mspdftgrad_free()
    endif
  endsubroutine mspdft_finalize

endmodule mspdft
