!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!
!!  dmrg-interface-utils is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  dmrg-interface-utils is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with dmrg-interface-utils. If not, see <http://www.gnu.org/licenses/>.

module qcmaquis_interface_measurements

  implicit none

  public compute_observables
  public extract_rdm
  public overlap_kernel

contains

! =================================================================================
!   transfer the measuring scripts (python) for later usage
! ---------------------------------------------------------------------------------
! Input  :  host_program driver
! Output :  QCMaquis measurements stored on the result*.h5 file
! =================================================================================
  subroutine compute_observables(pydriver)

    character(len=300), intent(in) :: pydriver

    call system(trim(pydriver)//" --executable=dmrg_meas"//' --tmpfull=$PWD/tmp'// &
                                " --output=dmrg.RDM_OUT "//" dmrg-input"           &
               )

  end subroutine compute_observables

! ======================================================================
!   Calculating 1-Density matrix, 2-Density matrix, and Spin-denstry
! ----------------------------------------------------------------------
! Input  :  maquis_model
!        :        result     maquis_name_results(iroot)   --  string
!        :    len_result                 length (integer)
!        :           spd     maquis_name_spin-density(iroot)
!        :       len_spd
!        :           DM1     maquis_name_1RDM(iroot)
!        :       len_DM1
!        :           DM2     maquis_name_2RDM(iroot)
!        :       len_DM2
!        :        CIONLY     only DMRG calculation, no orb-opt
! Output :  oneparticle.spd (then it will be renamed basing on i-th root)
!        :  oneparticle.rdm (as above)
!        :  twoparticle.rdm (as above)
! ======================================================================
  subroutine extract_rdm(result,len_result,DM1,len_DM1,DM2,len_DM2,SPD,len_SPD,cionly)

    integer                    :: len_result,len_DM1,len_DM2,len_spd
    character(len=len_DM1)     :: DM1
    character(len=len_DM2)     :: DM2
    character(len=len_SPD)     :: SPD
    character(len=len_result)  :: result
    logical                    :: cionly

    if(cionly)then
      call system("$MOLCAS/pytools/rdmsave_su2_onlyone.py "//result//" ")
    else
      call system("$MOLCAS/pytools/rdmsave_su2.py "//result//" ")
    end if
    call system("$MOLCAS/pytools/spdmsave_su2.py "//result//" ")

                    call system("mv "//"oneparticle.rdm "//DM1)
                    call system("mv "//"oneparticle.spd "//SPD)
    if(.not.cionly) call system("mv "//"twoparticle.rdm "//DM2)

  end subroutine extract_rdm

  subroutine overlap_kernel(pydriver,overlap_exe,maquis_name_states,overlap)

    character(len=300), intent(in)    :: pydriver
    character(len=300), intent(in)    :: overlap_exe
    character(len=2300),intent(inout) :: maquis_name_states(2)
    real*8,             intent(inout) :: overlap

    character(len=100)                :: overlap_between_states
    character(len=100)                :: overlap_string
    character(len=100)                :: overlap_trimmed
    integer                           :: lunit

    overlap_between_states = ""
    overlap_string         = ""
    lunit                  = 140

    call system(trim(pydriver)//' --overlap '//                          &
                " --executable="//trim(overlap_exe)//                    &
                " --output=overlap.DMRG.out "//' --tmpfull=$PWD/tmp'//   &
                ' --lhs='//trim(maquis_name_states(1))//' --rhs='//      &
                trim(maquis_name_states(2))                              &
               )

    !> read data from file and store result in output variable
    open(unit=lunit,form='formatted',status='old',action="read",position='rewind',file="overlap.DMRG.out")
    read(lunit,'(a)') overlap_between_states
    read(lunit,'(a)') overlap_string
    close(lunit,status="keep")

    overlap_trimmed = trim(overlap_string)

    overlap_trimmed = overlap_trimmed(index(overlap_trimmed,'=')+1:100)
    read(overlap_trimmed,*) overlap


#ifdef _DMRG_DEBUG_
    !> debug information
    write(6,*) '  MPS1:',trim(maquis_name_states(1))
    write(6,*) '  MPS2:',trim(maquis_name_states(2))

    write(6,*) overlap_between_states,overlap
    !> end of debug information
#endif

  end subroutine overlap_kernel

end module qcmaquis_interface_measurements

