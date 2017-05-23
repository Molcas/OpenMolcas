!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
subroutine bcast_2RDM(InFile)

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        function symlink(master, slave) result(ret) bind(c)
            use, intrinsic :: iso_c_binding
            character(c_char), dimension(*), intent(in) :: master, slave
            integer(c_int) :: ret
        end function
    end interface

    character(*) :: InFile
    character(1024) :: master
    integer(c_int) :: lmaster1

    Call prgmtranslate_master(InFile,master,lmaster1)
    write(6,*) 'Master File Name:'
    write(6,*) master(1:lmaster1)
    write(6,*) 'Slave File Name:'
    write(6,*) InFile

    call make_symlink(master,InFile)

contains

! This routine is used to make files generated into the master node available into all the slaves
! nodes via symbolic linking. The way it does it is by translating the names of the adequate paths.



   subroutine make_symlink(master, slave)

        integer(c_int) :: ret
        character(*),intent(in) :: master, slave

        ret = symlink(trim(master) // c_null_char, &
                      trim(slave) // c_null_char)

        if (ret == 0) then
            write(6,*) 'Sym link made successfully'
        else
            write(6,*) 'Symlinking failed', ret
        end if

      end subroutine

end subroutine

