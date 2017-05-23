************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
#ifdef _DMRG_
************************************************************************
*                                                                      *
*  Purpose:  read DMRG input in QCMaquis style                         *
*                                                                      *
************************************************************************
      subroutine socc_dmrg_rdinp(luinput,initial_occ,nrs2t,nroots)

      use qcmaquis_interface_utility_routines, only: lower_to_upper

      implicit none

      integer, intent(in)    :: luinput
      integer, intent(in)    :: nrs2t, nroots
      integer, intent(inout) :: initial_occ(nrs2t,nroots)

      integer                :: i, j, k, io, ls
      character(len=800)     :: string
      character(len=  1)     :: tag

      i = 1
      do
        if(i > nroots) exit

        string(1:800) = ""
        read(luinput,'(a)',IOSTAT=io) string
        if (is_iostat_end(io)) exit
        if (io>0) call SysAbendMsg('socc_dmrg_rdinp',
     &                 'problem reading QCMaquis SOCC input','')
        ls = len_trim(string)
        j = 0
        do k = 1, ls

          tag = string(k:k)
          if(tag == ' ' .or. tag == ',') cycle

          j = j + 1
          if(j > nrs2t) exit

          call lower_to_upper(tag)

          if(tag == 'U')then
            initial_occ(j,i) = 3
          else if(tag == 'D')then
            initial_occ(j,i) = 2
          else if(tag == '2')then
            initial_occ(j,i) = 4
          else if(tag == '0')then
            initial_occ(j,i) = 1
          else
            call SysAbendMsg('socc_dmrg_rdinp',
     &           'problem reading QCMaquis SOCC input:',
     &           'unknown orbital occupation: I know U, D, 2, 0')
          end if
        end do

        if(j /= nrs2t) call SysAbendMsg('socc_dmrg_rdinp',
     &                      'problem reading QCMaquis SOCC input:',
     &                      'wrong number of active orbitals')
        i = i + 1
      end do

      end subroutine socc_dmrg_rdinp
#else
      subroutine socc_dmrg_rdinp()
      implicit none
      end subroutine socc_dmrg_rdinp
#endif
