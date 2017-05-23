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
      subroutine read_abdata
      implicit none
#include "SysDef.fh"
#include "abtab.fh"
      character(*), parameter :: ABDATA_NAME = 'ABDATA'
      integer, parameter :: lu_abdata = 22
      logical :: found_abdata
*
      character(8) :: key
      integer :: i, itab, ipos, k, nerr
*
      call f_Inquire(ABDATA_NAME,found_abdata)
      if (.not.found_abdata) then
        call warningmessage(2,
     &              ' the abdata.ascii file does not exist.')
        call abend()
      end if
      call molcas_open(lu_abdata,ABDATA_NAME)

  10  read(lu_abdata,'(a8)') key
      if(key.ne.'NTAB1, N') goto 10
      read(lu_abdata,*) ntab1,ntab2,maxdeg
      nerr=0
      if(ntab2-ntab1+1.gt.mxsiz2) then
        call warningmessage(2,' mxsiz2 is too small in readab.')
        write(6,*)' recompile. needs mxsiz2=',ntab2-ntab1+1
        nerr=1
      end if
      if(maxdeg.gt.mxsiz1) then
        call warningmessage(2,' mxsiz1 is too small in readab.')
        write(6,*)' recompile. needs mxsiz1=',maxdeg
        nerr=1
      end if
      if(nerr.eq.1) call abend()
      ipos=0
      do i=ntab1,ntab2
  20    read(lu_abdata,'(a8)') key
        if(key.ne.'TAB POIN') goto 20
        ipos=ipos+1
        read(lu_abdata,*) itab, tvalue(ipos), p0(ipos)
        read(lu_abdata,*)
        read(lu_abdata,*)(atab(k,ipos),k=0,maxdeg)
        read(lu_abdata,*)
        read(lu_abdata,*)(btab(k,ipos),k=0,maxdeg)
      end do
*
      close (lu_abdata)
      end
