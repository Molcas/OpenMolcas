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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine Add_Info(Label,Array,nArray,iTol)
!***********************************************************************
!                                                                      *
!     written by:                                                      *
!     V.Veryazov                                                       *
!                                                                      *
!    parameters:                                                       *
!       Label - text label                                             *
!       Array(nArray) - values to check                                *
!       iTol - tolerance:                                              *
!               if positive: 10**(-iTol)                               *
!               else: -iTol                                            *
!                                                                      *
!***********************************************************************

use Para_Info, only: MyRank
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nArray
real(kind=wp), intent(in) :: Array(nArray)
integer(kind=iwp), intent(in) :: iTol
integer(kind=iwp) :: i, ia, iArray, icomma, iDum(1), iGeoData, iGeoInfo(2), ik, isfound, iuGeoData, j, k, l, length, LuDispEn, n0, &
                     nIntCoord, nlabel, Num
real(kind=wp) :: Aux(1)
logical(kind=iwp) :: Found
character(len=256) :: Collect, STMP, STRING, STRING2
character(len=120) :: Line
character(len=30) :: Junk
character(len=15) :: Energy_File
character(len=8) :: Tol
character(len=*), parameter :: GeoDataF = 'GEODATA'
integer(kind=iwp), external :: isFreeUnit
!integer(kind=iwp) :: irecl
!character(len=32) :: File_Name
!logical(kind=iwp) :: Exist, is_error
interface
  subroutine add_molcas_info(str,n) bind(C,name='add_molcas_info_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: str(*)
    integer(kind=MOLCAS_C_INT) :: n
  end subroutine add_molcas_info
end interface

!------------------------------------------------
! If this is a fake parallel run (e.g. inside the parallel loop of CASPT2_gradient,
! then do not add info - just return immidiately
#ifdef _MOLCAS_MPP_
if ((.not. King()) .and. (.not. Is_Real_Par())) return
#endif
! Num - is a number of exported variables from an array.
Num = 20
!File_Name = 'molcas_info'

!Lu_Info = 99

!---------------------------------------------------------------------*
! Check the file status                                               *
!---------------------------------------------------------------------*
call open_molcas_info()
!call f_Inquire(File_name, Exist)
!----------------------------------------------------------------------*
! Open existing file and position the record pointer to the end        *
!----------------------------------------------------------------------*
!if (Exist) then
!  call molcas_open_Ext2(Lu_info,file_name,'sequential','formatted',ios,.false.,irecl,'unknown',is_error)
!  if (ios /= 0) then
!    write(u6,*) 'Add_Info: cannot create info file'
!    write(u6,*) 'Check file permissions!'
!    call Abend()
!  end if
!  nLines = 0
!  do
!    read(Lu_Info,'(A)',iostat=ios) Line
!    if (ios < 0) exit
!    nLines = nLines+1
!  end do
!  rewind(Lu_Info)
!  do iLine=1,nLines
!    read(Lu_Info,'(A)') Line
!  end do
!# ifdef NAGFOR
!  ! FIXME: ugly hack to make NAG compiler happy
!  close(Lu_Info)
!  open(Lu_Info,file=file_name,position='append')
!# endif
!----------------------------------------------------------------------*
! Open new file                                                        *
!----------------------------------------------------------------------*
!else
!  call molcas_open_Ext2(Lu_info,file_name,'sequential','formatted',ios,.false.,irecl,'unknown',is_error)
!  if (ios /= 0) then
!    write(u6,*) 'Add_Info: cannot create a new file'
!    write(u6,*) 'Check file permissions!'
!    call Abend()
!  end if
!  do i=1,Len(Line)
!    Line(i:i)='#'
!  end do
!  write(Lu_Info,'(A)') Line
!  write(Lu_Info,'(A)') '# MOLCAS-Info_File Vers.No. 1.1'
!  write(Lu_Info,'(A)') Line
!end if
!----------------------------------------------------------------------*
! Append new information                                               *
!----------------------------------------------------------------------*
write(Tol,'(i8)') merge(8,iTol,iTol == 0)
nlabel = len(label)
Line = label
n0 = nlabel
do i=1,nlabel
  if (label(i:i) == ' ') Line(i:i) = '_'
end do
call upcase(Line)
!------------------ Some code for Geo-Environment //Jonas 2011
!----------------------------------------------------------------------*
! Check if we are in the GEO-Environment and append if energy          *
!----------------------------------------------------------------------*
call qpg_iArray('GeoInfo',Found,length)
if (Found) then
  call get_iArray('GeoInfo',iGeoInfo,2)
  if ((nArray == 1) .and. (iGeoInfo(1) == 1) .and. (Label(1:2) == 'E_')) then
    write(Energy_File,'(A,I4.4)') 'disp.energy',iGeoInfo(2)
    LuDispEn = 1
    LuDispEn = isfreeunit(LuDispEn)
    call Molcas_Open(LuDispEn,Energy_File)
    write(LuDispEn,'(F16.8)') Array(1)
    close(LuDispEn)
    iGeoData = 0
    iuGeoData = 10
    iuGeoData = isFreeUnit(iuGeoData)
    call DaName_WA(iuGeoData,GeoDataF)
    call iDaFile(iuGeoData,2,iDum,1,iGeoData)
    nIntCoord = iDum(1)
    iGeoData = iGeoInfo(2)*(nIntCoord+1)+1
    Aux(1) = Array(1) ! because the argument is inout in dDaFile
    call dDaFile(iuGeoData,1,Aux,1,iGeoData)
    call DaClos(iuGeoData)
  end if
end if
isfound = 1
if (MyRank == 0) then
  !---------------------------------------------------------------------

  ! If this label should not be checked, then just return immediately:
  STRING = ' '
  call GETENVF('MOLCAS_NOCHECK',STRING)
  call upcase(STRING)
  STMP = STRING
  isfound = 0
  do
    icomma = index(STMP,',')
    if (icomma /= 0) then
      STRING = ' '
      STRING = STMP(1:icomma-1)
      if (icomma < len(STMP)) then
        STMP = STMP(icomma+1:)
      else
        STMP = ''
      end if
    else
      STRING = STMP
      STMP = ' '
    end if
    k = 0
    l = len(STRING)
    do i=1,l
      if (STRING(i:i) == ' ') then
        if (k > 0) then
          if (STRING2(1:k) == Line(1:k)) then
            isfound = 1
            exit
          end if
        end if
        k = 0
      else
        k = k+1
        STRING2(k:k) = STRING(i:i)
      end if
    end do
    if (STMP == ' ') exit
  end do
end if

if (isfound /= 1) then
  do iArray=1,nArray
    nlabel = n0
    if (nArray /= 1) then
      write(Junk,'(a,i3,a)') '[',iArray-1,']'
      do j=1,5
        if (Junk(j:j) /= ' ') then
          nlabel = nlabel+1
          Line(nlabel:nlabel) = Junk(j:j)
        end if
      end do
    end if
    nlabel = nlabel+1
    Line(nlabel:nlabel) = '='
    nlabel = nlabel+1
    Line(nlabel:nlabel) = '"'
    ia = int(Array(iArray)+0.3_wp)
    if ((abs(Array(iArray)-ia) < 1.0e-7_wp) .and. (ia /= 0)) then
      write(Junk,'(I30)') ia
    else
      if (abs(Array(iArray)) > 1.0e-14_wp) then
        write(Junk,'(F30.12)') Array(iArray)
      else
        Junk = '0.0'
      end if
    end if
    do j=1,30
      if (Junk(j:j) /= ' ') then
        nlabel = nlabel+1
        Line(nlabel:nlabel) = Junk(j:j)
      end if
    end do
    nlabel = nlabel+1
    Line(nlabel:nlabel) = '"'
    if (iArray < Num) then
      !----------------------------------------------------------------*
      ! Export only head of Array                                      *
      !----------------------------------------------------------------*
      !write(Lu_Info,*)
      collect = ''
      collect(1:nlabel) = Line(1:nlabel)
      call add_molcas_info(collect,nlabel)
      !write(Lu_Info,'(a)') Line(1:nlabel)
      if (iArray == nArray) then
        collect = 'export '//Line(1:n0)
        call add_molcas_info(collect,n0+7)
        !write(Lu_Info,'(a,a)') 'export ',Line(1:n0)
      end if
    end if
    ik = 0
    do i=1,8
      if (Tol(i:i) /= ' ') then
        ik = ik+1
        Junk(ik:ik) = Tol(i:i)
      end if
    end do
    collect = '#> '//Line(1:nlabel)//'/'//Junk(1:ik)
    call add_molcas_info(collect,nlabel+ik+4)
    !write(Lu_Info,'(a,a,a,a)') '#> ',Line(1:nlabel),'/',Junk(1:ik)
  end do
end if
!----------------------------------------------------------------------*
! Close file                                                           *
!----------------------------------------------------------------------*
call close_molcas_info()
!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Add_Info
