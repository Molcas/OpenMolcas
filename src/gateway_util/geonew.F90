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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine GeoNew(lprint)
!***********************************************************************
!                                                                      *
! Object: to pick up the geometry from a special file. This will only  *
!         make any difference of there exist a file otherwise SEWARD   *
!         will use the geometry as specified by the standard input     *
!         file.                                                        *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             March 1991                                               *
!***********************************************************************

use RunFile_procedures, only: Get_Coord_New
use Basis_Info, only: dbsc, nCnttp
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: lprint
integer(kind=iwp) :: iCnt, iCnttp, iDC, iNuc, lBuf, nNuc
logical(kind=iwp) :: Exists
real(kind=wp), allocatable :: CN(:,:)

! Prologue

! Check if there is a data field called 'NewGeom'

call Get_Coord_New(CN,lBuf)

! Quit if the datafield 'NewGeom' is not available. However,
! if the field is available on RUNOLD pick it up there.

if (lBuf == 0) then

  ! Check RUNOLD

  call f_Inquire('RUNOLD',Exists)
  if (Exists) then
    call NameRun('RUNOLD')
    call Get_Coord_New(CN,lBuf)
    if (lBuf == 0) then
      nNuc = 0
      call NameRun('RUNFILE')
      return
    else
      call Get_iScalar('Unique atoms',nNuc)
      call NameRun('RUNFILE')
      if (lprint) then
        write(u6,*)
        write(u6,'(A)') '    Geometry read from RUNOLD'
        write(u6,*)
      end if
    end if
  else
    nNuc = 0
    return
  end if
else
  call Get_iScalar('Unique atoms',nNuc)
  if (lprint) then
    write(u6,*)
    write(u6,'(A)') '    Geometry read from RUNFILE'
    write(u6,*)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Replace coodinates read in subroutine input

iDC = 1
iNuc = 0
!call RecPrt('CN',' ',CN,3,nNuc)
outer: do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      dbsc(iCnttp)%Coor(1:3,iCnt) = CN(1:3,iDC)
      iDC = iDC+1
      iNuc = iNuc+1
      if (iNuc == nNuc) exit outer
    end do
  end if
end do outer
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end

call mma_deallocate(CN)

return

end subroutine GeoNew
