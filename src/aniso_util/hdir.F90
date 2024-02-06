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

subroutine hdir(nDir,nDirZee,dirX,dirY,dirZ,dir_weight,nP,nsymm,ngrid,nDirTot,dHX,dHY,dHZ,dHW)
! this routine generates the directions of the applied magnetic
! field according to Lebedev-Laikov grids using the the given parameters (nsymm, ngrid)

use Constants, only: Zero
use Definitions, only: u6

implicit none
integer :: nP, nDirTot, nDir, nDirZee, i, j
integer :: nsymm, ngrid
real(kind=8) :: dirX(nDir), dirY(nDir), dirZ(nDir), dir_weight(nDirZee,3)
real(kind=8) :: dHX(nDirTot), dHY(nDirTot), dHZ(nDirTot), dHW(nDirTot)
real(kind=8) :: X(nP), Y(nP), Z(nP), W(nP)

if ((nDirTot-nDir-nDirZee-nP) /= 0) then
  write(u6,'(A   )') 'the number of directions of applied magnetic field is not consistent:'
  write(u6,'(A,i4)') 'nDir    = ',nDir
  write(u6,'(A,i4)') 'nDirZee = ',nDirZee
  write(u6,'(A,i4)') 'nP      = ',nP
  write(u6,'(A,i4)') 'nDirTot = ',nDirTot
  write(u6,'(A,i4)') 'The rule is :'
  write(u6,'(A   )') 'nDir + nDirZee + nP = nDirTot'
  call xFlush(u6)
  call abend()
end if
! intialization
call dcopy_(nDirTot,[Zero],0,dHX(1),1)
call dcopy_(nDirTot,[Zero],0,dHY(1),1)
call dcopy_(nDirTot,[Zero],0,dHZ(1),1)
call dcopy_(nDirTot,[Zero],0,dHW(1),1)
call dcopy_(nP,[Zero],0,X(1),1)
call dcopy_(nP,[Zero],0,Y(1),1)
call dcopy_(nP,[Zero],0,Z(1),1)
call dcopy_(nP,[Zero],0,W(1),1)

!if (nDir > 0) then
!  call DCOPY_(nDir,dirX,1,dHX(1),1)
!  call DCOPY_(nDir,dirY,1,dHY(1),1)
!  call DCOPY_(nDir,dirZ,1,dHZ(1),1)
!end if

!if (nDirZee > 0) then
!  call DCOPY_(nDirZee,dir_weight(1:3,1),1,dHX(1+nDir),1)
!  call DCOPY_(nDirZee,dir_weight(1:3,2),1,dHY(1+nDir),1)
!  call DCOPY_(nDirZee,dir_weight(1:3,3),1,dHZ(1+nDir),1)
!end if

do i=1,nDir
  dHX(i) = dirX(i)
  dHY(i) = dirY(i)
  dHZ(i) = dirZ(i)
end do

do i=1,nDirZee
  dHX(i+nDir) = dir_weight(i,1)
  dHY(i+nDir) = dir_weight(i,2)
  dHZ(i+nDir) = dir_weight(i,3)
end do

call Lebedev_Laikov(nSymm,nGrid,nP,X,Y,Z,W)

do i=1,nP
  j = i+nDir+nDirZee
  dHX(j) = X(i)
  dHY(j) = Y(i)
  dHZ(j) = Z(i)
  dHW(j) = W(i)
end do

!call DCOPY_(nP,X(1),1,dHX(1+nDir+nDirZee),1)
!call DCOPY_(nP,Y(1),1,dHY(1+nDir+nDirZee),1)
!call DCOPY_(nP,Z(1),1,dHZ(1+nDir+nDirZee),1)
!call DCOPY_(nP,W(1),1,dHW(1+nDir+nDirZee),1)

return

end subroutine hdir
