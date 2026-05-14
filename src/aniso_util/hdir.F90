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

use Lebedev_quadrature, only: ld_by_rule
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDir, nDirZee, nP, nsymm, ngrid, nDirTot
real(kind=wp), intent(in) :: dirX(nDir), dirY(nDir), dirZ(nDir), dir_weight(nDirZee,3)
real(kind=wp), intent(out) :: dHX(nDirTot), dHY(nDirTot), dHZ(nDirTot), dHW(nDirTot)
real(kind=wp) :: X(nP), Y(nP), Z(nP), W(nP)

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

dHX(1:nDir) = dirX(:)
dHY(1:nDir) = dirY(:)
dHZ(1:nDir) = dirZ(:)

dHX(nDir+1:nDir+nDirZee) = dir_weight(:,1)
dHY(nDir+1:nDir+nDirZee) = dir_weight(:,2)
dHZ(nDir+1:nDir+nDirZee) = dir_weight(:,3)

dHW(1:nDir+nDirZee) = Zero

call ld_by_rule(nSymm,nGrid,X,Y,Z,W)

dHX(nDir+nDirZee+1:nDir+nDirZee+nP) = X(:)
dHY(nDir+nDirZee+1:nDir+nDirZee+nP) = Y(:)
dHZ(nDir+nDirZee+1:nDir+nDirZee+nP) = Z(:)
dHW(nDir+nDirZee+1:nDir+nDirZee+nP) = W(:)

return

end subroutine hdir
