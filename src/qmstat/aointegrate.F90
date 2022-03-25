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

subroutine AOIntegrate(nBaseQ,nBaseC,Ax,Ay,Az,iQ_Atoms,nAtomsCC,AOint,oV2,N,lmax,Inside)

use qmstat_global, only: CasOri, Cordst, iOrb, iPrint, iQn, nCent, nCnC_C, SavOri, V3
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nBaseQ, nBaseC, iQ_Atoms, nAtomsCC, N, lmax
real(kind=wp) :: Ax, Ay, Az, AOint(nBaseQ,nBaseC), oV2(nBaseC,iOrb(2))
logical(kind=iwp) :: Inside(iQ_Atoms,nAtomsCC)
#include "Molcas.fh"
integer(kind=iwp) :: iBaS, iOrS, j, k, m
real(kind=wp) :: Dummy(1), Dx, Dy, Dz, Rot(3,3), x, y, z
logical(kind=iwp) :: PrEne, PrOcc
character(len=30) :: Snack
real(kind=wp), allocatable :: V2(:,:)
character(len=LenIn8), allocatable :: BsLbl(:)

!----------------------------------------------------------------------*
! Call Transrot. There we compute the rotation matrix for the classical*
! water under consideration. Used later.                               *
!----------------------------------------------------------------------*
call TransRot(Cordst(:,N+1:N+3),N+1,Rot,Dx,Dy,Dz,Ax,Ay,Az)
if (iPrint >= 17) then
  write(u6,*)
  write(u6,*) 'ROTATION MATRIX, Molecule ',N/nCent
  write(u6,*) (Rot(1,k),k=1,3)
  write(u6,*) (Rot(2,k),k=1,3)
  write(u6,*) (Rot(3,k),k=1,3)
end if
!----------------------------------------------------------------------*
! Call OrbRot2. Given the rotation matrix (Rot) and the original MO-   *
! coefficients, we transform them to new MO-coefficients. V2 is on     *
! input the original MO-coefficients (stored in V3), and on output the *
! rotated.                                                             *
!----------------------------------------------------------------------*
call mma_allocate(V2,nBaseC,iOrb(2),label='V2')
do iOrS=1,iOrb(2) !Collect original MO-coeff.
  do iBaS=1,nBaseC
    V2(iBaS,iOrS) = V3(iBaS,iOrS)
  end do
end do
call OrbRot2(Rot,V2,iQn,iOrb(2),nBaseC,lMax,nCnC_C)
! Store the rotated in vector for later convenience.
oV2(:,:) = V2
if (iPrint >= 25) then !Optional print-out.
  PrOcc = .false.
  PrEne = .false.
  write(snack,'(A,I3)') 'Rotated orbitals for water ',N/nCent
  call mma_allocate(BsLbl,nBaseC,label='BsLbl')
  call NameRun('WRUNFIL')
  call Get_cArray('Unique Basis Names',BsLbl,LenIn8*nBaseC)
  call Primo(Snack,PrOcc,PrEne,Dummy(1),Dummy(1),1,[nBaseC],iOrb(2),BsLbl,Dummy,Dummy,V2,3)
  call mma_deallocate(BsLbl)
end if
call mma_deallocate(V2)
do m=1,lMax !New basis function origo defined.
  x = Zero
  y = Zero
  z = Zero
  do j=1,3
    x = x+Rot(1,j)*SavOri(j,m)
    y = y+Rot(2,j)*SavOri(j,m)
    z = z+Rot(3,j)*SavOri(j,m)
  end do
  CasOri(1,m) = x+Dx
  CasOri(2,m) = y+Dy
  CasOri(3,m) = z+Dz
end do
!----------------------------------------------------------------------*
! Compute overlap between the contracted basis functions on the water  *
! molecule presently studied and the QM-molecule.                      *
!----------------------------------------------------------------------*
AOInt(:,:) = Zero
call ContractOvl(AOint,nBaseQ,nBaseC,N,nCent,iQ_Atoms,nAtomsCC,iPrint,Inside)

return

end subroutine AOIntegrate
