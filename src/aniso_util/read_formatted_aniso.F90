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

subroutine read_formatted_aniso(input_file_name,nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, gElectron
use Definitions, only: wp, iwp, u6

implicit none
character(len=180), intent(in) :: input_file_name
integer(kind=iwp), intent(inout) :: nss, nstate
integer(kind=iwp), intent(out) :: multiplicity(nstate)
real(kind=wp), intent(out) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate), edmom(3,nstate,nstate), amfi(3,nstate,nstate)
complex(kind=wp), intent(out) :: U(nss,nss), MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss), DM(3,nss,nss), HSO(nss,nss)
integer(kind=iwp) :: j, j1, j2, l, LuAniso
real(kind=wp), allocatable :: tmpI(:,:), tmpR(:,:)
real(kind=wp), parameter :: g_e = -gElectron
integer(kind=iwp), external :: IsFreeUnit

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'Entering read_formatted_aniso'
call xFlush(u6)
#endif
! set to zero all arrays:
multiplicity(:) = 0
eso(:) = Zero
esfs(:) = Zero
angmom(:,:,:) = Zero
edmom(:,:,:) = Zero
amfi(:,:,:) = Zero
U(:,:) = cZero
MM(:,:,:) = cZero
MS(:,:,:) = cZero
ML(:,:,:) = cZero
DM(:,:,:) = cZero
HSO(:,:) = cZero
! read the file "aniso.input":
LuAniso = IsFreeUnit(81)
call molcas_open(LuAniso,trim(input_file_name))
! compatibility with the present version: of aniso_i.input file
read(LuAniso,*) nstate,nss
#ifdef _DEBUGPRINT_
write(u6,'(A,2I6)') 'nstate, nss:',nstate,nss
call xFlush(u6)
#endif
read(LuAniso,*) (eso(j),j=1,nss)
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'ESO:'
write(u6,'(5ES24.14)') (eso(j),j=1,nss)
#endif
read(LuAniso,*) (multiplicity(j),j=1,nstate)
#ifdef _DEBUGPRINT_
write(u6,'(A)') '(multiplicity(j),j=1,nstate)'
write(u6,'(50I3)') (multiplicity(j),j=1,nstate)
call xFlush(u6)
#endif

call mma_allocate(tmpR,nss,nss,'tmpR')
call mma_allocate(tmpI,nss,nss,'tmpI')
! magnetic moment
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  MM(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
  call xFlush(u6)
end do

! spin moment
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  MS(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

! spin-free energies
read(LuAniso,*) (esfs(j),j=1,nstate)
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'ESFS:'
write(u6,'(5ES24.14)') (esfs(j),j=1,nstate)
#endif

! U matrix
do j1=1,nss
  read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
end do

U(1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)

! angmom
do l=1,3
  do j1=1,nstate
    read(LuAniso,*) (angmom(l,j1,j2),j2=1,nstate)
  end do
end do

! compute the orbital moment
ML(:,1:nss,1:nss) = -MM(:,1:nss,1:nss)-g_e*MS(:,1:nss,1:nss)

! DMmom
do l=1,3
  do j1=1,nss
    read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
  end do
  DM(l,1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)
end do

! edmom
do l=1,3
  do j1=1,nstate
    read(LuAniso,*) (edmom(l,j1,j2),j2=1,nstate)
  end do
end do

! amfi
do l=1,3
  do j1=1,nstate
    read(LuAniso,*) (amfi(l,j1,j2),j2=1,nstate)
  end do
end do

! HSO matrix
do j1=1,nss
  read(LuAniso,*) (tmpR(j1,j2),tmpI(j1,j2),j2=1,nss)
end do

HSO(1:nss,1:nss) = cmplx(tmpR(:,:),tmpI(:,:),kind=wp)

call mma_deallocate(tmpR)
call mma_deallocate(tmpI)

close(LuAniso)

return

end subroutine read_formatted_aniso
