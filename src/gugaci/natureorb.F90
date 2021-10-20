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

subroutine natureorb(nsbas,nsall,nsdel,ngsm,den1,lden,cmo,lcmo,bsbl,lenb,cno,occ,nmo,pror)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nsbas(8), nsall(8), nsdel(8), ngsm, lden, lcmo, lenb, nmo
real(kind=wp), intent(in) :: den1(lden), cmo(lcmo), pror
real(kind=wp), intent(out) :: cno(lcmo), occ(nmo)
character, intent(in) :: bsbl(lenb)
integer(kind=iwp) :: i, im, j, nc, nc0, nc1, nc2, nc3, nc4, nc5, nsfrz(8)
real(kind=wp) :: val
integer(kind=iwp), allocatable :: nsort(:)
real(kind=wp), allocatable :: buff(:)
character(len=128), parameter :: header = 'MRCISD Natural orbital'

nc0 = 1
do im=1,ngsm
  if (nsall(im) == 0) cycle
  nc1 = nsall(im)*(nsall(im)+1)/2
  nc0 = nc0+nc1
end do
!open(100,file='dat2')
!open(200,file='dat1')
!do i=1,1293
!  read(100,*) den1(i)
!  write(200,'(f14.8)') den1(i)
!end do
!close(100)
!close(200)

call mma_allocate(buff,nmo**2,label='buff')
call mma_allocate(nsort,nmo,label='nsort')

val = 0
cno(:) = cmo(:)
nc0 = 1
nc1 = 1
nc2 = 1
do im=1,ngsm
  if (nsbas(im) == 0) cycle
  nsfrz(im) = nsbas(im)-nsall(im)-nsdel(im)
  ! For density matrix, we only have correlated orbital density matrix
  buff(:) = Zero
  nc = nsall(im)*(nsall(im)+1)/2
  ! Diagonalize density matrix in symmetry block im and transform MO
  ! Copy density matrix
  buff(1:nc) = den1(nc2:nc2+nc-1)
  nc3 = nc1+nsfrz(im)*nsbas(im)
  call jacob(buff,cno(nc3),nsall(im),nsbas(im))
  ! OCC num from diagonal element
  nc3 = nc0
  occ(nc3:nc3+nsfrz(im)-1) = Two
  nc3 = nc3+nsfrz(im)
  do i=1,nsall(im)
    nc = i*(i+1)/2
    occ(nc3) = buff(nc)
    nc3 = nc3+1
  end do

  ! Sort natural orbital in OCC num decreasing order
  ! Sort and copy correlated natural orb
  do i=1,nsall(im)
    nsort(i) = i
  end do
  nc = nc0+nsfrz(im)-1
  do i=1,nsall(im)
    nc3 = i
    do j=i+1,nsall(im)
      if (occ(nc+j) > occ(nc+i)) then
        val = occ(nc+j)
        occ(nc+j) = occ(nc+i)
        occ(nc+i) = val
        nc3 = j
      end if
    end do
    nsort(i) = nc3
  end do

  buff(:) = Zero
  nc3 = nsbas(im)**2
  buff(1:nc3) = cno(nc1:nc1+nc3-1)
  ! Copy sorted orbital
  nc = nc1+nsfrz(im)*nsbas(im)
  nc5 = 1+nsfrz(im)*nsbas(im)
  do i=1,nsall(im)
    j = nsort(i)
    nc3 = nc+(i-1)*nsbas(im)
    nc4 = nc5+(j-1)*nsbas(im)
    cno(nc3:nc3+nsbas(im)-1) = buff(nc4:nc4+nsbas(im)-1)
  end do

  ! Do nothing for deleted orbital
  !nc3 = 1+(nsfrz(im)+nsall(im))*nsbas(im)
  !nc = nsdel(im)*nsbas(im)-1
  !cno(nc3:nc3+nc) = buff(nc3:nc3+nc)

  nc0 = nc0+nsbas(im)
  nc1 = nc1+nsbas(im)**2
  nc2 = nc2+nsall(im)*(nsall(im)+1)/2
end do

call mma_deallocate(buff)
call mma_deallocate(nsort)

!nc0 = 1
!do im=1,ngsm
!  if (nsbas(im) == 0) cycle
!  write(u6,*) ' '
!  write(u6,*) 'SYM',im
!  write(u6,'(10(1xf8.4))') occ(nc0:nc0+nsbas(im)-1)
!  nc0 = nc0+nsbas(im)
!end do

call primo(Header,.true.,.false.,pror,Zero,ngsm,nsbas,nsbas,bsbl,[val],occ,cno,-1)

return

end subroutine natureorb
