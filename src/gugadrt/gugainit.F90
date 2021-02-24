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

subroutine gugainit()
! default value for performing ci calculation

use gugadrt_global, only: ludrt, ng_sm, nlsm_all, nlsm_bas, nlsmddel, nlsmedel
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: maxmolcasorb = 5000
integer(kind=iwp) :: i, idisk, idum(1), idx, idx_idisk(64), lenrd, lucimo, luonemo, nbas(8), nc, ncone(64), ndel(8), nfro(8), &
                     norb(8), nsym
character(len=8) :: fncimo, fndrt, fnonemo
real(kind=wp) :: dum(1)
real(kind=wp), allocatable :: cmo(:)
character, allocatable :: bsbl(:)

fnonemo = 'TRAONE'
fndrt = 'CIDRT'
fncimo = 'CIMO'
luonemo = 30
ludrt = 31
lucimo = 32

call daname(lucimo,fncimo)
call daname(ludrt,fndrt)
!call molcas_open(ludrt,fndrt)
! open file traone to read orbital informations
call daname(luonemo,fnonemo)
!call molcas_open(luonemo,fnonemo)
idisk = 0
call idafile(luonemo,2,ncone,64,idisk)
call ddafile(luonemo,2,dum,1,idisk)
!ecor = dum(1)
call idafile(luonemo,2,idum,1,idisk)
nsym = idum(1)
call idafile(luonemo,2,nbas,8,idisk)
call idafile(luonemo,2,norb,8,idisk)
call idafile(luonemo,2,nfro,8,idisk)
call idafile(luonemo,2,ndel,8,idisk)
lenrd = 2*4*maxmolcasorb
call mma_allocate(bsbl,lenrd,label='bsbl')
call cdafile(luonemo,2,bsbl,lenrd,idisk)
nc = 0
do i=1,nsym
  nc = nc+nbas(i)**2
end do
call mma_allocate(cmo,nc,label='cmo')
call ddafile(luonemo,2,cmo,nc,idisk)

idx = 0
call idafile(lucimo,1,idx_idisk,64,idx)
idx_idisk(2) = idx
call cdafile(lucimo,1,bsbl,lenrd,idx)
idx_idisk(3) = idx
call ddafile(lucimo,1,cmo,nc,idx)
idx_idisk(4) = idx
idx = 0
call idafile(lucimo,1,idx_idisk,64,idx)

call mma_deallocate(bsbl)
call mma_deallocate(cmo)

call daclos(lucimo)
call daclos(luonemo)
nlsm_bas(1:8) = nbas(1:8)
nlsmddel(1:8) = nfro(1:8)
nlsmedel(1:8) = ndel(1:8)

!write(u6,'(a4,1x,8(2x,i8))') 'ncon',ncone(1:8)
!write(u6,*) 'idisk : ', idisk
!write(u6,'(a4,1x,f18.9)') 'ecor',ecor
!write(u6,'(a4,1x,i8)') 'nsym',nsym
!write(u6,'(a4,1x,8(2x,i8))') 'nbas',nbas(1:8)
!write(u6,'(a4,1x,8(2x,i8))') 'norb',norb(1:8)
!write(u6,'(a4,1x,8(2x,i8))') 'nfro',nfro(1:8)
!write(u6,'(a4,1x,8(2x,i8))') 'ndel',ndel(1:8)

ng_sm = nsym
nlsm_all(1:8) = norb(1:8)

return

end subroutine gugainit
