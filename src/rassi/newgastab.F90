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

subroutine NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM,ICASE)

use rassi_global_arrays, only: REST1, REST2
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSYM, NGAS, NGASORB(NSYM,NGAS), NGASLIM(2,NGAS), ICASE
integer(kind=iwp) :: IGAS, ISUM, ISYM, ITYPE, KORB, KREST, LPOS, NSIZE
integer(kind=iwp), pointer :: REST(:)

NSIZE = 4+(NGAS+1)*(NSYM+1)+2*NGAS
ITYPE = 91
select case (ICASE)
  case (1)
    call mma_allocate(REST1,NSIZE,Label='REST1')
    REST => REST1(:)
  case (2)
    call mma_allocate(REST2,NSIZE,Label='REST2')
    REST => REST2(:)
  case default
    write(u6,*) 'NEWGASTAB: Illegal ICASE value'
    write(u6,*) 'ICASE=',ICASE
    call Abend()
    REST => REST1(:) !dummy
end select
REST(1) = NSIZE
REST(2) = ITYPE
REST(3) = NGAS
REST(4) = NSYM
KORB = 5
!TEST write(u6,*) ' In NEWGASTAB. NGASORB array is:'
!TEST do igas=1,ngas
!TEST   write(u6,'(1x,8i5)') (ngasorb(isym,igas),isym=1,nsym)
!TEST end do
!TEST write(u6,*) ' In NEWGASTAB. NGASLIM array is:'
!TEST write(u6,'(1x,20i3)') (ngaslim(1,igas),igas=1,ngas)
!TEST write(u6,'(1x,20i3)') (ngaslim(2,igas),igas=1,ngas)

do IGAS=1,NGAS
  ISUM = 0
  do ISYM=1,NSYM
    LPOS = KORB+ISYM+(NSYM+1)*IGAS
    REST(LPOS) = 2*NGASORB(ISYM,IGAS)
    ISUM = ISUM+REST(LPOS)
  end do
  LPOS = KORB+0+(NSYM+1)*IGAS
  REST(LPOS) = ISUM
end do
do ISYM=0,NSYM
  ISUM = 0
  do IGAS=1,NGAS
    LPOS = KORB+ISYM+(NSYM+1)*IGAS
    ISUM = ISUM+REST(KORB+ISYM+(NSYM+1)*IGAS)
  end do
  LPOS = KORB+ISYM
  REST(LPOS) = ISUM
end do
KREST = KORB+(NGAS+1)*(NSYM+1)
REST(KREST:KREST+2*NGAS-1) = pack(NGASLIM,.true.)

nullify(REST)

end subroutine NEWGASTAB
