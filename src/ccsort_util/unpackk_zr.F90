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

subroutine unpackk_zr(i,vint,ndimv1,ndimv2,ndimv3,key)
! this routine expands integrals packed in i-th TEMP file
! to vint(j,k,l) = <i,j|k,l>
!
! i      - value of pivot index (I)
! vint   - array of integrals (O)
! ndimv1 - first dimension of vint (norb(symj)) (I)
! ndimv2 - second dimension of vint (norb(symk)) (I)
! ndimv3 - third dimension of vint (norb(syml)) (I)
! key    - reduced storing key (I)
!          = 0 if symj is not syml
!          = 1 if symj = syml

use ccsort_global, only: iokey, jh, kh, lh, lrectemp, lunpublic, nrectemp, nsize, tmpnam, valh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: i, ndimv1, ndimv2, ndimv3, key
real(kind=wp), intent(out) :: vint(ndimv1,ndimv2,ndimv3)
integer(kind=iwp) :: daddr, ihelp, ires, length, nhelp, nrec
integer(kind=iwp), allocatable :: iBuf(:)
integer(kind=iwp), parameter :: constj = 1024**2, constk = 1024

! set vint=0

vint(:,:,:) = Zero

! open corresponding TEMP file

if (iokey == 1) then
  ! Fortran IO
  call molcas_binaryopen_vanilla(lunpublic,tmpnam(i))
  !open(unit=lunpublic,file=tmpnam(i),form='unformatted')
else
  ! MOLCAS IO
  call daname(lunpublic,tmpnam(i))
  daddr = 0
end if

call mma_allocate(iBuf,nsize)

do nrec=1,nrectemp(i)

  if (nrec /= nrectemp(i)) then
    length = nsize
  else
    length = lrectemp(i)
  end if

  if (iokey == 1) then
    ! Fortran IO
    read(lunpublic) valh(1:length),iBuf(1:length)
  else
    ! MOLCAS IO
    call ddafile(lunpublic,2,valh,length,daddr)
    call idafile(lunpublic,2,iBuf,length,daddr)
  end if

  ! get indexes jh,kh,lh and value valh from packed form

  do nhelp=1,length
    ihelp = iBuf(nhelp)
    jh(nhelp) = int(ihelp/constj,kind=kind(jh))
    ires = ihelp-constj*jh(nhelp)
    kh(nhelp) = int(ires/constk,kind=kind(jh))
    lh(nhelp) = int(ires-constk*kh(nhelp),kind=kind(lh))
  end do

  if (key == 0) then
    do nhelp=1,length
      vint(jh(nhelp),kh(nhelp),lh(nhelp)) = valh(nhelp)
    end do
  else
    do nhelp=1,length
      vint(jh(nhelp),kh(nhelp),lh(nhelp)) = valh(nhelp)
      vint(lh(nhelp),kh(nhelp),jh(nhelp)) = valh(nhelp)
    end do
  end if

end do

call mma_deallocate(iBuf)

if (iokey == 1) then
  ! Fortran IO
  close(lunpublic)
else
  ! Molcas IO
  call daclos(lunpublic)
end if

return

end subroutine unpackk_zr
