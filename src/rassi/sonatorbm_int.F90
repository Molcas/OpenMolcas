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

subroutine SONATORBM_INT(DENS,CHARPROP,IC_,CHARTYPE,ASS,BSS,iOpt,ROTMAT,PROPVALXR,PROPVALYR,PROPVALZR,PROPVALXI,PROPVALYI,PROPVALZI)

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use rassi_data, only: NBTRI
use Constants, only: Zero
use Definitions, only: u6

implicit none
real*8 DENS(6,NBTRI)
character(len=*) CHARPROP
integer IC_
character(len=*) CHARTYPE
integer ASS, BSS, iOpt
real*8 ROTMAT(3,3), PROPVALXR, PROPVALYR, PROPVALZR, PROPVALXI, PROPVALYI, PROPVALZI
integer IDUM(1)
real*8, allocatable :: IP(:), IPX(:), IPY(:), IPZ(:)
integer ITYPE, NIP, IC_End, IC_Str, IC, JOPT, ICMP, IRC, I, ISCHK

! NOW DO INTEGRATION WITH AO MATRICES
! FOR THE EXPECTATION VALUE

! Get the proper type of the property
ITYPE = 0
if (CHARTYPE == 'HERMSING') ITYPE = 1
if (CHARTYPE == 'ANTISING') ITYPE = 2
if (CHARTYPE == 'HERMTRIP') ITYPE = 3
if (CHARTYPE == 'ANTITRIP') ITYPE = 4
if (ITYPE == 0) then
  write(u6,*) 'RASSI/SONATORB internal error.'
  write(u6,*) 'Erroneous property type:',trim(CHARTYPE)
  call ABEND()
end if

! ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
! The extra 4 elements correspond to the nuclear contribution
! and the origin of the operator
NIP = 4+NBTRI
call mma_allocate(IP,NIP,Label='IP')
if (iOpt == 1) then
  call mma_allocate(IPX,NIP,Label='IPX')
  call mma_allocate(IPY,NIP,Label='IPY')
  call mma_allocate(IPZ,NIP,Label='IPZ')
  IPX(:) = Zero
  IPY(:) = Zero
  IPZ(:) = Zero
end if

if (iOpt == 1) then
  IC_End = 3
  IC_Str = 1
else
  IC_End = IC_
  IC_Str = IC_
end if
do IC=IC_Str,IC_End ! loop over reading X,Y, and Z AO Integrals

  ! Get info from the stored integrals
  ! JOPT controls what is read.
  ! JOPT=1 Read the size information
  ! JOPT=0 Read the property
  ! JOPT=6 Read the property, skipping the nuclear contribution and the origin
  ! (see OneDat module)
  JOPT = ibset(0,sOpSiz)
  ICMP = IC
  call iRDONE(IRC,JOPT,CHARPROP,ICMP,IDUM,ISCHK)

  ! Actually read the integral
  JOPT = 0
  call RDONE(IRC,JOPT,CHARPROP,ICMP,IP,ISCHK)

  if (IRC /= 0) then
    write(u6,*)
    write(u6,'(6X,A)') '*** ERROR IN SUBROUTINE SONATORB ***'
    write(u6,'(6X,A)') '  FAILED IN READING FROM  ONEINT'
    write(u6,'(6X,A,A)') '  LABEL     = ',trim(CHARPROP)
    write(u6,'(6X,A,I2)') '  COMPONENT = ',IC
    write(u6,*)
    call ABEND()
  end if

  if (iOpt == 1) then
    ! note reordering
    call DAXPY_(NIP,ROTMAT(IC,1),IP,1,IPX,1)
    call DAXPY_(NIP,ROTMAT(IC,2),IP,1,IPY,1)
    call DAXPY_(NIP,ROTMAT(IC,3),IP,1,IPZ,1)
  end if

end do ! end loop over reading X,Y, and Z AO Integrals

PROPVALXR = Zero
PROPVALYR = Zero
PROPVALZR = Zero
PROPVALXI = Zero
PROPVALYI = Zero
PROPVALZI = Zero

! The integral is NBTRI matrix
! The property is NBTRI matrix
! We only work with half the matrix. Therefore, this would
! have a factor of 2. However, the factor of 1/2 was missing
! in SONATORB.F from the symmetric/antisymmetric equations
if ((ITYPE == 1) .or. (ITYPE == 3)) then
  if (iOpt == 1) then
    do I=1,NBTRI
      PROPVALXR = PROPVALXR+IPX(I)*DENS(1,I)
      PROPVALYR = PROPVALYR+IPY(I)*DENS(2,I)
      PROPVALZR = PROPVALZR+IPZ(I)*DENS(3,I)

      PROPVALXI = PROPVALXI+IPX(I)*DENS(4,I)
      PROPVALYI = PROPVALYI+IPY(I)*DENS(5,I)
      PROPVALZI = PROPVALZI+IPZ(I)*DENS(6,I)
    end do
  else
    do I=1,NBTRI
      PROPVALXR = PROPVALXR+IP(I)*DENS(1,I)
      PROPVALYR = PROPVALYR+IP(I)*DENS(2,I)
      PROPVALZR = PROPVALZR+IP(I)*DENS(3,I)

      PROPVALXI = PROPVALXI+IP(I)*DENS(4,I)
      PROPVALYI = PROPVALYI+IP(I)*DENS(5,I)
      PROPVALZI = PROPVALZI+IP(I)*DENS(6,I)
    end do
  end if
else
  if (iOpt == 1) then
    do I=1,NBTRI
      PROPVALXI = PROPVALXI+IPX(I)*DENS(1,I)
      PROPVALYI = PROPVALYI+IPY(I)*DENS(2,I)
      PROPVALZI = PROPVALZI+IPZ(I)*DENS(3,I)

      PROPVALXR = PROPVALXR-IPX(I)*DENS(4,I)
      PROPVALYR = PROPVALYR-IPY(I)*DENS(5,I)
      PROPVALZR = PROPVALZR-IPZ(I)*DENS(6,I)
    end do
  else
    do I=1,NBTRI
      PROPVALXI = PROPVALXI+IP(I)*DENS(1,I)
      PROPVALYI = PROPVALYI+IP(I)*DENS(2,I)
      PROPVALZI = PROPVALZI+IP(I)*DENS(3,I)

      PROPVALXR = PROPVALXR-IP(I)*DENS(4,I)
      PROPVALYR = PROPVALYR-IP(I)*DENS(5,I)
      PROPVALZR = PROPVALZR-IP(I)*DENS(6,I)
    end do
  end if
end if

write(u6,*)
write(u6,*) '************************************'
write(u6,*) 'SONATORB EXPECTATION VALUES'
write(u6,*) ' PROPERTY: ',trim(CHARPROP)
write(u6,*) ' TYPE: ',trim(CHARTYPE)
write(u6,*) ' STATE (K,L): ',ASS,BSS
write(u6,*) '************************************'
write(u6,*) 'Property: Re(X): ',PROPVALXR
write(u6,*) 'Property: Re(Y): ',PROPVALYR
write(u6,*) 'Property: Re(Z): ',PROPVALZR
write(u6,*) 'Property: Im(X): ',PROPVALXI
write(u6,*) 'Property: Im(Y): ',PROPVALYI
write(u6,*) 'Property: Im(Z): ',PROPVALZI
write(u6,*) '************************************'

! Free up un-needed space
call mma_deallocate(IP)
if (iOpt == 1) then
  call mma_deallocate(IPX)
  call mma_deallocate(IPY)
  call mma_deallocate(IPZ)
end if

end subroutine SONATORBM_INT

