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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine VECSORT(NSYM,NBAS,NORB,CMO,OCC,INDT,NNWORD,NEWORD,iErr)
! Sorting routine: sort CMO, OCC according to INDT

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSYM, NBAS(NSYM), NORB(NSYM), NNWORD
real(kind=wp), intent(inout) :: CMO(*), OCC(*)
integer(kind=iwp), intent(inout) :: INDT(*), iErr
! PAM 2012: Have VecSort return an array NEWORD with new orbital indices.
! If typedef info in an orbital file is used to change the order
! of orbitals, this is done by a call to VecSort. If user, as a result,
! needs to alter orbital indices ( e.g. in the supersymmetry array
! IXSYM) he needs to know how orbitals have changed order.
! The mapping is (New orbital index)=NEWORD(Old orbital index).
integer(kind=iwp), intent(out) :: NEWORD(NNWORD)
! PAM 2012: End of update
integer(kind=iwp) :: i, ii, iii, ip, isfirst, ISYM, iType, j, kcmo, kocc, m, MBAS, NeedSort, nw
real(kind=wp) :: q
real(kind=wp), allocatable :: TCMO(:)

MBAS = NBAS(1)
do i=2,nSym
  MBAS = max(MBAS,NBAS(i))
end do
call mma_allocate(TCMO,MBAS,label='TCMO')

kcmo = 0
kocc = 0
iii = 0

! PAM 2012: NewOrd update
! If NNWORD is > 0, this indicates the caller wish to get back
! a reindicing array. Then this must be large enough:
if (NNWORD > 0) then
  nw = 0
  do ISYM=1,NSYM
    do i=1,NORB(ISYM)
      nw = nw+1
    end do
  end do
  if (nw > NNWORD) then
    call Abend()
  end if

  do nw=1,NNWORD
    NEWORD(nw) = nw
  end do
end if
! PAM 2012: End of update

do ISYM=1,NSYM
  ! Check Do we need make sort?
  NeedSort = 0
  !write(u6,*) 'indt'
  !write(u6,*) (indt(i+iii),i=1,norb(isym))
  do I=1,NORB(ISYM)
    if (IndT(i+iii) == 0) then
      iErr = 1
    end if
    if (i > 1) then
      if (IndT(i+iii) < IndT(i-1+iii)) NeedSort = 1
    end if
  end do
  !write(u6,*) 'NeedSort=',NeedSort
  if (NeedSort == 1) then
    ! Start sort
    ! we do have only a few types of orbitals, so sorting is a simple...
    do iType=1,7
      ip = 0
      isfirst = 0
      do i=1,NORB(ISYM)
        if (isfirst == 0) then
          if (IndT(i+iii) > iType) then
            isfirst = 1
          end if
          if (IndT(i+iii) <= iType) then
            ip = i
          end if
        end if
        if ((isfirst == 1) .and. (IndT(i+iii) == iType)) then
          ! We need to shift CMO, Occ
          m = IndT(i+iii)
          q = Occ(i+KOCC)
          ! PAM 2012: NewOrd update
          if (NNWORD > 0) nw = NEWORD(i+KOCC)
          ! PAM 2012: End of update
          do ii=1,NBAS(ISYM)
            TCMO(ii) = CMO((i-1)*NORB(ISYM)+KCMO+ii)
          end do
          do j=i,ip+2,-1
            IndT(j+iii) = IndT(j-1+iii)
            Occ(j+KOCC) = Occ(j-1+KOCC)
            ! PAM 2012: NewOrd update
            if (NNWORD > 0) NEWORD(j+KOCC) = NEWORD(j-1+KOCC)
            ! PAM 2012: End of update
            do ii=1,NBAS(ISYM)
              CMO((j-1)*NORB(ISYM)+KCMO+ii) = CMO((j-2)*NORB(ISYM)+KCMO+ii)
            end do
          end do

          IndT(ip+1+iii) = m
          Occ(ip+1+KOCC) = q
          ! PAM 2012: NewOrd update
          if (NNWORD > 0) NEWORD(ip+1+KOCC) = nw
          ! PAM 2012: End of update
          do ii=1,NBAS(ISYM)
            CMO((ip)*NORB(ISYM)+KCMO+ii) = TCMO(ii)
          end do

          ip = ip+1
        end if
      end do
    end do
    !write(u6,*) 'sorted:'
    !write(u6,'(10i2)') (IndT(i+iii),i=1,NORB(ISYM))
    !write(u6,'(4ES19.12)') (Occ(i+KOCC),i=1,NORB(ISYM))

    !do ii=1,NBAS(ISYM)
    !  write(u6,*)
    !  write(u6,'(4ES19.12)') (CMO(i+KOCC+(ii-1)*NBAS(ISYM)),i=1,NORB(ISYM))
    !end do

    ! End sort
  end if
  ! Next symmetry
  KCMO = KCMO+NBAS(ISYM)*NORB(ISYM)
  KOCC = KOCC+NORB(ISYM)
  iii = iii+NORB(ISYM)

end do
call mma_deallocate(TCMO)

return

end subroutine VECSORT
