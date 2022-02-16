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

subroutine DENSCT(AREF)

use mrci_global, only: CSPCK, DMO, ICPF, INDX, INTSY, JREFX, ITRANS, LUEIG, LUREST, NBAST, NBTRI, NCONF, NRROOT, NVMAX, NVSQ, TDMO
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: AREF(*)
integer(kind=iwp) :: I, IDDMO, IDREST, J
integer(kind=iwp), allocatable :: IDC(:)
real(kind=wp), allocatable :: ASCR2(:), BSCR2(:), CI(:), FSCR2(:), SGM(:)

IDREST = 0
IDDMO = 0
call mma_allocate(CI,NCONF,label='CI')
call mma_allocate(SGM,NCONF,label='SGM')
call mma_allocate(ASCR2,NVMAX**2,label='ASCR2')
call mma_allocate(BSCR2,NVMAX**2,label='BSCR2')
call mma_allocate(FSCR2,NVSQ,label='FSCR2')
call mma_allocate(IDC,NRROOT,label='IDC')
do I=1,NRROOT
  IDC(I) = IDREST
  call dDAFILE(LUREST,2,CI,NCONF,IDREST)
  call FZERO(DMO,NBTRI)
  !PAM04 if (ICPF /= 0) call DCORR(JREFX,HWork(LAREF),CSPCK,INTSY,INDX,DMO)
  if (ICPF /= 0) call DCORR(JREFX,AREF,CSPCK,INTSY,INDX,DMO)
  !PAM04 call FIJD(INTSY,INDX,CI,DMO,JREFX,HWork(LAREF))
  call FIJD(INTSY,INDX,CI,DMO,JREFX,AREF)
  !PAM04 call AID(INTSY,INDX,CI,DMO,
  call AID(INTSY,INDX,CI,DMO,ASCR2,BSCR2,FSCR2)
  !PAM04 call ABD(CSPCK,INTSY,INDX,CI,DMO,JREFX)
  call ABD(CSPCK,INTSY,INDX,CI,DMO,ASCR2,BSCR2,FSCR2,JREFX)
  !PAM04 call dDAFILE(LUEIG,1,DMO,NBTRI,IDDMO)
  call dDAFILE(LUEIG,1,DMO,NBTRI,IDDMO)
end do
if (ITRANS /= 0) then
  do I=2,NRROOT
    IDREST = IDC(I)
    call dDAFILE(LUREST,2,CI,NCONF,IDREST)
    do J=1,I-1
      IDREST = IDC(J)
      call dDAFILE(LUREST,2,SGM,NCONF,IDREST)
      !PAM04 call FZERO(TDMO,NBAST**2)
      call FZERO(TDMO,NBAST**2)
      !PAM04 call FIJTD(INTSY,INDX,CI,SGM,TDMO)
      call FIJTD(INTSY,INDX,CI,SGM,TDMO)
      !PAM04 call AITD(INTSY,INDX,CI,TDMO,ASCR2,BSCR2,
      call AITD(INTSY,INDX,CI,SGM,TDMO,ASCR2,BSCR2,FSCR2)
      !PAM04 call ABTD(CSPCK,INTSY,INDX,TDMO,ASCR2,BSCR2,
      call ABTD(CSPCK,INTSY,INDX,CI,SGM,TDMO,ASCR2,BSCR2,FSCR2)
      !PAM04 call dDAFILE(LUEIG,1,TDMO,NBAST**2,IDDMO)
      call dDAFILE(LUEIG,1,TDMO,NBAST**2,IDDMO)
    end do
  end do
end if
call mma_deallocate(CI)
call mma_deallocate(SGM)
call mma_deallocate(ASCR2)
call mma_deallocate(BSCR2)
call mma_deallocate(FSCR2)
call mma_deallocate(IDC)

return

end subroutine DENSCT
