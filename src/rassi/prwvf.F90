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

subroutine PRWVF(IORBTAB,ISSTAB,IFSBTAB,PRTHR,CI)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: u6

implicit none
integer IORBTAB(*), ISSTAB(*), IFSBTAB(*)
real*8 CI(*), PRTHR
integer IBLKDET, NPRINTED
integer NASPRT, KSPART, KSSTTB, KSBSMRS
integer NFSB, KSTARR, NSSTP, ISUM, ISST, NSBS
integer IFSB, KPOS, ISPART, NBLKDET
integer ISSTARR(50), NDIMARR(50)
integer ICREL, ICI, INDX, ISORB, N, I, ISBS, IMORS
character(len=80) DETTXT
integer, allocatable :: SBSET(:)

! The orbital table:
NASPRT = IORBTAB(9)
KSPART = IORBTAB(10)
! The substring table:
NSSTP = ISSTAB(7)
KSSTTB = 15
KSBSMRS = ISSTAB(11)
! The FS blocks of the SGM wave function:
NFSB = IFSBTAB(3)
KSTARR = 8
! Make an array with nr of earlier substrings for each
! substring type:
call mma_allocate(SBSET,NSSTP,Label='SBSET')
ISUM = 0
do ISST=1,NSSTP
  SBSET(ISST) = ISUM
  NSBS = ISSTAB(KSSTTB+5*(ISST-1))
  ISUM = ISUM+NSBS
end do
! Loop over FS blocks of the SGM wave function
NPRINTED = 0
do IFSB=1,NFSB
  KPOS = KSTARR+(NASPRT+2)*(IFSB-1)
  do ISPART=1,NASPRT
    ISSTARR(ISPART) = IFSBTAB(KPOS-1+ISPART)
  end do
  NBLKDET = IFSBTAB(KPOS+NASPRT)
  IBLKDET = IFSBTAB(KPOS+NASPRT+1)
  ! Dimension of each substring type:
  do ISPART=1,NASPRT
    ISST = ISSTARR(ISPART)
    NSBS = ISSTAB(KSSTTB+5*(ISST-1))
    NDIMARR(ISPART) = NSBS
  end do
  do ICREL=1,NBLKDET
    ICI = IBLKDET-1+ICREL
    if (abs(CI(ICI)) >= PRTHR) then
      ! Get occupation array in the form of string DETTXT:
      INDX = ICREL-1
      ISORB = 0
      do ISPART=1,NASPRT
        N = NDIMARR(ISPART)
        I = mod(INDX,N)
        INDX = INDX-N*I
        ISST = ISSTARR(ISPART)
        ISBS = SBSET(ISST)+I+1
        IMORS = ISSTAB(KSBSMRS+2*(ISBS-1))
        N = IORBTAB(KSPART-1+ISPART)
        call MORSWRITE(IMORS,DETTXT(ISORB+1:ISORB+N))
        ISORB = ISORB+N
      end do
      write(u6,'(1x,a,5x,f16.8)') DETTXT(1:ISORB),CI(ICI)
      NPRINTED = NPRINTED+1
    end if
  end do
  ! End of loop over FS blocks
end do
if (NPRINTED == 0) write(u6,*) ' (PRWVF: Nothing worth printing)'
call mma_deallocate(SBSET)

end subroutine PRWVF
