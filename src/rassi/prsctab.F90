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

subroutine PRSCTAB(SCTAB,TRANS)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: SCTAB(*)
real(kind=wp), intent(in) :: TRANS(*)
integer(kind=iwp) :: IBLK, IOPEN, ITYPE, KSPCPL, KSPDET, LTRANS, MAXOP, MINOP, MLTPL, MS2, N, NBLK, NCPL, ND, NSIZE, NTRANS
integer(kind=iwp), external :: ngene

write(u6,*)
write(u6,*) '------------------------------------------'
write(u6,*) ' Spin Coupling Table printout'
write(u6,*) '------------------------------------------'
NSIZE = SCTAB(1)
ITYPE = SCTAB(2)
MLTPL = SCTAB(3)
MS2 = SCTAB(4)
MINOP = SCTAB(5)
MAXOP = SCTAB(6)
LTRANS = SCTAB(7)
NTRANS = SCTAB(8)
write(u6,'(1x,A,I16)') ' Table size       :',NSIZE
write(u6,'(1x,A,I16)') ' Table type ID    :',ITYPE
write(u6,'(1x,A,I16)') ' Spin multiplicity:',MLTPL
write(u6,'(1x,A,I16)') ' Spin projection  :',MS2
write(u6,'(1x,A,I16)') ' Open shells; min :',MINOP
write(u6,'(1x,A,I16)') ' Open shells; max :',MAXOP
write(u6,'(1x,A,I16)') ' Transf data; addr:',LTRANS
write(u6,'(1x,A,I16)') ' Transf data; wrds:',NTRANS
! Number of (non-trivial) values of IOPEN:
N = 0
do IOPEN=MINOP,MAXOP
  if (NGENE(IOPEN,MLTPL) > 0) N = N+1
end do
if (N == 0) then
  write(u6,*)
  write(u6,*) ' There is no such spin-coupling scheme.'
  write(u6,*)
else
  write(u6,'(1x,A,I9)') '   Nr of schemes  :',N
  NBLK = MAXOP-MINOP+1
  do IBLK=1,NBLK
    IOPEN = SCTAB(9+(IBLK-1)*6)
    NCPL = SCTAB(10+(IBLK-1)*6)
    if (NCPL /= 0) then
      ND = SCTAB(11+(IBLK-1)*6)
      KSPCPL = SCTAB(12+(IBLK-1)*6)
      KSPDET = SCTAB(13+(IBLK-1)*6)
      LTRANS = SCTAB(14+(IBLK-1)*6)
      write(u6,*) '------------------------------------------'
      write(u6,'(1x,A,I16)') ' Nr of open shells  :',IOPEN
      write(u6,'(1x,A,I16)') ' Nr of proto-CSF    :',NCPL
      write(u6,'(1x,A,I16)') ' Nr of proto-SD     :',ND
      write(u6,'(1x,A,I16)') ' Addr of proto-CSF  :',KSPCPL
      write(u6,'(1x,A,I16)') ' Addr of proto-SD   :',KSPDET
      write(u6,'(1x,A,I16)') ' Addr of transf matr:',LTRANS
      write(u6,*) ' proto-CSF''s:'
      call PRPCSF(IOPEN,NCPL,SCTAB(KSPCPL))
      write(u6,*) ' proto-SD''s:'
      call PRPDET(IOPEN,ND,SCTAB(KSPDET))
      write(u6,*) ' Transformation matrix:'
      call PRPTRA(ND,NCPL,TRANS)
    end if
  end do
end if
write(u6,*) '------------------------------------------'

end subroutine PRSCTAB
