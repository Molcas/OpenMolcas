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

subroutine RDMGRD(JOB,IDISP,LABEL,STYPE,ISYMP,NARRAY,ARRAY)
! Purpose: Read in the derivatives of 1-electron integrals
! of some operator, with respect to some displacement IDISP.
! ISYMP is the symmetry irrep label of the derivatives.

use rassi_aux, only: ipglob
use Cntrl, only: LuMck, MINAME, NJOB
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NBASF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: JOB, IDISP, ISYMP, NARRAY
character(len=8) :: LABEL, STYPE
real(kind=wp) :: ARRAY(NARRAY)
integer(kind=iwp) :: I, IA1, IA2, IAOFF(8), IOPT, IRC, IS, ISCODE, ISUM, ITOFF(8), J, JS, LT, NBI, NBIJ, NBJ, NTEMP
real(kind=wp) :: F
real(kind=wp), allocatable :: TEMP(:)

if ((JOB < 1) .or. (JOB > NJOB)) then
  write(u6,*) ' RASSI/RDMGRD: Invalid JOB parameter.'
  write(u6,*) ' JOB:',JOB
  call ABEND()
end if

if (IPGLOB > 3) then
  write(u6,*) ' RDMGRD called for JOB=',JOB
  write(u6,*) ' perturbed by displacement nr.',IDISP
  write(u6,*) ' MckInt file name:',MINAME(JOB)
  write(u6,*) ' Operator name LABEL=',LABEL
  write(u6,*) ' Symmetry type STYPE=',STYPE
  write(u6,*) ' Irrep label   ISYMP=',ISYMP
  write(u6,*) ' Length NARRAY=',NARRAY
end if

! Open MCKINT file:
IRC = -1
IOPT = 0
call OPNMCK(IRC,IOPT,MINAME(JOB),LUMCK)
if (IRC /= 0) then
  write(u6,*) 'RASSI/RDMGRD: Failed to open '//MINAME(JOB)
  write(u6,*) 'Unit nr LUMCK=',LUMCK
  write(u6,*) 'Option code IOPT=',IOPT
  write(u6,*) 'Return code IRC =',IRC
  call ABEND()
end if

! Addressing integral blocks in the buffer:
ISUM = 0
do IS=1,nIrrep
  JS = MUL(IS,ISYMP)
  if (IS >= JS) then
    ITOFF(IS) = ISUM
    ITOFF(JS) = ISUM
    NBI = NBASF(IS)
    NBJ = NBASF(JS)
    NBIJ = NBI*NBJ
    if (IS == JS) NBIJ = (NBIJ+NBI)/2
    ISUM = ISUM+NBIJ
  end if
end do
NTEMP = ISUM
! Read MCKINT file:
IOPT = 0
ISCODE = 2**(ISYMP-1)
! Get temporary buffer to read data by RDMCK calls
call mma_allocate(TEMP,NTEMP,Label='TEMP')
! Read 1-electron integral derivatives:
IRC = NTEMP
call dRDMCK(IRC,IOPT,LABEL,IDISP,TEMP,ISCODE)
if (IRC /= 0) then
  write(u6,*) 'RDMGRD: RDMGRD failed to read '//MINAME(JOB)
  write(u6,*) '  Displacement IDISP=',IDISP
  write(u6,*) '    Option code IOPT=',IOPT
  write(u6,*) '    Data label LABEL=',LABEL
  write(u6,*) 'Symmetry code ISCODE=',ISCODE
  write(u6,*) '    Return code IRC =',IRC
  call ABEND()
end if

! Addressing integral blocks in ARRAY:
ISUM = 0
do IS=1,nIrrep
  JS = MUL(IS,ISYMP)
  IAOFF(IS) = ISUM
  NBI = NBASF(IS)
  NBJ = NBASF(JS)
  NBIJ = NBI*NBJ
  ISUM = ISUM+NBIJ
end do
if (ISUM > NARRAY) then
  write(u6,*) 'RASSI/RDMGRD: Output ARRAY has insufficient length.'
  write(u6,*) ' Input parameter NARRAY=',NARRAY
  write(u6,*) ' Needed size       ISUM=',ISUM
  call ABEND()
end if
! Move buffer integrals into ARRAY in proper format:
do IS=1,nIrrep
  NBI = NBASF(IS)
  if (NBI <= 0) goto 11
  if (ISYMP == 1) then
    LT = 1+ITOFF(IS)
    IA1 = 1+IAOFF(IS)
    call SQUARE(TEMP(LT),ARRAY(IA1),1,NBI,NBI)
    if (STYPE(1:4) == 'ANTI') then
      do J=1,NBI-1
        do I=J+1,NBI
          ARRAY(IA1-1+I+NBI*(J-1)) = -ARRAY(IA1-1+J+NBI*(I-1))
        end do
      end do
    end if
  else
    JS = MUL(IS,ISYMP)
    if (IS < JS) goto 11
    NBJ = NBASF(JS)
    if (NBJ <= 0) goto 11
    LT = 1+ITOFF(IS)
    IA1 = 1+IAOFF(IS)
    IA2 = 1+IAOFF(JS)
    call DCOPY_(NBI*NBJ,TEMP(LT),1,ARRAY(IA1),1)
    F = One
    if (STYPE(1:4) == 'ANTI') F = -F
    do I=1,NBI
      do J=1,NBJ
        ARRAY(IA2-1+J+NBJ*(I-1)) = F*ARRAY(IA1-1+I+NBI*(J-1))
      end do
    end do
  end if
11 continue
end do
! Get rid of temporary buffer
call mma_deallocate(TEMP)

! Close MCKINT file:
IRC = -1
IOPT = 0
call CLSMCK(IRC,IOPT)
if (IRC /= 0) then
  write(u6,*) 'RASSI/RDMGRD: Failed to close '//MINAME(JOB)
  write(u6,*) 'Unit nr LUMCK=',LUMCK
  write(u6,*) 'Option code IOPT=',IOPT
  write(u6,*) 'Return code IRC =',IRC
  call ABEND()
end if

end subroutine RDMGRD
