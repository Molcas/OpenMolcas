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

subroutine RDMCCI(JOB,IDISP,LABEL,ISYMP,NARRAY,ARRAY)
! Purpose: Read in the derivatives of CI array derivatives
! from MCKINT file, with respect to some displacement IDISP.
! ISYMP is the symmetry irrep label of the derivatives.

use rassi_aux, only: ipglob
use stdalloc, only: mma_allocate, mma_deallocate
use Cntrl, only: NJOB, NCONF1, MINAME
use cntrl, only: LuMck
use Definitions, only: u6

implicit none
integer JOB, IDISP, ISYMP, nArray
character(len=8) LABEL
real*8 ARRAY(NARRAY)
real*8, allocatable :: TEMP(:)
integer IRC, IOPT, NTEMP, ISCODE

if ((JOB < 1) .or. (JOB > NJOB)) then
  write(u6,*) ' RDMCI: Invalid JOB parameter.'
  write(u6,*) ' JOB:',JOB
  call ABEND()
end if

if (IPGLOB >= 3) then
  write(u6,*) ' RDMCCI called for JOB=',JOB
  write(u6,*) ' perturbed by displacement nr.',IDISP
  write(u6,*) ' MckInt file name:',MINAME(JOB)
  write(u6,*) ' Irrep label   ISYMP=',ISYMP
  write(u6,*) ' Length NARRAY=',NARRAY
end if

! Open MCKINT file:
IRC = -1
IOPT = 0
call OPNMCK(IRC,IOPT,MINAME(JOB),LUMCK)
if (IRC /= 0) then
  write(u6,*) 'RDMCCI: Failed to open '//MINAME(JOB)
  write(u6,*) 'Unit nr LUMCK=',LUMCK
  write(u6,*) 'Option code IOPT=',IOPT
  write(u6,*) 'Return code IRC =',IRC
  call ABEND()
end if

! Read MCKINT file:
NTEMP = NCONF1
IOPT = 0
ISCODE = 2**(ISYMP-1)
! Get temporary buffer to read data by RDMCK calls
call mma_allocate(TEMP,NTEMP,Label='TEMP')
! Read 1-electron integral derivatives:
IRC = NTEMP
call dRDMCK(IRC,IOPT,LABEL,IDISP,TEMP,ISCODE)
if (IRC /= 0) then
  write(u6,*) 'RDMCCI: RDMCCI failed to read '//MINAME(JOB)
  write(u6,*) '  Displacement IDISP=',IDISP
  write(u6,*) '    Option code IOPT=',IOPT
  write(u6,*) 'Symmetry code ISCODE=',ISCODE
  write(u6,*) '    Return code IRC =',IRC
  call ABEND()
end if

if (NTEMP > NARRAY) then
  write(u6,*) 'RDMCCI: Output ARRAY has insufficient length.'
  write(u6,*) ' Input parameter NARRAY=',NARRAY
  write(u6,*) ' Needed size       NTEMP=',NTEMP
  call ABEND()
end if
! Move buffer integrals into ARRAY in proper format:
call DCOPY_(NTEMP,TEMP,1,ARRAY,1)
! Get rid of temporary buffer
call mma_deallocate(TEMP)

! Close MCKINT file:
IRC = -1
IOPT = 0
call CLSMCK(IRC,IOPT)
if (IRC /= 0) then
  write(u6,*) 'RDMCCI: Failed to close '//MINAME(JOB)
  write(u6,*) 'Unit nr LUMCK=',LUMCK
  write(u6,*) 'Option code IOPT=',IOPT
  write(u6,*) 'Return code IRC =',IRC
  call ABEND()
end if

end subroutine RDMCCI
