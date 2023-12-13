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

subroutine t3reainput()
! this routine does:
! 1) read INPDAT file, produced by REORG with mmul,nsym,noa,nob,nva,nvb,norb,eps
! 2) read input file for NIT3 to read (parameters transported through cmm common)

! ####################
! Due to the merging of CC input files to one, to avoid conflicts
! Denominators in CCT3 has become T3Denominators
! and Shift has become T3Shift
!                         (JR) Lund 2003
!   ! title   - jobtitle
!   1-ntit rows with 72 characters
!   no default
!   ! ntit    - number of rowws in jobtitle
!   ntit is limited to 10 CGG From now 1
!   default (1)
!   ! typt3   - type of T3 cpntribution
!   0 - CCSD
!   1 - CCSD+T(CCSD)
!   2 - CCSD(T) Ragh
!   3 - CCSD(T) Bart
!   default (3)
!   ! typden  - type of denominator (division of fok)
!   0 - diagonal elements
!   1 - average of faa and fbb
!   2 - orbital energies
!   default (0)
!   ! keysa   - Spin adaptation key
!   0 - no adaptation
!   1 - T2 DDVV adaptation
!   2 - T2 DDVV + T1 DV adaptation
!   3 - full T1 and T2 adaptation (only for doublets)
!   4 - full T2 adaptation without SDVS (only for doublets)
!   (default=0)
!   ! filerst - name for CCSD results containing file
!   (default=RSTART)
!   ! mchntyp - type of machine in matrix multiplication
!   1 - C=A*B is faster or comparable with C=AT*B
!   2 - C=AT*B is faster
!   (default=1)
!   ! slim    - limitation for usieng C=AT*B
!   no default (suitable=Two)
!   ! shifhto - shift for occupied
!   (default=Zero)
!   ! shifhtv - shift for virtuals
!   (default=Zero)
!   ! maxspace - maximal allowed work space
!   (default=0 - unlimited)
!   ! fullprint - level of printing control key
!   (default=0)
!   ! noop - No Operation key
!   (default=no)
!   & iokey - I/O control key
!     1 - Fortran I/O system
!     2 - MOLCAS DA IO system
!   (default=2)
!   & mhkey - Matrix handling control key
!     1 - ESSL routines
!     2 - Fortran I/O system
!   (default=1)
!   .....   - can be added

use CCT3_global, only: dimm, eps, filerst, fullprint, ijsegkey, imax, imin, iokey, ispin, jmax, jmin, keysa, lsym, maxorb, &
                       maxspace, mchntyp, mhkey, mmul, noa, nob, noop, norb, nshf, nsym, nva, nvb, slim, shifto, shiftv, symimax, &
                       symimin, symjmax, symjmin, typden, typt3
use Constants, only: Zero, One
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: Lu, LuSpool, nactel, nhelp
character(len=80) :: LINE
character(len=72) :: title
integer(kind=iwp), external :: isFreeUnit

#include "macros.fh"

!1 read INPDAT

Lu = isFreeUnit(1)
!open(unit=1,file='INPDAT',form='unformatted')
call molcas_binaryopen_vanilla(Lu,'INPDAT')
read(Lu) nactel,ispin,nsym,lsym,mmul,noa,nob,nva,nvb,norb,eps
unused_var(nactel)
close(Lu)

!2 def dimm

dimm(1,1:nsym) = noa(1:nsym)
dimm(2,1:nsym) = nob(1:nsym)
dimm(3,1:nsym) = nva(1:nsym)
dimm(4,1:nsym) = nvb(1:nsym)
dimm(5,1:nsym) = norb(1:nsym)

!3 define nshf

do nhelp=1,maxorb
  nshf(nhelp) = (nhelp-1)*(nhelp-2)/2
end do

!4 define defaults

typt3 = 3
typden = 0
keysa = 0
filerst = 'RSTART'
mchntyp = 1
slim = One
shifto = Zero
shiftv = Zero
maxspace = 0
!GG fullprint = 0
noop = 0
iokey = 1
mhkey = 1
ijsegkey = 0
symimin = 1
symjmin = 1
symimax = nsym
symjmax = nsym
imin = 0
jmin = 0
imax = 0
jmax = 0

!5 read input file

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
TITLE = ' '
do
  read(LuSpool,'(A80)') LINE
  call UPCASE(LINE)
  if (index(LINE,'&CCSDT') /= 0) exit
end do
do
  read(LuSpool,'(A80)') LINE
  if (LINE(1:1) == '*') cycle
  call UPCASE(LINE)

  if (LINE(1:4) == 'TITL') then
    read(LuSpool,'(A72)') TITLE
    unused_var(title)
  else if (LINE(1:4) == 'TRIP') then
    read(LuSpool,*) typt3
  else if (LINE(1:4) == 'T3DE') then
    read(LuSpool,*) typden
  else if (LINE(1:4) == 'ADAP') then
    read(LuSpool,*) keysa
    if ((keysa > 4) .or. (keysa < 0)) then
      keysa = 0
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, keysa was changed to 0'
      end if
    end if
    if ((keysa /= 0) .and. (typden == 0)) then
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, typden is incompatible with SA'
      end if
    end if
  else if (LINE(1:4) == 'LOAD') then
    read(LuSpool,*) filerst
  else if (LINE(1:4) == 'MACH') then
    read(LuSpool,*) mchntyp,slim
    if ((mchntyp < 1) .or. (mchntyp > 2)) then
      mchntyp = 1
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, mchntyp was changed to 1'
      end if
    end if
  else if (LINE(1:4) == 'T3SH') then
    read(LuSpool,*) shifto,shiftv
  else if (LINE(1:4) == 'PRIN') then
    read(LuSpool,*) fullprint
  else if (LINE(1:4) == 'NOOP') then
    noop = 1
  else if (LINE(1:4) == 'IOKE') then
    read(LuSpool,*) iokey
    if ((iokey < 0) .or. (iokey > 2)) then
      iokey = 2
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, I/O key out of range'
        write(u6,*) ' parameter iokey changed to 2'
      end if
    end if
  else if (LINE(1:4) == 'MHKE') then
    read(LuSpool,*) mhkey
    if ((mhkey < 0) .or. (mhkey > 2)) then
      mhkey = 1
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, Matrix handling key is out of range'
        write(u6,*) ' parameter iokey changed to 1'
      end if
    end if
  else if (LINE(1:4) == 'IJSE') then
    ijsegkey = 1
    read(LuSpool,*) symimin,imin,symjmin,jmin,symimax,imax,symjmax,jmax
  else if (LINE(1:4) == 'END ') then
    exit
  end if
end do

call Close_LuSpool(LuSpool)
return

end subroutine t3reainput
