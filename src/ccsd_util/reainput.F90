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

subroutine reainput()
! this routine does:
! 1) read INPDAT file, produced by REORG with mul,nsym,noa,nob,nva,nvb,norb,eps
! 2) read input file for CCSD to read (parameters transported through cmm common)
! title   - jobtitle
!   1-ntit rows with 72 characters
!   no default
! ntit    - number of rows in jobtitle
!   ntit is limited to 10 !!! NOT ANY MORE !!! ntit is limited to 1
!   default (1)
! maxiter - maximum number of iterations
!   default (10)
! typt3   - type of T3 cpntribution
!   0 - CCSD
!   1 - CCSD+T(CCSD)
!   2 - CCSD(T)
!   default (0)
! typden  - type of denominator (division of fok)
!   0 - diagonal elements
!   1 - average of faa and fbb
!   2 - orbital energies
!   default (0)
! firstext- first iteration wnere extrapolation is used
!   (default-no) less than cycext is not possible
! cycext  - cycle of extrapolation
!   (default-no) limited to 2-4
! ccnonv  - energy convergence criterion
!   (default=1.0e-6)
! keysa   - Spin adaptation key
!   0 - no adaptation
!   1 - T2 DDVV adaptation
!   2 - T2 DDVV + T1 DV adaptation
!   3 - full T1 and T2 adaptation (only for doublets)
!   4 - full T2 adaptation without SDVS (only for doublets)
!   (default=0)
! keyrst  - restart key
!   0 - no saving restart informations (amplitudes)
!   1 - save amplitudes
!   2 - start from previous informations
!   (default=1)
! filerst - name for restart informations file
!   (default=RSTART)
! mchntyp - type of machine in matrix multiplication
!   1 - C=A*B is faster or comparable with C=AT*B
!   2 - C=AT*B is faster
! (default=1)
! slim    - limitation for usieng C=AT*B
!   no default (suitable=2.0)
! shifhto - shift for occupied
!   (default=0.0)
! shifhtv - shift for virtuals
!   (default=0.0)
! maxspace - maximal allowed work space
!   (default=0 - unlimited)
! fullprint - level of printing contrlo key
!   (default=0)
! noop  - no operation key
!   (default=no)
! iokey - I/O control key
!   1 - Fortran I/O system
!   2 - MOLCAS DA IO system
!   (default=2)
! mhkey - Matrix handling control key
!   1 - ESSL routines
!   2 - Fortran I/O system
! noccsd - key to supress CCSD run
!   (default=no)c     (default=1)
! .....  - can be added

use ccsd_global, only: ccconv, cycext, dimm, eps, Escf, filerst, firstext, fullprint, iokey, ispin, keyrst, keysa, lsym, maxiter, &
                       maxorb, maxspace, mchntyp, mhkey, mmul, nactel, noa, nob, noccsd, noop, norb, nshf, nsym, ntit, nva, nvb, &
                       slim, shifto, shiftv, title, typden, typt3, yesext
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: f_iostat, f_recl, LuSpool, nhelp
character(len=80) :: LINE
logical(kind=iwp) :: is_error

!1 read INPDAT

call molcas_open_ext2(1,'INPDAT','sequential','unformatted',f_iostat,.false.,f_recl,'unknown',is_error)
!open(unit=1,file='INPDAT',form='unformatted')
read(1) nactel,ispin,nsym,lsym,mmul,noa,nob,nva,nvb,norb,eps,Escf
close(1)

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

maxiter = 30
typt3 = 0
ntit = 1
typden = 0
yesext = 0
firstext = 0
cycext = 0
ccconv = 1.0e-7_wp
keysa = 0
keyrst = 1
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
noccsd = 0

!5 read input file

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
do
  read(LuSpool,'(A80)') LINE
  call UPCASE(LINE)
  if (index(LINE,'&CCSDT') /= 0) exit
end do
TITLE = ' '
do
  read(LuSpool,'(A80)') LINE
  if (LINE(1:1) == '*') cycle
  call UPCASE(LINE)

  if (LINE(1:4) == 'TITL') then
    read(LuSpool,'(A72)') TITLE
  else if (LINE(1:4) == 'ITER') then
    read(LuSpool,*) maxiter
  else if (LINE(1:4) == 'DENO') then
    read(LuSpool,*) typden
    if ((typden < 0) .or. (typden > 2)) then
      typden = 2
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, Invalid type of denominators'
        write(u6,*) ' parameter typden changed to 2'
      end if
    end if
  else if (LINE(1:4) == 'EXTR') then
    yesext = 1
    read(LuSpool,*) firstext,cycext
    if ((cycext < 2) .or. (cycext > 4)) then
      cycext = 4
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, Size of DIIS procedure out of range'
        write(u6,*) ' parameter cycext changed to 4'
      end if
    end if
    if (firstext < cycext) then
      firstext = cycext
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, First DIIS iteration is smaller then DIIS size'
        write(u6,*) ' parameter firstext was changed to:',firstext
      end if
    end if
  else if (LINE(1:4) == 'ACCU') then
    read(LuSpool,*) ccconv
  else if (LINE(1:4) == 'ADAP') then
    read(LuSpool,*) keysa
    if ((keysa > 4) .or. (keysa < 0)) then
      keysa = 0
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, Adaptation key out of range'
        write(u6,*) ' parameter keysa changed to 0'
      end if
    end if
    if ((keysa /= 0) .and. (typden == 0)) then
      typden = 2
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, typden is incompatible with SA'
        write(u6,*) ' type of denominators changed to 2 - Orb. energies'
      end if
    end if
  else if (LINE(1:4) == 'REST') then
    read(LuSpool,*) keyrst
    if ((keyrst < 0) .or. (keyrst > 2)) then
      keyrst = 1
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, Restart key out of range'
        write(u6,*) ' parameter keyrst changed to 1'
      end if
    end if
    read(LuSpool,*) filerst
  else if (LINE(1:4) == 'MACH') then
    read(LuSpool,*) mchntyp,slim
    if ((mchntyp < 1) .or. (mchntyp > 2)) then
      mchntyp = 1
      if (fullprint >= 0) then
        write(u6,*) ' Warning!!!, Machinetype out of range'
        write(u6,*) ' parameter mchtyp changed to 1'
      end if
    end if
  else if (LINE(1:4) == 'SHIF') then
    read(LuSpool,*) shifto,shiftv
  else if (LINE(1:4) == 'PRIN') then
    read(LuSpool,*) fullprint
    if ((fullprint < 0) .or. (fullprint > 3)) then
      fullprint = 0
      write(u6,*) ' Warning!!!, Printing key out of range'
      write(u6,*) ' parameter fullprint changed to 0'
    end if
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
        write(u6,*) ' Warning!!!, Matrix handling key out of range'
        write(u6,*) ' parameter mhkey changed to 1'
      end if
    end if
  else if (LINE(1:4) == 'NOSD') then
    noccsd = 1
  else if (LINE(1:4) == 'END ') then
    exit
  end if
end do

call Close_LuSpool(LuSpool)

return

end subroutine reainput
