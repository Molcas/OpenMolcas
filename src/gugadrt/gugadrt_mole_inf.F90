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

subroutine gugadrt_mole_inf()

use gugadrt_global, only: iprint, iref_occ, logic_mr, logic_mrelcas, lsm_inn, ludrt, max_ref, n_electron, n_ref, ng_sm, nlsm_all, &
                          nlsm_bas, nlsm_ext, nlsmddel, nlsmedel, norb_act, norb_all, norb_dbl, norb_dz, norb_ext, norb_frz, &
                          norb_inn, ns_sm, nstart_act, spin
use Symmetry_Info, only: Mul
use Constants, only: Zero, Two
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: err, i, icmd, idisk, idum(1), im, iorb, ispin, itmpstr(72), j, jcmd, ln1, lsm, lsmtmp(8), ms_ref, mult, &
                     nact_sm, nactel, nde, ndisk, ne_act, neact, ngsm, nlsm_act(8), nlsm_dbl(8), nlsm_frz(8), nlsm_inn(8), &
                     noidx(8), norb1, norb2, norb_all_tmp, ntit
logical(kind=iwp) :: log_debug, skip
character(len=4) :: command
character(len=72) :: line
character(len=132) :: modline
! copy from molcas, bsuo, jun. 30, 2009
integer(kind=iwp), parameter :: ncmd = 18, mxtit = 10
character(len=*), parameter :: cmd(ncmd) = ['TITL','ELEC','SPIN','SYMM','ACTI','PRIN','REFE','FIRS','INAC','CIAL', &
                                            'VALE','INTE','NOCO','ONEO','EXTR','NONI','NACT','END ']

n_electron = 0
nactel = 0
spin = Zero
ns_sm = 1
call rdnlst(5,'GUGADRT')
ntit = 0
skip = .false.
jcmd = 0
do
  if (.not. skip) then
    read(5,'(a)',iostat=err) line
    if (err < 0) call error(1)
    command = line(1:8)
    call upcase(command)
    if (command(1:1) == '*') cycle
    if (command == ' ') cycle
    jcmd = 0
    do icmd=1,ncmd
      if (command == cmd(icmd)) jcmd = icmd
    end do
  end if
  skip = .false.
  select case (jcmd)
    case (1)
      !---  process title    command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        command = line(1:8)
        call upcase(command)
        if (command(1:1) == '*') cycle
        jcmd = 0
        do icmd=1,ncmd
          if (command == cmd(icmd)) jcmd = icmd
        end do
        if (jcmd /= 0) then
          skip = .true.
          exit
        end if
        ntit = ntit+1
      end do

    case (2)
      !---  process electron command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=err) n_electron
      if (err > 0) call error(2)

    case (3)
      !---  process spin     command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=err) ispin
      if (err > 0) call error(2)
      spin = (ispin-1)/Two

    case (4)
      !---  process symmetry command ----------------------------------------*
      !write (u6,*)'input_guga: keyword symmetry is obsolete and ignored!
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=err) ns_sm
      if (err > 0) call error(2)

    case (5)
      !---  process active   command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      modline = line//' 0 0 0 0 0 0 0 0'
      read(modline,*,iostat=err) (nlsm_act(i),i=1,8)
      if (err > 0) call error(2)

    case (6)
      !---  process print    command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=err) iprint
      if (err > 0) call error(2)

    case (7)
      !---  process referenc command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=err) n_ref,ln1
      if (err > 0) call error(2)
      if (n_ref > max_ref) then
        write(u6,*) ' Warnning! Program could not deal with so many reference states!'
        write(u6,*) ' Maximum number of reference states is',max_ref
        write(u6,*) ' Set number of reference states into ',max_ref
        n_ref = max_ref
      end if
      if (ln1 == 0) then
        logic_mr = .true.
      else
        do i=1,n_ref
          read(5,'(80i1)',iostat=err) iref_occ(1:ln1,i)
          if (err < 0) call error(1)
          if (err > 0) call error(2)
        end do
        logic_mr = .true.
      end if

    case (8)
      !---  process first    command ----------------------------------------*

    case (9)
      !---  process inactive command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      modline = line//' 0 0 0 0 0 0 0 0'
      read(modline,*,iostat=err) (nlsm_dbl(i),i=1,8)
      if (err > 0) call error(2)

    case (10)
      !---  process ciall    command ----------------------------------------*
      !do
      !  read(5,'(a)',iostat=err) line
      !  if (err < 0) call error(1)
      !  if (line(1:1) /= '*') exit
      !end do
      !read(line,*,iostat=err) ns_sm
      !if (err > 0) call error(2)
      logic_mrelcas = .true.

    case (11)
      !---  process valence  command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      modline = line//' 0 0 0 0 0 0 0 0'
      read(modline,*,iostat=err) (nlsm_ext(i),i=1,8)
      if (err > 0) call error(2)

    case (12)
      !---  process interact command ----------------------------------------*

    case (13)
      !---  process nocorr   command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      modline = line//' 0 0 0 0 0 0 0 0'
      read(modline,*,iostat=err)  !(ncor(i),i=1,8)
      if (err > 0) call error(2)
      !bsuo, jun. 30, 2009 - neglect them

    case (14)
      !---  process oneocc   command ----------------------------------------*
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      modline = line//' 0 0 0 0 0 0 0 0'
      read(modline,*,iostat=err) ! (ione(i),i=1,8)
      if (err > 0) call error(2)
      !bsuo, jun. 30, 2009 - neglect them

    case (15)
      !---  process extract  command ----------------------------------------*
      write(u6,*) 'input: extract option is redundant and is ignored!'

    case (16)
      !---  process non-interact command -------------------------------------
      write(u6,*) 'input: non-interact option is redundant and is ignored!'

    case (17)
      !---  process nactel       command -------------------------------------
      do
        read(5,'(a)',iostat=err) line
        if (err < 0) call error(1)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=err) nactel
      if (err > 0) call error(2)

    case (18)
      !---  the end of the input is reached, print the title ----------------*
      exit

    case default
      write(u6,*) 'input: illegal keyword'
      write(u6,'(a,a)') 'command=',command
      call abend()

  end select

end do

! frozen orbital(dbl, ext) have been delete in mo trans step, so we negl here.
nlsm_frz(1:ng_sm) = 0

norb_frz = 0
norb_act = 0
norb_dz = 0
norb_all = 0
do i=1,ng_sm
  norb_frz = norb_frz+nlsm_frz(i)
  norb_dbl = norb_dbl+nlsm_dbl(i)
  norb_act = norb_act+nlsm_act(i)
  nlsm_inn(i) = nlsm_frz(i)+nlsm_dbl(i)+nlsm_act(i)
  nlsm_ext(i) = nlsm_all(i)-nlsm_inn(i)
  norb_all = norb_all+nlsm_all(i)
end do
norb_inn = norb_frz+norb_dbl+norb_act
norb_ext = norb_all-norb_inn
norb_dz = norb_dbl+norb_frz

nstart_act = norb_dz+1

iorb = 0
do i=1,ng_sm
  do j=1,nlsm_frz(i)
    iorb = iorb+1
    lsm_inn(iorb) = i
  end do
end do
do i=1,ng_sm
  do j=1,nlsm_dbl(i)
    iorb = iorb+1
    lsm_inn(iorb) = i
  end do
end do
do i=1,ng_sm
  do j=1,nlsm_act(i)
    iorb = iorb+1
    lsm_inn(iorb) = i
  end do
end do

nact_sm = 1
do im=1,norb_inn
  if (lsm_inn(im) > nact_sm) nact_sm = lsm_inn(im)
end do

norb_all_tmp = 0
do ngsm=1,ng_sm
  norb_all_tmp = norb_all_tmp+nlsm_all(ngsm)
end do
if (norb_all_tmp /= norb_all) then
  write(u6,*) '  input num.of orbital err! check again!'
  call abend()
end if

lsm = 0
do im=ng_sm,1,-1
  if (nlsm_ext(im) > 0) then
    lsm = im
    exit
  end if
end do
lsm_inn(norb_inn+1) = lsm

lsmtmp(1:8) = 0
if (logic_mr) then
  do i=1,n_ref
    itmpstr(1:norb_act) = iref_occ(1:norb_act,i)
    do j=1,norb_dz
      iref_occ(j,i) = 2
    end do
    iref_occ(norb_dz+1:norb_dz+norb_act,i) = itmpstr(1:norb_act)
  end do
end if

do i=1,ng_sm !norb_inn
  lsmtmp(i) = 0
end do
do i=1,norb_inn
  lsmtmp(lsm_inn(i)) = lsmtmp(lsm_inn(i))+1
end do

!============ chck input data block ====================
if (norb_frz > norb_dz) then
  write(u6,*) ' check input data: norb_frz  error'
  call abend()
end if
! check number of electrons
if (nactel /= 0) then
  if (n_electron /= 0) then
    ne_act = n_electron-2*norb_dz
    if (nactel /= ne_act) then
      write(u6,*) 'Input error, Error in checking ELECtron and NACTel!'
      call abend()
    end if
  else
    ne_act = nactel
    n_electron = 2*norb_dz+nactel
  end if
else
  if (n_electron /= 0) then
    ne_act = n_electron-2*norb_dz
    nactel = ne_act
  else
    ne_act = 0
    write(u6,*) 'Input error, you must input ELECtron or NACTel!'
    call abend()
  end if
end if
nde = 2*norb_dz
if (nde > n_electron) then
  write(u6,'(1x,42a)') 'check input date: number of electrons error'
  write(u6,'(1x,a20,1x,i4)') 'number of electrons ',n_electron
  write(u6,'(1x,a36,1x,i4)') 'number of doubly occupied electrons ',nde
  call abend()
end if

norb1 = norb_inn+norb_ext
norb2 = 0
do i=1,ng_sm
  norb2 = norb2+nlsm_all(i)
end do
if (norb1 /= norb2) then
  write(u6,*) ' check input data: orb_number error'
  write(u6,*) norb1,norb2
  call abend()
end if

if (logic_mr) then
  write(u6,*) ' Refrence states'
  !ne_act = n_electron-2*norb_dz
  do i=1,n_ref
    write(u6,'(80i1)') iref_occ(1:norb_inn,i)
    ms_ref = 1
    neact = 0
    do j=norb_dz+1,norb_inn
      if (iref_occ(j,i) == 1) ms_ref = Mul(ms_ref,lsm_inn(j))
      if (iref_occ(j,i) == 1) neact = neact+1
      if (iref_occ(j,i) == 2) neact = neact+2
    end do
    if (ms_ref /= ns_sm .or. neact /= ne_act) then
      write(u6,*) '  input ref_conf  err! check again! ',i
      call abend()
    end if
  end do
end if
!============ block end =============================
!****************************************************
log_debug = .false.
if (log_debug) then
  write(u6,1001)
  write(u6,1002) norb_all,norb_frz,norb_dz,norb_inn,norb_ext
  write(u6,1002) nlsm_all(1:ng_sm)
  write(u6,1002) nlsm_frz(1:ng_sm)
  write(u6,1002) nlsm_dbl(1:ng_sm)
  write(u6,1002) nlsm_inn(1:ng_sm)
  write(u6,*) 'lsm inn'
  write(u6,1002) lsm_inn(1:norb_inn)
end if
!*****************************************************

! merge into molcas
! write date into cidrt for ci calculation
noidx = 0
idisk = 0
call idafile(ludrt,1,noidx,2,idisk)
! group symmetry
idum(1) = ng_sm
call idafile(ludrt,1,idum,1,idisk)
! state symmetry
idum(1) = ns_sm
call idafile(ludrt,1,idum,1,idisk)
! number of roots to be cal
idum(1) = 1
call idafile(ludrt,1,idum,1,idisk)
! number of corelation electrons
idum(1) = n_electron
call idafile(ludrt,1,idum,1,idisk)
! number of active electrons
idum(1) = nactel
call idafile(ludrt,1,idum,1,idisk)
! spin symmetry of the state, 2s+1
idum(1) = ispin
call idafile(ludrt,1,idum,1,idisk)
! dbl orb
call idafile(ludrt,1,nlsm_dbl,8,idisk)
! act orb
call idafile(ludrt,1,nlsm_act,8,idisk)
! all correlated orb
call idafile(ludrt,1,nlsm_all,8,idisk)
! num. basis
call idafile(ludrt,1,nlsm_bas,8,idisk)
! method to choose ref. state, 4 for cas, 2 for rst
if (logic_mr) then
  idum(1) = 2
  call idafile(ludrt,1,idum,1,idisk)
  idum(1) = n_ref
  call idafile(ludrt,1,idum,1,idisk)
  ! reference states
  !iref_occ(norb_dz+1:norb_dz+norb_act,i)=itmpstr(1:norb_act)
  do i=1,n_ref
    call idafile(ludrt,1,iref_occ(1,i),norb_inn,idisk)
  end do
end if
if (logic_mrelcas) then
  idum(1) = 4
  call idafile(ludrt,1,idum,1,idisk)
end if
! number of ref. states, if logic_mr = .true.
noidx(2) = idisk
idisk = 0
call idafile(ludrt,1,noidx,2,idisk)
noidx = 0
do i=1,ng_sm
  noidx(i) = i
end do
mult = nint(2*spin)+1
write(u6,*) '-----------------------------------------------'
write(u6,*) '    ci orbitals information'
ndisk = 0
do i=1,8
  ndisk = ndisk+nlsm_bas(i)
end do
write(u6,*) '    num. of basis:          ',ndisk
write(u6,*) '    num. of orbitals:       ',norb_all
write(u6,*) '    num. of active-orbitals:',norb_act
write(u6,*) '    num. of electrons:      ',n_electron
write(u6,*) '    multiplicity:           ',mult
write(u6,*) '    symmetry:               ',ns_sm
write(u6,*)
write(u6,*) '    oribtials per-symmtry'
write(u6,1003) noidx(1:ng_sm)
write(u6,1004) nlsmddel(1:ng_sm)
write(u6,1005) nlsm_dbl(1:ng_sm)
write(u6,1006) nlsm_act(1:ng_sm)
write(u6,1007) nlsm_ext(1:ng_sm)
write(u6,1008) nlsmedel(1:ng_sm)
write(u6,*) '-----------------------------------------------'

! following information will be drts
! notice we do not close file ludrt here, it will be closed after drt
! information is written into this file
!---------------------------------------------------------------------

return

1001 format(1x,'norb all group sm')
1002 format(8(1x,i3))
1003 format(16x,8(i4))
1004 format(8x,'frozen  ',8(i4))
1005 format(8x,'double  ',8(i4))
1006 format(8x,'active  ',8(i4))
1007 format(8x,'virtual ',8(i4))
1008 format(8x,'frz ext ',8(i4))

contains

subroutine error(code)
  integer(kind=iwp), intent(in) :: code
  select case (code)
    case (1)
      write(u6,*) 'input: end of input file encountered'
    case (2)
      write(u6,*) 'input: error while reading input!'
  end select
  write(u6,'(a,a)') 'last command: ',command
  call abend()
end subroutine error

end subroutine gugadrt_mole_inf
