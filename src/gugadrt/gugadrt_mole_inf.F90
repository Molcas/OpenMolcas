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

use Definitions, only: iwp, u6

implicit none
#include "gendrt.fh"
#include "files_gugadrt.fh"
#include "Sysdrt.fh"
#include "mcorb.fh"
#include "refstate.fh"
integer(kind=iwp) :: i, icmd, idisk, im, iml, imr, im_lr_sta, iorb, ispin, itmp, itmpstr(72), j, jcmd, l, ln1, lr, lsmtmp(maxgdm), &
                     ms_ref, mul, nact_sm, nactel, nde, ndisk, ne_act, neact, ngsm, norb1, norb2, norb_all_tmp, ntit
logical(kind=iwp) :: log_debug
character(len=4) :: command
character(len=72) :: line
character(len=132) :: modline
! copy from molcas, bsuo, jun. 30, 2009
integer(kind=iwp), parameter :: ncmd = 18, mxtit = 10
character(len=4), parameter :: cmd(ncmd) = ['TITL','ELEC','SPIN','SYMM','ACTI','PRIN','REFE','FIRS','INAC','CIAL', &
                                            'VALE','INTE','NOCO','ONEO','EXTR','NONI','NACT','END ']

n_electron = 0
nactel = 0
spin = 0.d0
ns_sm = 1
call rdnlst(5,'GUGADRT')
ntit = 0
10 read(5,'(a)',end=991) line
command = line(1:8)
call upcase(command)
if (command(1:1) == '*') goto 10
if (command == ' ') goto 10
jcmd = 0
do icmd=1,ncmd
  if (command == cmd(icmd)) jcmd = icmd
end do
20 goto(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800) jcmd
write(u6,*) 'input: illegal keyword'
write(u6,'(a,a)') 'command=',command
call abend()

!---  process title    command ----------------------------------------*
100 continue
read(5,'(a)',end=991) line
command = line(1:8)
call upcase(command)
if (command(1:1) == '*') goto 100
jcmd = 0
do icmd=1,ncmd
  if (command == cmd(icmd)) jcmd = icmd
end do
if (jcmd /= 0) goto 20
ntit = ntit+1
goto 100

!---  process electron command ----------------------------------------*
200 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 200
read(line,*,err=992) n_electron
goto 10

!---  process spin     command ----------------------------------------*
300 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 300
read(line,*,err=992) ispin
spin = (ispin-1)/2.d0
goto 10

!---  process symmetry command ----------------------------------------*
400 continue
!write (u6,*)'input_guga: keyword symmetry is obsolete and ignored!
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 400
read(line,*,err=992) ns_sm
goto 10

!---  process active   command ----------------------------------------*
500 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 500
modline = line//' 0 0 0 0 0 0 0 0'
read(modline,*,err=992) (nlsm_act(i),i=1,8)
goto 10

!---  process print    command ----------------------------------------*
600 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 600
read(line,*,err=992) iprint
goto 10

!---  process referenc command ----------------------------------------*
700 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 700
read(line,*,err=992) n_ref,ln1
if (n_ref > max_ref) then
  write(u6,*) ' Warnning! Program could not deal with so many reference states!'
  write(u6,*) ' Maximum number of reference states is',max_ref
  write(u6,*) ' Set number of reference states into ',max_ref
  n_ref = max_ref
end if
if (ln1 == 0) then
  logic_mr = .true.
  goto 10
end if
do i=1,n_ref
  read(5,'(80i1)',end=991,err=992) iref_occ(1:ln1,i)
end do
logic_mr = .true.
goto 10

!---  process first    command ----------------------------------------*
800 continue
goto 10

!---  process inactive command ----------------------------------------*
900 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 900
modline = line//' 0 0 0 0 0 0 0 0'
read(modline,*,err=992) (nlsm_dbl(i),i=1,8)
goto 10

!---  process ciall    command ----------------------------------------*
1000 continue
!      read(5,'(a)',end=991) line
!      if ( line(1:1).eq.'*' ) goto 1000
!      read(line,*,err=992) ns_sm
logic_mrelcas = .true.
goto 10

!---  process valence  command ----------------------------------------*
1100 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 1100
modline = line//' 0 0 0 0 0 0 0 0'
read(modline,*,err=992) (nlsm_ext(i),i=1,8)
goto 10

!---  process interact command ----------------------------------------*
1200 continue
goto 10

!---  process nocorr   command ----------------------------------------*
1300 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 1300
modline = line//' 0 0 0 0 0 0 0 0'
read(modline,*,err=992)  !(ncor(i),i=1,8)
!bsuo, jun. 30, 2009 - neglect them
goto 10

!---  process oneocc   command ----------------------------------------*
1400 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 1400
modline = line//' 0 0 0 0 0 0 0 0'
read(modline,*,err=992) ! (ione(i),i=1,8)
!bsuo, jun. 30, 2009 - neglect them
goto 10

!---  process extract  command ----------------------------------------*
1500 write(u6,*) 'input: extract option is redundant and is ignored!'
goto 10

!---  process non-interact command -------------------------------------
1600 continue
write(u6,*) 'input: non-interact option is redundant and is ignored!'
goto 10

!---  process nactel       command -------------------------------------
1700 continue
read(5,'(a)',end=991) line
if (line(1:1) == '*') goto 200
read(line,*,err=992) nactel
goto 10

!---  the end of the input is reached, print the title ----------------*
1800 continue

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
ngw1(1) = 0
ngw2(1) = 0
ngw2(2) = 0
ngw3(1) = 0
ngw3(2) = 0
ngw3(3) = 0
do i=1,norb_all
  ngw2(i+2) = ngw2(i+1)+i
  ngw3(i+3) = ngw3(i+2)+ngw2(i+2)
  ngw4(i+4) = ngw4(i+3)+ngw3(i+3)
end do
nabc = norb_ext-2+ngw2(norb_ext-1)+ngw3(norb_ext)

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

do l=1,norb_inn
  lr = norb_all-l+1
  lsm(lr) = lsm_inn(l)
end do

lr = 0
int_dd_offset(1:8,1:8) = 0
do im=1,ng_sm
  im_lr_sta = 0
  do iml=1,ng_sm
    imr = mul_tab(im,iml)
    if (imr > iml) cycle
    int_dd_offset(iml,imr) = im_lr_sta
    int_dd_offset(imr,iml) = im_lr_sta
    if (iml == imr) im_lr_sta = im_lr_sta+nlsm_ext(iml)*(nlsm_ext(iml)-1)/2
    if (iml /= imr) im_lr_sta = im_lr_sta+nlsm_ext(iml)*nlsm_ext(imr)
  end do
  do l=1,nlsm_ext(im)
    lr = lr+1
    lsm(lr) = im
  end do
end do
lsm_inn(norb_inn+1) = lsm(norb_ext)

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

itmp = 0
do i=1,ng_sm
  ibsm_ext(i) = itmp+1
  itmp = itmp+nlsm_ext(i)
  iesm_ext(i) = itmp
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
  write(u6,'(1x,42a)') 'check input date: number of elctrons error'
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
      if (iref_occ(j,i) == 1) ms_ref = mul_tab(ms_ref,lsm_inn(j))
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
  write(u6,*) 'lsm all',norb_all
  write(u6,1002) (lsm(i),i=norb_all,1,-1)
end if
!*****************************************************

! merge into molcas
! write date into cidrt for ci calculation
noidx = 0
idisk = 0
call idafile(ludrt,1,noidx,2,idisk)
! group symmetry
call idafile(ludrt,1,[ng_sm],1,idisk)
! state symmetry
call idafile(ludrt,1,[ns_sm],1,idisk)
! number of roots to be cal
call idafile(ludrt,1,[mroot],1,idisk)
! number of corelation electrons
call idafile(ludrt,1,[n_electron],1,idisk)
! number of active electrons
call idafile(ludrt,1,[nactel],1,idisk)
! spin symmetry of the state, 2s+1
call idafile(ludrt,1,[ispin],1,idisk)
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
  call idafile(ludrt,1,[2],1,idisk)
  call idafile(ludrt,1,[n_ref],1,idisk)
  ! reference states
  !iref_occ(norb_dz+1:norb_dz+norb_act,i)=itmpstr(1:norb_act)
  do i=1,n_ref
    call idafile(ludrt,1,iref_occ(1,i),norb_inn,idisk)
  end do
end if
if (logic_mrelcas) then
  call idafile(ludrt,1,[4],1,idisk)
end if
! number of ref. states, if logic_mr = .true.
noidx(2) = idisk
idisk = 0
call idafile(ludrt,1,noidx,2,idisk)
noidx = 0
do i=1,ng_sm
  noidx(i) = i
end do
mul = nint(2*spin)+1
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
write(u6,*) '    multiplicity:           ',mul
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

991 continue
write(u6,*) 'input: end of input file encountered'
write(u6,'(a,a)') 'last command: ',command
call abend()
992 continue
write(u6,*) 'input: error while reading input!'
write(u6,'(a,a)') 'last command: ',command
call abend()

1001 format(1x,'norb all group sm')
1002 format(8(1x,i3))
1003 format(16x,8(i4))
1004 format(8x,'frozen  ',8(i4))
1005 format(8x,'double  ',8(i4))
1006 format(8x,'active  ',8(i4))
1007 format(8x,'virtual ',8(i4))
1008 format(8x,'frz ext ',8(i4))

end subroutine gugadrt_mole_inf
