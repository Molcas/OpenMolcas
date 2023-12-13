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
! Copyright (C) 2007, Bingbing Suo                                     *
!***********************************************************************

! 28 dec 2007 -bsuo- the initial vectors for davidson diagonalization are
!                    the vectors obtained by diagnal h0 space will be us
!                    initial vector for mrci calculation.
! 18 jul 2007 -bsuo- multi-root calculation is revised to perform calcul
!                    nci_dim*mroot>max_civector, the roots are calculated
!                    by one, old subroutine for davidson diagonalization
!                    deleted. the subroutines to init the initial vector
!                    revised.
! 15 mar 2007 -bsuo- new method for multi-root calculation finished
! 23 feb 2007 -bsuo- revised by suo bing for multi-root calculation

! this module contains some routines to calculate the segmentation
! value of the partial loops, information of the molecular orbital
! and the davidson diagonalization calculation

subroutine mole_inf()

use gugaci_global, only: cm_cri, logic_calpro, logic_inivec_read, maxciiter, mroot, pror, vthrealp, vthreen, vthreresid
                         !, logic_tdav
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp), parameter :: ncmd = 9
integer(kind=iwp) :: icmd, istatus, jcmd, ntit
character(len=4) :: command
character(len=72) :: line
#ifndef MOLPRO
#define _END_ 'END '
#endif
#ifdef _XIANEST_
#define _END_ '$END'
#endif
logical(kind=iwp) :: skip
character(len=4), parameter :: cmd(ncmd) = ['TITL','NRRO','MAXI','CPRO','PTHR','CONV','PROR','REST',_END_]

#ifndef MOLPRO
call rdnlst(u5,'GUGACI')
#endif

! init some program control logic variables
! traditional davidson diagonalization method is used
!logic_tdav = .true.
logic_inivec_read = .false.
logic_calpro = .false.
cm_cri = 0.05_wp
pror = 1.0e-4_wp
! set the default convergence threshold
vthreen = 1.0e-8_wp
vthrealp = 1.0e-6_wp
vthreresid = 1.0e-8_wp
mroot = 1
maxciiter = 30
ntit = 0
skip = .false.
jcmd = 0
do
  if (.not. skip) then
    read(u5,'(a)',iostat=istatus) line
    if (istatus < 0) call error(istatus)
    command = line(1:8)
    call upcase(command)
    if (command(1:1) == '*') cycle
#   ifdef _XIANEST_
    if (command(1:4) == '$MRC') cycle
#   endif
    if (command == ' ') cycle
    jcmd = 0
    do icmd=1,ncmd
      if (command == cmd(icmd)) jcmd = icmd
    end do
  end if
  skip = .false.
  select case (jcmd)
    case (1)
      !---  process title    command -----------------------------------
      do
        read(u5,'(a)',iostat=istatus) line
        if (istatus < 0) call error(istatus)
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
      !---  process nrroot command -------------------------------------
      do
        read(u5,'(a)',iostat=istatus) line
        if (istatus < 0) call error(istatus)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=istatus) mroot
      if (istatus > 0) call error(istatus)

    case (3)
      !---  process Maxiterations command ------------------------------
      do
        read(u5,'(a)',iostat=istatus) line
        if (istatus < 0) call error(istatus)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=istatus) maxciiter
      if (istatus > 0) call error(istatus)

    case (4)
      !---  process Calculate property command -------------------------
      logic_calpro = .true.

    case (5)
      !---  process Thresh print command -------------------------------
      do
        read(u5,'(a)',iostat=istatus) line
        if (istatus < 0) call error(istatus)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=istatus) cm_cri
      if (istatus > 0) call error(istatus)

    case (6)
      !---  process Convergence threshold command ----------------------
      do
        read(u5,'(a)',iostat=istatus) line
        if (istatus < 0) call error(istatus)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=istatus) vthreen,vthrealp,vthreresid
      if (istatus > 0) call error(istatus)

    case (7)
      !---  process print orbital command ------------------------------
      do
        read(u5,'(a)',iostat=istatus) line
        if (istatus < 0) call error(istatus)
        if (line(1:1) /= '*') exit
      end do
      read(line,*,iostat=istatus) pror
      if (istatus > 0) call error(istatus)

    case (8)
      !---  process restart command ------------------------------------

    case (9)
      !--- End of GUGACI input -----------------------------------------
      exit

    case default
      write(u6,*) 'input: illegal keyword'
      write(u6,'(a,a)') 'command=',command
#     ifndef MOLPRO
      call abend()
#     endif
  end select

end do

call mole_inf_molcas()

return

contains

subroutine error(code)

  integer(kind=iwp), intent(in) :: code

  if (code < 0) then
    write(u6,*) 'input: end of input file encountered'
  else
    write(u6,*) 'input: error while reading input!'
  end if
  write(u6,'(a,a)') 'last command: ',command
# ifndef MOLPRO
  call abend()
# endif

end subroutine error

end subroutine mole_inf

subroutine mole_inf_molcas()

use gugaci_global, only: ibsm_ext, iesm_ext, int_dd_offset, iref_occ, logic_assign_actorb, logic_mr, lsm, lsm_inn, lsmorb, LuDrt, &
                         n_ref, nabc, ng_sm, ngw2, ngw3, ngw4, nlsm_all, nlsm_bas, nlsm_dbl, nlsm_ext, nlsm_frz, noidx, norb_act, &
                         norb_all, norb_dbl, norb_dz, norb_ext, norb_frz, norb_inn, ns_sm, nstart_act, spin !, logic_mrelcas
use Symmetry_Info, only: Mul
use Constants, only: Half
use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: i, idisk, idum(1), idx, im, im_lr_sta, iml, imr, imrcas_case, iorb, ispin, itmp, j, l, lr, nact_sm, &
                     nlsm_act(mxSym), nlsm_inn(mxSym), lsmtmp(mxSym), ngsm, ni, norb_all_tmp

!open(nf1,file='drt.inp')
!read(nf1,*)
!read(nf1,*) mroot,mth_eigen,kin
!read(nf1,*) n_electron,spin,ng_sm,ns_sm,cm_cri
!read(nf1,*) nlsm_all(1:ng_sm)
!read(nf1,*) nlsm_frz(1:ng_sm)
!read(nf1,*) nlsm_dbl(1:ng_sm)
!read(nf1,*) nlsm_act(1:ng_sm)
!read(nf1,*) log_thre
! merge into molcas
! write date into cidrt for ci calculation
noidx = 0
idisk = 0
call idafile(ludrt,2,noidx,2,idisk)
! group symmetry
call idafile(ludrt,2,idum,1,idisk)
ng_sm = idum(1)
! state symmetry
call idafile(ludrt,2,idum,1,idisk)
ns_sm = idum(1)
! number of roots to be cal
call idafile(ludrt,2,idum,1,idisk)
! number of corelation electrons
call idafile(ludrt,2,idum,1,idisk)
!n_electron = idum(1)
! number of active electrons
call idafile(ludrt,2,idum,1,idisk)
!nactel = idum(1)
! spin symmetry of the state, 2s+1
call idafile(ludrt,2,idum,1,idisk)
ispin = idum(1)
! dbl orb
call idafile(ludrt,2,nlsm_dbl,8,idisk)
! act orb
call idafile(ludrt,2,nlsm_act,8,idisk)
! all correlated orb
call idafile(ludrt,2,nlsm_all,8,idisk)
! num. basis
call idafile(ludrt,2,nlsm_bas,8,idisk)
spin = (ispin-1)*Half
nlsm_frz(:) = 0

norb_frz = 0
norb_dbl = 0
norb_act = 0
norb_dz = 0
norb_all = 0
noidx = 0
idx = 0
ni = 0
do i=1,ng_sm
  norb_frz = norb_frz+nlsm_frz(i)
  norb_dbl = norb_dbl+nlsm_dbl(i)
  norb_act = norb_act+nlsm_act(i)
  nlsm_inn(i) = nlsm_frz(i)+nlsm_dbl(i)+nlsm_act(i)
  nlsm_ext(i) = nlsm_all(i)-nlsm_inn(i)
  norb_all = norb_all+nlsm_all(i)
  noidx(i) = idx
  idx = idx+nlsm_all(i)
  lsmorb(ni+1:idx) = i
  ni = idx
end do
norb_inn = norb_frz+norb_dbl+norb_act
norb_ext = norb_all-norb_inn
norb_dz = norb_dbl+norb_frz

nstart_act = norb_dz+1
ngw2(:) = 0
ngw3(:) = 0
ngw4(:) = 0
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
# ifndef MOLPRO
  call abend()
# endif
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
    imr = Mul(im,iml)
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

logic_mr = .false.
!logic_mrelcas = .false.
logic_assign_actorb = .false.

call idafile(ludrt,2,idum,1,idisk)
imrcas_case = idum(1)
if (imrcas_case == 2) then
  logic_mr = .true.
  call idafile(ludrt,2,idum,1,idisk)
  n_ref = idum(1)
  do i=1,n_ref
    call idafile(ludrt,2,iref_occ(1,i),norb_inn,idisk)
  end do
end if
!if (imrcas_case == 4) then
!  logic_mrelcas = .true.
!end if

do i=1,ng_sm   !norb_inn
  lsmtmp(i) = 0
end do
do i=1,norb_inn
  lsmtmp(lsm_inn(i)) = lsmtmp(lsm_inn(i))+1
end do

itmp = 0
ibsm_ext(:) = 0
iesm_ext(:) = 0
do i=1,ng_sm
  ibsm_ext(i) = itmp+1
  itmp = itmp+nlsm_ext(i)
  iesm_ext(i) = itmp
end do

!***********************************************************************
#ifdef _DEBUG
idebug = 1
if (idebug == 1) then
  write(u6,1001)
  write(u6,1002) norb_all,norb_frz,norb_dz,norb_inn,norb_ext
  write(u6,1002) nlsm_all(1:ng_sm)
  write(u6,1002) nlsm_frz(1:ng_sm)
  write(u6,1002) nlsm_dbl(1:ng_sm)
  write(u6,1002) nlsm_inn(1:ng_sm)
  write(u6,1002) nlsm_ext(1:ng_sm)
  write(u6,*) 'lsm inn, imrcas_case',imrcas_case
  write(u6,1002) lsm_inn(1:norb_inn)
  write(u6,*) 'lsm all',norb_all
  write(u6,1002) (lsm(i),i=norb_all,1,-1)
end if
1001 format(1x,'all frz dz inn ext')
1002 format(8(1x,i3))
#endif
!***********************************************************************

return

end subroutine mole_inf_molcas
