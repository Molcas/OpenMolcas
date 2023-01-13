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

subroutine int_sort()

use gugaci_global, only: intind_abkk, intind_iaqq, intspace_abkk, LuCiInt, maxintseg, norb_dz, norb_ext, norb_frz, norb_inn, &
                         viasum_0, viasum_1, vijkk_0sum, vijkk_1sum, vint_ci
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ia, ia0, idisk, idorbint, idum(1), idx(4), intpos, intposbase, intspace, ivalue, lra, lri, lrk, numb
real(kind=wp) :: etime, stime, time
real(kind=wp), external :: c_time

idx = 0
stime = c_time()

!estimate the largest block of integrals, notice we do not use symmetry
!here
if (norb_ext > norb_inn) then
  ! (ab|cd)
  numb = norb_ext*(norb_ext+1)/2
  numb = numb*(numb+1)/2+norb_ext*numb
else
  ! (ii|ab)
  numb = norb_inn*(norb_inn+1)/2
  numb = numb*norb_ext*(norb_ext+1)/2
end if
! (ia|bc)+(ij|ka)
lra = norb_inn*norb_ext*(norb_ext+1)*norb_ext/2+norb_inn*(norb_inn+1)*norb_inn*norb_ext/2
if (lra > numb) numb = lra

call mma_allocate(vint_ci,numb,label='vint_ci')
vint_ci(:) = Zero
numb = 0
write(u6,900)
call int_index(numb)
#ifdef MOLPRO
call intrd()
#else
call intrd_molcas()
#endif
! --------     vtint(*)-->vint_ci(*)    ------------
numb = 1
call int_sort_inn(numb)         !_ext_0  (ijkl,ijkk,iijj)
idisk = 0
call idafile(luciint,1,idx,4,idisk)
idx(1) = idisk
numb = numb-1
idum(1) = numb
call idafile(luciint,1,idum,1,idisk)
call ddafile(luciint,1,vint_ci,numb,idisk)
write(u6,902) numb
maxintseg = numb

vint_ci(:) = Zero
numb = 1
call int_sort_inn_1(numb)       !_ext_3_1  (iabc and iaqq)
call int_sort_inn_3(numb)       !_ext_1    (ijka)
idx(2) = idisk
numb = numb-1
idum(1) = numb
call idafile(luciint,1,idum,1,idisk)
call ddafile(luciint,1,vint_ci,numb,idisk)
write(u6,903) numb
if (numb > maxintseg) then
  maxintseg = numb
end if

! sum over iaqq
viasum_0(1:norb_inn,1:norb_ext) = Zero
viasum_1(1:norb_inn,1:norb_ext) = Zero
do lri=norb_frz+1,norb_inn
  ia0 = (lri-1)*norb_ext
  do lra=1,norb_ext
    ia = ia0+lra
    intposbase = intind_iaqq(ia)
    if (intposbase == 0) cycle
    do lrk=1,norb_dz
      idorbint = lrk*2-2
      intpos = intposbase+idorbint
      viasum_0(lri,lra) = viasum_0(lri,lra)+vint_ci(intpos)
      viasum_1(lri,lra) = viasum_1(lri,lra)+vint_ci(intpos+1)
    end do
  end do
end do

vint_ci(:) = Zero
numb = 1
call int_sort_inn_2(numb)   !_ext_2_1  (ijcc, ijab, abkk)
idx(3) = idisk
numb = numb-1
idum(1) = numb
call idafile(luciint,1,idum,1,idisk)
call ddafile(luciint,1,vint_ci,numb,idisk)
write(u6,904) numb
if (numb > maxintseg) then
  maxintseg = numb
end if

! sum over (abkk)
intspace = intspace_abkk(1)
intpos = intind_abkk(1)
vijkk_0sum(1:intspace) = Zero
vijkk_1sum(1:intspace) = Zero
do lrk=1,norb_dz
  do ivalue=1,intspace
    vijkk_0sum(ivalue) = vijkk_0sum(ivalue)+vint_ci(intpos)
    vijkk_1sum(ivalue) = vijkk_1sum(ivalue)+vint_ci(intpos+1)
    intpos = intpos+2
  end do
end do

vint_ci(:) = Zero
numb = 1
call int_sort_ext(numb)         !_ext_4_3_2  (abcd,abcc,aabb)
idx(4) = idisk
numb = numb-1
idum(1) = numb
call idafile(luciint,1,idum,1,idisk)
call ddafile(luciint,1,vint_ci,numb,idisk)
write(u6,905) numb
if (numb > maxintseg) then
  maxintseg = numb
end if

idisk = 0
call idafile(luciint,1,idx,4,idisk)

write(u6,910)
call mma_deallocate(vint_ci)

etime = c_time()
time = etime-stime
write(u6,912) time

return

900 format(1x,'start integral sorting',/,1x,'start reading integral file',/)
!901 format(1x,'end reading integral file',/)
902 format(1x,'number of integrals in internal space is ',i8)
903 format(1x,'number of integrals which have one or three indices in external space is ',i8)
904 format(1x,'number of integrals which have two indices in external space is ',i8)
905 format(1x,'number of integrals in external space is ',i8)
910 format(1x,/,1x,'end of integral sorting')
912 format(1x,'total wall clock time for integral sorting=',f8.2,' s',/)

end subroutine int_sort

!subroutine blocks()
!
!use gugaci_global, only: ng_sm, nlsm_all
!use Symmetry_Info, only: Mul
!use Definitions, only: iwp, u6
!
!implicit none
!integer(kind=iwp) :: i, iblktb(5,106), ip, ipq, iq, iqm, ir, irm, is, ism, ispq, ispqr, j, nblock, nint1, nint12, nint2, nintb, &
!                     npp, npq, nrr, nrs
!
!! ----- 1 - electron integrals -----
!nint1 = 0
!do ip=1,ng_sm
!  nint1 = nint1+((nlsm_all(ip)+1)*nlsm_all(ip))/2
!end do
!
!! ----- 2 - electron integrals -----
!nint2 = 0
!nblock = 0
!do ip=1,ng_sm
!  npp = (nlsm_all(ip)+1)*nlsm_all(ip)/2
!  nintb = npp*(npp-1)/2
!  nint2 = nint2+nintb
!  nblock = nblock+1
!  do i=1,4
!    iblktb(i,nblock) = ip
!  end do
!  iblktb(5,nblock) = nint2
!end do
!do ip=1,ng_sm
!  npp = (nlsm_all(ip)+1)*nlsm_all(ip)/2
!  irm = ip
!  do ir=1,irm-1
!    nrr = (nlsm_all(ir)+1)*nlsm_all(ir)/2
!    nintb = npp*nrr
!    nint2 = nint2+nintb
!    nblock = nblock+1
!    iblktb(1,nblock) = ip
!    iblktb(2,nblock) = ip
!    iblktb(3,nblock) = ir
!    iblktb(4,nblock) = ir
!    iblktb(5,nblock) = nint2
!  end do
!end do
!do ip=1,ng_sm
!  iqm = ip
!  do iq=1,iqm-1
!    ipq = nlsm_all(ip)*nlsm_all(iq)
!    nintb = ipq*(ipq-1)/2
!    nint2 = nint2+nintb
!    nblock = nblock+1
!    iblktb(1,nblock) = ip
!    iblktb(2,nblock) = iq
!    iblktb(3,nblock) = ip
!    iblktb(4,nblock) = iq
!    iblktb(5,nblock) = nint2
!  end do
!end do
!do ip=4,8
!  if (nlsm_all(ip) == 0) cycle
!  iqm = ip
!  do iq=1,iqm-1
!    if (nlsm_all(iq) == 0) cycle
!    npq = nlsm_all(ip)*nlsm_all(iq)
!    ispq = Mul(ip,iq)
!    irm = ip
!    do ir=1,irm-1
!      if (nlsm_all(ir) == 0) cycle
!      ispqr = Mul(ispq,ir)
!      ism = ir
!      if (ip == ir) ism = iq
!      do is=1,ism-1
!        if (nlsm_all(is) == 0) cycle
!        if (is /= ispqr) cycle
!        nrs = nlsm_all(ir)*nlsm_all(is)
!        nintb = npq*nrs
!        nint2 = nint2+nintb
!        nblock = nblock+1
!        iblktb(1,nblock) = ip
!        iblktb(2,nblock) = iq
!        iblktb(3,nblock) = ir
!        iblktb(4,nblock) = is
!        iblktb(5,nblock) = nint2
!      end do
!    end do
!  end do
!end do
!nint12 = nint1+nint2
!
!write(u6,100) nint1,nint2,nint12
!!write(12,100) nint1,nint2,nint12
!write(u6,200) nint1
!!write(12,200) nint1
!write(u6,300) (j,(iblktb(i,j),i=1,5),j=1,nblock)
!!write(12,300) (j,(iblktb(i,j),i=1,5),j=1,nblock)
!
!return
!
!100 format(' ',1x/2x,'number of 1-electron integrals  :',i9/2x,'number of 2-electron integrals  :',i9/2x, &
!           'total number of integrals       :',i9)
!200 format(' ',1x/2x,'1-electron blocks  :   1 to',i8/2x,29('*'))
!300 format(' ',1x/2x,'2-electron block description  :'/2x,40('*')/2x,50(3('(',i3,')',4i2,i8,3x)/2x))
!
!end subroutine blocks

!subroutine ff(i,j)
!
!use Definitions, only: iwp
!
!implicit none
!integer(kind=iwp), intent(out) :: i
!integer(kind=iwp), intent(in) :: j
!integer(kind=iwp) :: iq, j0
!
!i = 0
!do iq=1,j
!  j0 = iq*(iq-1)/2+iq
!  if (j0 == j) then
!    i = iq
!    exit
!  end if
!end do
!
!return
!
!end subroutine ff

subroutine int_sort_ext(ii)         !_ext_4_3_2

use gugaci_global, only: ibsm_ext, iesm_ext, ip2_aa_ext_base, ip2_dd_ext_base, ip3_abd_ext_base, ip4_abcd_ext_base, jp2, jp3, &
                         ng_sm, nlsm_ext, norb_ext, norb_number, np3_abd_ext, vint_ci, voint
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: ii
integer(kind=iwp) :: ia, iaend, iasta, ib, ibend, ibsta, ic, icend, icsta, id, idend, idsta, isma, ja, jb, lra, lrb, lrc, lrd, &
                     lsma, lsmb, lsmc, lsmcd, lsmd, lsmtmp(8), na
real(kind=wp), external :: vfutei

! _002_aa_
ip2_aa_ext_base = ii
do lsma=1,ng_sm
  iasta = ibsm_ext(lsma)
  ibsta = iasta+1
  ibend = iesm_ext(lsma)
  do ib=ibsta,ibend
    jb = norb_number(ib)
    do ia=iasta,ib-1
      ja = norb_number(ia)
      vint_ci(ii) = voint(jb,ja)          !(ja|jb)
      !write(10,'(2x,a20,i8,f16.8)') '  int_sort_002',ii,vint_ci(ii)
      ii = ii+1
    end do
  end do
end do
! _002_dd_
ip2_dd_ext_base = ii
!write(10,'(2x,a20,i8)') ' ip2_dd_ext_base=',ii
do ib=2,norb_ext
  jb = norb_number(ib)
  do ia=1,ib-1
    ja = norb_number(ia)
    vint_ci(ii) = voint(ja,jb)   !(ja,jb|ja,jb)
    !write(u6,'(i2,2x,i2)') ja,jb
    !write(u6,'(2x,a20,i8,f16.8)') '  int_sort_002',ii,vint_ci(ii)
    ii = ii+1
  end do
end do
!_003
ip3_abd_ext_base = ii
!write(10,'(2x,a20,i8)') ' ip3_abd_ext_base=',ii

do ic=1,norb_ext
  lrc = norb_number(ic)
  !ip3_abd_ext_base(ic) = ii
  do lsma=1,ng_sm
    iasta = ibsm_ext(lsma)
    !iaend = iesm_ext(lsma)-1
    ibsta = iasta+1
    ibend = iesm_ext(lsma)
    do ib=ibsta,ibend
      lrb = norb_number(ib)
      do ia=iasta,ib-1
        lra = norb_number(ia)
        vint_ci(ii) = vfutei(lra,lrc,lrb,lrc)
        vint_ci(ii+1) = vfutei(lra,lrb,lrc,lrc)
        !write(10,'(2x,4i6,i8,2f16.8)') lra,lrb,lrc,lrc,ii,vint_ci(ii),vint_ci(ii+1)
        ii = ii+2
      end do
    end do
  end do
end do
np3_abd_ext = 0
do isma=1,ng_sm
  na = nlsm_ext(isma)
  if (na <= 0) cycle
  np3_abd_ext = np3_abd_ext+na*(na-1)
end do
!_004
!==================================
!    la<lb<lc<ld
!===================================
! 1. ext _abcd_
!ip4_abcd_ext_base = ii
do lsmd=1,ng_sm
  idsta = ibsm_ext(lsmd)
  idend = iesm_ext(lsmd)
  do lsmc=1,lsmd
    icsta = ibsm_ext(lsmc)
    icend = iesm_ext(lsmc)
    if (lsmc == lsmd) idsta = idsta+1
    lsmcd = Mul(lsmc,lsmd)
    do lsmb=1,lsmc
      lsma = Mul(lsmb,lsmcd)
      if (lsma > lsmb) cycle
      ibsta = ibsm_ext(lsmb)
      ibend = iesm_ext(lsmb)
      iasta = ibsm_ext(lsma)
      iaend = iesm_ext(lsma)
      if (lsmd == lsmb) idsta = idsta+1
      if (lsmc == lsmb) icsta = icsta+1
      if (lsmd == lsma) idsta = idsta+1
      if (lsmc == lsma) icsta = icsta+1
      if (lsmb == lsma) ibsta = ibsta+1
      lsmtmp = lsma+jp2(lsmb)+jp3(lsmc)
      ip4_abcd_ext_base(lsmtmp) = ii
      !write(10,'(2x,a20,2i8)') ' ip4_abcd_ext_base=',lsmtmp,ii
      do id=idsta,idend
        lrd = norb_number(id)
        do ic=icsta,min(id-1,icend)
          lrc = norb_number(ic)
          do ib=ibsta,min(ic-1,ibend)
            lrb = norb_number(ib)
            do ia=iasta,min(ib-1,iaend)
              lra = norb_number(ia)
              vint_ci(ii) = vfutei(lra,lrb,lrc,lrd)       !lra>lrb>lrc>lrd
              vint_ci(ii+1) = vfutei(lra,lrc,lrb,lrd)
              vint_ci(ii+2) = vfutei(lra,lrd,lrc,lrb)
              ii = ii+3
            end do
          end do
        end do
      end do

    end do
  end do
end do

end subroutine int_sort_ext

!===================== ext_2_1 sta =====================================

subroutine int_sort_inn_2(ii)

use gugaci_global, only: ibsm_ext, iesm_ext, intind_abkk, intspace_abkk, lsm_inn, ng_sm, norb_frz, norb_inn, norb_number, vint_ci
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: ii
integer(kind=iwp) :: ibend, ibsta, ira, irb, ismab, lmb, lra, lrb, lri, lrj, lsmi, lsmij, lsmj
real(kind=wp), external :: vfutei

! (ijcc,ijab)
do ismab=1,ng_sm
  !ip2_ab_inn_base(ismab) = ii
  do lri=norb_frz+1,norb_inn-1
    lsmi = lsm_inn(lri)
    do lrj=lri+1,norb_inn
      lsmj = lsm_inn(lrj)
      lsmij = Mul(lsmi,lsmj)
      if (ismab /= lsmij) cycle
      call int_ext_2_1(lri,lrj,lsmij,ii)
    end do
  end do
end do

! (abkk)
do lri=1,norb_inn
  intind_abkk(lri) = ii
  intspace_abkk(lri) = 0
  do lmb=1,ng_sm
    ibsta = ibsm_ext(lmb)
    ibend = iesm_ext(lmb)
    do irb=ibsta,ibend
      lrb = norb_number(irb)
      do ira=ibsta,irb-1
        lra = norb_number(ira)
        intspace_abkk(lri) = intspace_abkk(lri)+1
        vint_ci(ii) = vfutei(lra,lri,lrb,lri)
        vint_ci(ii+1) = vint_ci(ii)-2*vfutei(lra,lrb,lri,lri)
        !write(u6,'(2x,4i6,i8,3f16.8)')lra,lrb,lri,lri,ii,vint_ci(ii),vint_ci(ii+1)
        ii = ii+2
      end do
    end do
  end do
end do

return

end subroutine int_sort_inn_2

subroutine int_ext_2_1(lri,lrj,lsmij,ii)

use gugaci_global, only: ibsm_ext, iesm_ext, intind_ijab, intind_ijcc, intspace_ijab, intspace_ijcc, ng_sm, ngw2, norb_ext, &
                         norb_frz, norb_number, vint_ci
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lrj, lsmij
integer(kind=iwp), intent(inout) :: ii
integer(kind=iwp) :: ic, icend, icsta, id, idend, idsta, ij, jc, jd, lrc, lsmc, lsmd
real(kind=wp), external :: vfutei

!write(10,*) '     start intind_ijab',lri,lrj,ii
ij = lri-norb_frz+ngw2(lrj-norb_frz)
if (lsmij == 1) then       !ijcc
  !write(10,*) '     start intind_ijcc',ii
  intind_ijcc(ij) = ii
  intspace_ijcc(ij) = 0
  do ic=1,norb_ext
    intspace_ijcc(ij) = intspace_ijcc(ij)+1
    lrc = norb_number(ic)
    vint_ci(ii) = vfutei(lrj,lrc,lri,lrc)
    vint_ci(ii+1) = vfutei(lrj,lri,lrc,lrc)
    !write(10,'(2x,4i6,i8,3f16.8)') lrj,lri,lrc,lrc,ii,vint_ci(ii),vint_ci(ii+1)
    ii = ii+2
  end do
end if
intind_ijab(ij) = ii
intspace_ijab(ij) = 0
do lsmc=1,ng_sm
  lsmd = Mul(lsmij,lsmc)
  if (lsmd > lsmc) cycle
  icsta = ibsm_ext(lsmc)
  icend = iesm_ext(lsmc)
  idsta = ibsm_ext(lsmd)
  idend = iesm_ext(lsmd)
  if (lsmc == lsmd) icsta = icsta+1
  do ic=icsta,icend
    jc = norb_number(ic)
    do id=idsta,min(idend,ic-1)
      intspace_ijab(ij) = intspace_ijab(ij)+1
      jd = norb_number(id)
      vint_ci(ii) = vfutei(jd,jc,lrj,lri)
      vint_ci(ii+1) = vfutei(jd,lrj,jc,lri)
      vint_ci(ii+2) = vfutei(jd,lri,lrj,jc)
      !write(10,'(2x,4i6,i8,3f16.8)') jd,jc,lrj,lri,ii,vint_ci(ii),vint_ci(ii+1),vint_ci(ii+2)
      ii = ii+3
    end do
  end do
end do

return

end subroutine int_ext_2_1

!===================== ext_2_1 end =====================================

!===================== ext_1 sta ========================================

subroutine int_sort_inn_3(ii)

use gugaci_global, only: ibsm_ext, iesm_ext, intind_ijka, lsm_inn, ngw2, ngw3, norb_frz, norb_inn, norb_number, vint_ci
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(inout) :: ii
integer(kind=iwp) :: ia, iaend, iasta, ijk, ja, lri, lrj, lrk, lsmd, lsmi, lsmij, lsmj, lsmk
real(kind=wp), external :: vfutei

!write(10,'(2x,a20,i8)') ' int_inn_3_base=',ii
!write(10,*) '     start intind_ijka',ii
do lri=norb_frz+1,norb_inn-2
  lsmi = lsm_inn(lri)
  do lrj=lri+1,norb_inn-1
    lsmj = lsm_inn(lrj)
    lsmij = Mul(lsmi,lsmj)
    do lrk=lrj+1,norb_inn
      lsmk = lsm_inn(lrk)
      lsmd = Mul(lsmij,lsmk)
      ijk = lri-norb_frz+ngw2(lrj-norb_frz)+ngw3(lrk-norb_frz)
      intind_ijka(ijk) = ii
      iasta = ibsm_ext(lsmd)
      iaend = iesm_ext(lsmd)
      do ia=iasta,iaend
        ja = norb_number(ia)
        vint_ci(ii) = vfutei(ja,lrk,lrj,lri)
        vint_ci(ii+1) = vfutei(ja,lrj,lrk,lri)
        vint_ci(ii+2) = vfutei(ja,lri,lrk,lrj)
        !write(10,'(2x,4i6,i8,3f16.8)')ja,lrk,lrj,lri,ii,vint_ci(ii),vint_ci(ii+1),vint_ci(ii+2)
        ii = ii+3
      end do
    end do
  end do
end do

end subroutine int_sort_inn_3

!===================== _ext_1 end ======================================

subroutine int_ext_3_2_1(lri,lsmi,ii)

use gugaci_global, only: ibsm_ext, iesm_ext, intind_iabc, intind_iaqq, nabc, ng_sm, ngw2, ngw3, norb_ext, norb_inn, norb_number, &
                         vint_ci
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lri, lsmi
integer(kind=iwp), intent(inout) :: ii
integer(kind=iwp) :: iabc, iabc0, ib, ib0, ibend, ibsta, icend, icsta, idend, idsta, irb, irc, ird, lrb, lrc, lrd, lsmb, lsmbc, &
                     lsmc, lsmd, nb, nc, nd
real(kind=wp), external :: vfutei

!write(10,*) '   int_ext_3_2_1'
iabc0 = (lri-1)*nabc
!write(10,*) '     start intind_iabc',ii
do lsmd=1,ng_sm
  lsmbc = Mul(lsmi,lsmd)
  do lsmc=1,lsmd
    lsmb = Mul(lsmbc,lsmc)
    if (lsmb > lsmc) cycle
    idsta = ibsm_ext(lsmd)
    idend = iesm_ext(lsmd)
    icsta = ibsm_ext(lsmc)
    icend = iesm_ext(lsmc)
    ibsta = ibsm_ext(lsmb)
    ibend = iesm_ext(lsmb)
    if (lsmb == lsmc) icsta = icsta+1
    if (lsmc == lsmd) idsta = idsta+1
    if (lsmb == lsmd) idsta = idsta+1
    do ird=idsta,idend
      lrd = norb_number(ird)
      nd = ngw3(ird)
      do irc=icsta,min(icend,ird-1)
        lrc = norb_number(irc)
        nc = ngw2(irc)
        do irb=ibsta,min(ibend,irc-1)
          lrb = norb_number(irb)
          nb = irb
          iabc = iabc0+nb+nc+nd
          intind_iabc(iabc) = ii
          vint_ci(ii) = vfutei(lrb,lrc,lrd,lri)
          vint_ci(ii+1) = vfutei(lrb,lrd,lrc,lri)
          vint_ci(ii+2) = vfutei(lrb,lri,lrd,lrc)
          ii = ii+3
        end do
      end do
    end do
  end do
end do

!write(10,*) '     start intind_iaqq',ii
ib0 = (lri-1)*norb_ext
lsmb = lsmi                       !3_ibqq
ibsta = ibsm_ext(lsmb)
ibend = iesm_ext(lsmb)
do irb=ibsta,ibend
  lrb = norb_number(irb)
  ib = ib0+irb
  intind_iaqq(ib) = ii
  do lrc=1,norb_inn
    vint_ci(ii) = vfutei(lrb,lrc,lri,lrc)
    vint_ci(ii+1) = vfutei(lrb,lri,lrc,lrc)
    ii = ii+2
  end do
  do irc=norb_ext,1,-1
    lrc = norb_number(irc)
    vint_ci(ii) = vfutei(lrb,lrc,lri,lrc)
    vint_ci(ii+1) = vfutei(lrb,lri,lrc,lrc)
    !write(u6,'(2x,4i6,i8,3f16.8)')lrb,lri,lrc,lrc,ii,vint_ci(ii),vint_ci(ii+1)
    ii = ii+2
  end do
end do

end subroutine int_ext_3_2_1

subroutine int_sort_inn_1(ii)

use gugaci_global, only: lsm_inn, norb_frz, norb_inn
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: ii
integer(kind=iwp) :: lri, lsmi

call determine_para_array_for_int1ind()

do lri=norb_frz+1,norb_inn
  lsmi = lsm_inn(lri)
  call int_ext_3_2_1(lri,lsmi,ii)
end do

end subroutine int_sort_inn_1

!===================== ext_3_2_1 end ===================

function list3(i,j,k)

use gugaci_global, only: loij, ngw2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: list3
integer(kind=iwp), intent(in) :: i, j, k
integer(kind=iwp) :: nij

nij = i+ngw2(j)
list3 = loij(nij)+2*(k-1)

return

end function list3

function list4(ld,lc,lb,la)

use gugaci_global, only: loijk, ncibl, ngw2, ngw3
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: list4
integer(kind=iwp), intent(in) :: ld, lc, lb, la
integer(kind=iwp) :: lra, njkl

lra = ncibl(la)
njkl = ld+ngw2(lc)+ngw3(lb)
list4 = loijk(njkl)+3*(lra-1)

return

end function list4

subroutine int_sort_inn(numb)

use gugaci_global, only: loij, loijk, lsm_inn, ncibl, ngw2, ngw3, norb_inn, vint_ci
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: numb
integer(kind=iwp) :: i, j, k, la, lb, lc, ld, lra, ms, msa, msb, msc, mscd, msd, msob(120), nij, njkl, nolra
real(kind=wp), external :: vfutei

msob = 0
do lra=norb_inn,1,-1
  ms = lsm_inn(lra)
  msob(ms) = msob(ms)+1
  ncibl(lra) = msob(ms)
end do
!=======================================================================
!write(u6,'(1x,14i3)') (ncibl(i),i=1,norb_inn)

numb = 1
do i=1,norb_inn-1
  do j=i+1,norb_inn
    if (lsm_inn(i) /= lsm_inn(j)) cycle

    nij = i+ngw2(j)
    loij(nij) = numb
    do k=1,norb_inn
      vint_ci(numb) = vfutei(j,k,i,k)
      vint_ci(numb+1) = vfutei(j,i,k,k)
      !write(u6,'(2x,4i6,i8,3f16.8)') i,j,k,k,numb,vint_ci(numb),vint_ci(numb+1)
      numb = numb+2
    end do
  end do
end do
!write(u6,*) 'num_inn',numb
!write(u6,*) loij(1:8)
!call abend()

!=======================================================================
! la<lb<lc<ld
do ld=1,norb_inn-3
  do lc=ld+1,norb_inn-2
    msd = lsm_inn(ld)
    msc = lsm_inn(lc)
    mscd = Mul(msd,msc)
    do lb=lc+1,norb_inn-1
      msb = lsm_inn(lb)
      msa = Mul(mscd,msb)

      njkl = ld+ngw2(lc)+ngw3(lb)
      loijk(njkl) = numb

      nolra = 0
      do la=norb_inn,lb+1,-1
        if (lsm_inn(la) /= msa) cycle
        nolra = nolra+1
        !list = loijk(njkl)+3*(nolra-1)

        !write(u6,'(2x,4i3,2i7)')  la,lb,lc,ld,list,numb

        vint_ci(numb) = vfutei(la,lc,lb,ld)        !tmp stop
        vint_ci(numb+1) = vfutei(la,lb,lc,ld)
        vint_ci(numb+2) = vfutei(la,ld,lc,lb)
        !write(10,'(2x,4i6,i8,3f16.8)')la,lb,lc,ld,numb,vint_ci(numb),vint_ci(numb+1),vint_ci(numb+2)
        numb = numb+3
      end do
    end do
  end do
end do

return
!=======================================================================

end subroutine int_sort_inn
