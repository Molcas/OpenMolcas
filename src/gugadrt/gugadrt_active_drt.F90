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

subroutine gugadrt_active_drt()

use gugadrt_global, only: iseg_sta, iseg_downwei, jpad_upwei, jd, jj, js, jt, jv, logic_mr, logic_mrelcas, max_node, mxnode, &
                          nci_dim, ng_sm, no, norb_act, norb_inn, nu_ad, nu_ae
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, im, indd, iseg_dim(25), jde, jdim, jdn, jds, jji, jp, jpe, jpn, jsim, jtim, ndd, ndi, ndim
integer(kind=iwp), allocatable :: iin(:)

iseg_dim = 0
nci_dim = 0
if (norb_act == 0) then
  !====================  norb_act=0 ========================
  iseg_sta(1) = 0
  iseg_dim(1) = 1
  do im=1,ng_sm
    jdim = nu_ae(1+im)
    jtim = nu_ae(9+im)
    jsim = nu_ae(17+im)
    jd(im) = jdim
    jt(im) = jtim
    js(im) = jsim
    iseg_dim(jdim) = iseg_downwei(jdim)*jpad_upwei(jdim)
    iseg_dim(jtim) = iseg_downwei(jtim)*jpad_upwei(jtim)
    iseg_dim(jsim) = iseg_downwei(jsim)*jpad_upwei(jsim)
    if (iseg_dim(jdim) == 0) then
      jd(im) = 0
      nu_ad(jdim) = 0
      nu_ae(jdim) = 0
    end if
    if (iseg_dim(jtim) == 0) then
      jt(im) = 0
      nu_ad(jtim) = 0
      nu_ae(jtim) = 0
    end if
    if (iseg_dim(jsim) == 0) then
      js(im) = 0
      nu_ad(jsim) = 0
      nu_ae(jsim) = 0
    end if
  end do
  do jp=2,25
    iseg_sta(jp) = iseg_sta(jp-1)+iseg_dim(jp-1)
  end do
  nci_dim = iseg_sta(25)+iseg_dim(25)
  iseg_sta(26) = nci_dim
else
  !====================  norb_act<>0 ========================
  if (logic_mr) call gugadrt_rst(ndd,indd)
  if (logic_mrelcas) call gugadrt_rcas(ndd,indd) !npp=3

  nu_ae(1) = jv
  do im=1,ng_sm
    nu_ae(1+im) = jd(im)
    nu_ae(9+im) = jt(im)
    nu_ae(17+im) = js(im)
  end do
  jds = 1
  jde = mxnode
  jpe = no(norb_inn+1)

  call mma_allocate(iin,[0,max_node],label='iin')
  iin(1:max_node) = 0
  jp = jv
  iin(1:jpe) = 0
  iin(0) = 0
  iin(jp) = 1
  iseg_sta(1) = 0 ! for node_ext
  do jpn=jpe,1,-1
    do i=1,4
      jji = jj(i,jpn)
      if (iin(jji) == 0) cycle
      iin(jpn) = iin(jpn)+iin(jji)
    end do
  end do

  do jdn=jds,jde
    if (nu_ad(jdn) == 0) cycle
    ndi = iin(jdn)*iseg_downwei(1)*jpad_upwei(jdn)
    iseg_dim(1) = iseg_dim(1)+ndi
  end do

  do im=1,ng_sm
    jp = jd(im)
    iseg_sta(1+im) = nci_dim ! for node_d_ext
    if (jp == 0) cycle
    iin(1:jpe) = 0
    iin(0) = 0
    iin(jp) = 1
    do jpn=jpe,1,-1
      do i=1,4
        jji = jj(i,jpn)
        if (iin(jji) == 0) cycle
        iin(jpn) = iin(jpn)+iin(jji)
      end do
    end do
    do jdn=jds,jde
      if (nu_ad(jdn) == 0) cycle
      ndi = iin(jdn)*iseg_downwei(1+im)*jpad_upwei(jdn)
      iseg_dim(1+im) = iseg_dim(1+im)+ndi
    end do
  end do

  do im=1,ng_sm
    jp = jt(im)
    iseg_sta(9+im) = nci_dim ! for node_t_ext
    if (jp == 0) cycle
    iin(1:jpe) = 0
    iin(0) = 0
    iin(jp) = 1
    do jpn=jpe,1,-1
      do i=1,4
        jji = jj(i,jpn)
        if (iin(jji) == 0) cycle
        iin(jpn) = iin(jpn)+iin(jji)
      end do
    end do
    do jdn=jds,jde
      if (nu_ad(jdn) == 0) cycle
      ndi = iin(jdn)*iseg_downwei(9+im)*jpad_upwei(jdn)
      iseg_dim(9+im) = iseg_dim(9+im)+ndi
    end do
  end do
  do im=1,ng_sm
    jp = js(im)
    iseg_sta(17+im) = nci_dim ! for node_s_ext
    if (jp == 0) cycle
    iin(1:jpe) = 0
    iin(0) = 0
    iin(jp) = 1
    do jpn=jpe,1,-1
      do i=1,4
        jji = jj(i,jpn)
        if (iin(jji) == 0) cycle
        iin(jpn) = iin(jpn)+iin(jji)
      end do
    end do
    do jdn=jds,jde
      if (nu_ad(jdn) == 0) cycle
      ndi = iin(jdn)*iseg_downwei(17+im)*jpad_upwei(jdn)
      iseg_dim(17+im) = iseg_dim(17+im)+ndi
    end do
  end do
  do jp=2,25
    iseg_sta(jp) = iseg_sta(jp-1)+iseg_dim(jp-1)
  end do
  nci_dim = iseg_sta(25)+iseg_dim(25)
  iseg_sta(26) = nci_dim
  call mma_deallocate(iin)
end if
! to the end of dbl,act,ext parts
call gugadrt_dbl_downwalk()
!write(u6,*) '  end of drt,nci_dim= ',nci_dim
!write(u6,*) 'number of cfss: ',nci_dim
write(u6,*)
write(u6,*) '-----------------------------------------------'
write(u6,*) '    csf information'
write(u6,*) '    num. of configurations:        ',nci_dim
ndim = iseg_dim(1)
write(u6,*) '    num. of valence states:        ',ndim
ndim = 0
do i=2,9
  ndim = ndim+iseg_dim(i)
end do
write(u6,'(5x,a32,1x,i12)') 'num. of doublet couple singles: ',ndim
ndim = 0
do i=10,17
  ndim = ndim+iseg_dim(i)
end do
write(u6,'(5x,a32,1x,i12)') 'num. of triplet couple doubles: ',ndim
ndim = 0
do i=18,25
  ndim = ndim+iseg_dim(i)
end do
write(u6,'(5x,a32,1x,i12)') 'num. of singlet couple doubles: ',ndim
write(u6,*) '-----------------------------------------------'

return

end subroutine gugadrt_active_drt
