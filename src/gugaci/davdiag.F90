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
! Copyright (C) 2007,2009, Bingbing Suo                                *
!***********************************************************************
! Jul. 3, 2009 -bsuo- subroutines are used in davidson diagonalization

subroutine cidiagonalize(mxvec)
!******************************************************
! this subroutine does diagonalization of
! ci matrix for multi-root mrci program
! 26 feb 2007 - write by suo bing
!------------------------------------------------------

#include "drt_h.fh"
#include "files_gugaci.fh"
#include "scratch.fh"
#include "thresh.fh"
dimension idxvec(max_iter)
dimension valpha(max_root), vcien(max_root), vcienold(max_root), difeci(max_root), vresid(max_root)
!dimension vad(max_vector),th(max_vector)
logical log_convergence, log_muliter, restart
data dzero/0.d0/,depc/1.0d-7/
data vortho_criterion/1.d-8/
!data venergy_criterion/1.d-8/,valpha_criterion/1.d-7/,vresid_criterion/1.d-8/

venergy_criterion = vthreen
valpha_criterion = vthrealp
vresid_criterion = vthreresid
write(6,*) ' threshold for convergence is set as'
write(6,*) venergy_criterion,valpha_criterion,vresid_criterion

!write(6,*) vthreen,vthrealp,vthreresid
log_convergence = .false.
log_muliter = .false.
if (mroot*nci_dim <= mxvec) log_muliter = .true.
mcroot = mroot
!*********************************************************************
! for debug use only

id = 1
if (id == 2) then

  indx(1) = 0
  indx(2) = nci_dim
  vector1(1:nci_dim) = 0.d0
  vector2(1:nci_dim) = 0.d0
  !call cielement()
  vector1 = 0.d0; vector2 = 0.d0
  vector1(55) = 1.d0
  call readint(2,vint_ci)
  call vd_drt_ci_new()
  call dv_drt_ci_new()
  !call matrix_vector_multi_parallel_drt(sechc)
  do i=1,3450
    write(11,'(i8,1x,f18.8)') i,vector2(i)
  end do
  call abend()
end if
!***********************************************************************

mth_eigen = 0
numroot = mroot
do
  restart = .false.
  sc0 = c_time()
  vector1 = 0.d0
  vector2 = 0.d0
  if (log_muliter) then
    msroot = 1
    mcroot = mroot
    call mrcibasis(nci_dim,mroot,mjn,indx,vector1,vector2,vcien,mth_eigen,mroot)
    !call mrcibasis_init(nci_dim,mroot,mjn,indx,vector1,vector2,vcien,mth_eigen,mroot)
  else
    msroot = 1
    mcroot = 1
    if (mth_eigen == 0) then
      mth_eigen = 1
    else
      mth_eigen = mth_eigen+1
    end if
    if (mth_eigen == 1) then
      call mrcibasis(nci_dim,mroot,mjn,indx,vector1,vector2,vcien,mth_eigen,1)
    else
      msroot = 1
      mcroot = 1
      call mrcibasis_rest(nci_dim,numroot,mjn,indx,vector1,vector2,vcien,mth_eigen,1)
    end if
  end if
  !call abend()
  sc1 = c_time()
  sct = sc1-sc0
  write(6,890) 2*mroot,sct

  idxvec(1) = mroot
  idxvec(2) = 2*mroot
  idxvec(3:max_iter) = 0
  vcienold(1:max_root) = dzero
  vresid(1:max_root) = dzero
  valpha(1:max_root) = dzero
  iiter = 2
  kval = 2*mroot
  msroot = 1
  mcroot = mroot
  mtsta = 1
  iciter = 0
  irset = 0

  write(6,901)
  do while (.not. log_convergence)
    iciter = iciter+1
    sc0 = c_time()
    call matrmkmul(kval,mtsta,iiter,idxvec,irset)
    call hotred(max_kspace,kval,vp,vd,ve,vu)
    call qlcm(max_kspace,kval,vd,ve,vu)
    !if (iciter == 3) call abend()

    vcienold(1:mroot) = vcien(1:mroot)
    valpha(mtsta:mroot) = dzero
    do m=mtsta,mroot
      vcien(m) = vd(m)
      valpha(m) = abs(vu(kval,m))
      difeci(m) = abs(vcien(m)-vcienold(m))
    end do

    call compute_residual_vector_mroot(kval,mtsta,iiter,idxvec,vresid,vcien)
    !write(6,900) iciter
    !write(6,901)
    if (iciter == 1) then
      sechc = 0.d0
      scvp = 0.d0
    end if
    do mt=1,mroot
      write(6,902) iciter,mt,vcien(mt),difeci(mt),valpha(mt),vresid(mt),sechc,scvp
    end do
    call xflush(6)

    !write(6,903) sechc,scvp
    !write(nf2,903) sechc,scvp
    mtsta0 = mtsta
    do mt=mtsta0,mroot
      if (((valpha(mt) < valpha_criterion) .or. (vresid(mt) < vresid_criterion)) .and. (abs(difeci(mt)) < venergy_criterion)) then
        if (mt == mtsta) then
          mtsta = mtsta+1
          mcroot = mcroot-1
        end if
        if (mtsta == mroot+1) then
          if (log_muliter) then
            log_convergence = .true.
            call get_eigvector(mtsta0,vcien,valpha,diffci,vresid,mth_eigen,log_muliter)
            return
          else
            if (mth_eigen < numroot) then
              call get_eigvector(mtsta0,vcien,valpha,diffci,vresid,mth_eigen,log_muliter)
              log_convergence = .false.
              restart = .true.
            else
              log_convergence = .true.
              call get_eigvector(mtsta0,vcien,valpha,diffci,vresid,mth_eigen,log_muliter)
              return
            end if
          end if
        end if
      end if
    end do

    if (iciter > maxciiter) then
      write(6,*) ' warning! mrci fail to converged! program stop!'
      write(6,'(a30,1x,i3)') ' number of converged roots is=',mtsta0
      mtsta0 = mroot
      call get_eigvector(mtsta0,vcien,valpha,diffci,vresid,mth_eigen,log_muliter)
      return
    end if

    if (mtsta /= mtsta0) then
      nd = nci_dim*(mroot-mtsta0+1)
      call read_ml(lucidia,1,vector2,nd,2)
      ! write converged cm into file 7
      do mt=mtsta0,mtsta-1
        mtidx = indx(mt-mtsta0+1)
        call write_ml(lucivec,1,vector2(mtidx+1:mtidx+nci_dim),nci_dim,mt)
      end do
      nd = (mroot-mtsta0+1)*nci_dim
      vector2(1:nd) = vector1(1:nd)
      vector1(1:nd) = dzero
      mi = indx(mtsta-mtsta0+1)
      nd = (mroot-mtsta+1)*nci_dim
      do i=1,nd
        vector1(i) = vector2(i+mi)
      end do
    end if
    !write(6,910) mtsta-1

    !*******************************************************************
    ! reset kspace

    if (kval+mcroot > max_kspace-1) then
      write(6,911)

      nd = mroot*nci_dim
      vector1(1:nd) = dzero
      nd = (mroot-mtsta0+1)*nci_dim
      call read_ml(lucidia,1,vector2,nd,2)
      nda = nci_dim*mroot
      call read_bv(lucitv1,1,vector1,nda)
      !rewind(7)
      do mt=1,mtsta-1
        mtidx = indx(mt)
        call read_ml(lucivec,1,vector1(mtidx+1:mtidx+nci_dim),nci_dim,mt)
        !read(7) vector1(mtidx+1:mtidx+nci_dim)
      end do

      nda = nci_dim*mroot
      nd = (mroot-mtsta+1)*nci_dim
      mtidx = indx(mtsta)
      mid = indx(mtsta-mtsta0+1)
      vector1(mtidx+1:mtidx+nd) = vector2(mid+1:mid+nd)
      call write_bv(lucitv1,1,vector1,nda)
      do mt=1,mtsta-1
        vet = vcien(mt)
        mtidx = indx(mt)
        do l=1,nci_dim
          vector1(l+mtidx) = vet*vector1(l+mtidx)
        end do
      end do

      nda = mroot*nci_dim
      vector1(1:nda) = dzero
      do jiter=1,iiter
        if (jiter == 1) then
          irts = 1
        else
          irts = idxvec(jiter-1)+1
        end if
        irte = idxvec(jiter)
        nd = (irte-irts+1)*nci_dim
        call read_bv(lucitv2,jiter,vector2,nd)
        do mt=mtsta,mroot
          mtidx = indx(mt)
          do irot=irts,irte
            itidx = indx(irot-irts+1)
            vuim = vu(irot,mt)
            do l=1,nci_dim
              vector1(mtidx+l) = vector1(mtidx+l)+vuim*vector2(itidx+l)
            end do
          end do
        end do
      end do

      nda = mroot*nci_dim
      vector2(1:nda) = vector1(1:nda)
      call write_bv(lucitv2,1,vector2,nda)

      ! start new trial vector
      nda = nci_dim*mroot
      nd = nci_dim*(mroot-mtsta+1)
      mtidx = indx(mtsta)
      vector2(1:nda) = dzero
      vector2(1:nd) = vector1(mtidx+1:mtidx+nd)

      nd = (mroot-mtsta0+1)*nci_dim
      call read_ml(lucidia,1,vector1,nd,2)
      !rewind(nf22)
      !read(nf22) vector1(1:nd)
      do mt=mtsta,mroot
        mtidx = indx(mt-mtsta+1)
        mtid = indx(mt-mtsta0+1)
        do l=1,nci_dim
          vector1(l+mtidx) = vector1(mtid+l)
        end do
      end do

      do mt=mtsta,mroot
        mtidx = indx(mt-mtsta+1)
        venergy = vcien(mt)
        do l=1,nci_dim
          vector1(mtidx+l) = venergy*vector1(mtidx+l)-vector2(mtidx+l)
        end do
      end do

      !call read_bv(nf8,1,vector2,nci_dim)
      call read_ml(lucidia,1,vector2,nci_dim,1)
      do mt=mtsta,mroot
        mtidx = indx(mt-mtsta+1)
        do l=1,nci_dim
          depff = vector2(l)-vcien(mt)
          if (abs(depff) <= depc) depff = depc
          vector1(mtidx+l) = vector1(mtidx+l)/depff
        end do
      end do

      iiter = 1
      kval = mroot
      mcroot = mroot-mtsta+1
      call orthog(kval,iiter,mtsta,idxvec)

      idxvec(iiter) = kval

      mcroot = mroot-mtsta+1
      iiter = 2
      if (.not. log_muliter) then
        if (mth_eigen > 1) then
          call orthogwconvec()
        end if
      end if
      nd = (mroot-mtsta+1)*nci_dim
      call write_bv(lucitv1,iiter,vector1,nd)
      idxvec(iiter) = kval
      mn = mroot*(mroot+1)/2
      vp(1:mn) = dzero
      mn = 0
      do m=1,mroot
        mn = mn+m
        vp(mn) = vcien(m)
      end do
      irset = 1
    else

      !*****************************************************************
      ! compute revised new appoximate vector
      !
      !call read_ml(lucidia,1,vector2,nci_dim,1)
      !irset=1
      !do mt=mtsta,mroot
      !  mtidx = indx(mt-mtsta+1)
      !  !logic_tda v= .false.
      !  if (logic_tdav) then
      !    ! traditional davidson diagonalization method is used
      !    do l=1,nci_dim
      !      depff = vector2(l)-vcien(mt)
      !      !if (abs(depff) <= depc) depff = depc
      !      vector1(mtidx+l) = vector1(mtidx+l)/depff
      !    end do
      !  else
      !    ! generalized davidson diagonalization method is used
      !    call gdavdiag(mt,mtidx,vcien(mt))
      !  end if
      !end do

      mcroot = mroot-mtsta+1
      call orthog(kval,iiter,mtsta,idxvec)
      if (.not. log_muliter) then
        if (mth_eigen > 1) then
          call orthogwconvec()
        end if
      end if
      idxvec(iiter) = kval
      nd = (mroot-mtsta+1)*nci_dim
      call write_bv(lucitv1,iiter,vector1,nd)
    end if

    sc1 = c_time()
    scvp = sc1-sc0

    nd = nci_dim*mroot
    vector2(1:nd) = dzero
    ! start h*c

    call read_ml(lucidia,1,vector2,nci_dim,1)
    do mt=mtsta+1,mroot
      mtidx = indx(mt-mtsta+1)
      do l=1,nci_dim
        vector2(mtidx+l) = vector2(l)
      end do
    end do

    do mt=mtsta,mroot
      mtidx = indx(mt-mtsta+1)
      do l=1,nci_dim
        ni = mtidx+l
        vector2(ni) = vector2(ni)*vector1(ni)
      end do
    end do

    call matrix_vector_multi_parallel_drt(sechc)
    nd = nci_dim*(mroot-mtsta+1)
    call write_bv(lucitv2,iiter,vector2,nd)

  end do
  if (.not. restart) exit
end do

890 format(/,1x,'number of initial trial vectors is',i3,/,1x,'total wall clock time=',f9.2,' seconds')
!900 format(/,1x,'no.',i3,1x,'iter',/)
901 format(2x,'NITER',1x,'NROOT',3x,'TOTAL ENERGY',4x,'ENERGY DIFF',4x,'VALPHA',7x,'VRESIDE',1x,'T HC(s)',1x,'T KSPACE(s)')
902 format(2(2x,i3),2x,f16.9,1x,f12.9,1x,f12.8,1x,f12.8,1x,2(f8.2,1x))
!903 format(/,1x,'total wall time for h*c=',f8.2,1x,'seconds',/1x,'total wall time for k space calculation=',f8.2,1x,'seconds')
!910   format(/,1x,'number of converged roots is ',i4)
911 format(/,1x,'number of kspace exceeds maxium kspace dimension,',/,1x,'kspace is reseted ',/)
!...end of dav_diagonalize

end subroutine cidiagonalize

subroutine mrcibasis_rest(ndim,mroot,mjn,indx,vb1,vb2,vcien,mth_eigen,ncivec)
!***********************************************************************
! this subroutine is revised by suo bing. the initial trial
! vectors are calculated.
! on entry:
!-----------------------------------------------------------------------
!     ndim  - dimension of ci space
!     mroot - number of roots are calculated
!     mjn
!     indx  - index of mth vector in vb1 vector
!     vb1   - vector1
!     vb2   - vector2
!
!  on out:
!-----------------------------------------------------------------------
!     vb1   - trial vectors

implicit real*8(a-h,o-z)
#include "ci_parameter.fh"
#include "files_gugaci.fh"
!common /file_descript/nf1,nf2,nf3,nf4,nf7,nf8,nf9,nf10,nf11,nf13,nf15,nf20,nf21,nf22,nf23
!logical logic_tdav,logic_inivec_read
!#include "scratch.fh"
#include "program_control.fh"
data dzero/0.d0/
dimension vb1(ncivec*ndim), vb2(ncivec*ndim)
dimension indx(max_kspace), mjn(2*max_root), vcien(mroot)

!write(6,*) 'indx',indx(1:2*mroot)
!call read_bv(nf8,1,vb2,ndim)
call read_ml(lucidia,1,vb2,ndim,1)

i = mth_eigen
write(6,'(2x,2i8,f18.8)') i,mjn(i),vb2(mjn(i)),mth_eigen

! initialize vb1-vector1 and th-vector2 to zero
numdim = ndim
do m=1,numdim
  vb1(m) = dzero
  vb2(m) = dzero
end do

j = mth_eigen
if (.not. logic_inivec_read) then
  ij = indx(1)
  vb1(ij+mjn(j)) = 1.0d0
else
  call read_bv(nf23,j,vb1,ndim)
end if
!rewind(nf7)
do i=1,mth_eigen-1
  !read(nf7) vb2(1:ndim)
  call read_ml(lucivec,1,vb2,ndim,i)
  call orth_ab(ndim,vb1,vb2)
end do
call norm_a(ndim,vb1)

write(6,'(2x,2i8,f18.8)') i,mjn(i),vb1(mjn(i))

call read_ml(lucidia,1,vb2,ndim,1)
!call read_bv(nf8,1,vb2,ndim)
!vcien(1) = vb2(mjn(j))
do i=1,ndim
  vb2(i) = vb1(i)*vb2(i)
end do

!write(6,*) 'bbs debug 1'

call matrix_vector_multi_parallel_drt(sechc)
vsum = 0.d0
do i=1,ndim
  vsum = vsum+vb1(i)*vb2(i)
end do
vcien(1) = vsum

call write_bv(lucitv1,1,vb1,ndim)
call write_bv(lucitv2,1,vb2,ndim)

!write(6,*) 'bbs debug 2'

vad = 0.d0
idx = 0
outer: do m=1,ndim
  do k=1,mth_eigen
    if (m == mjn(k)) cycle outer
  end do
  if (abs(vb2(m)) > vad) then
    vad = abs(vb2(m))
    idx = m
  end if
end do outer
mjntm = idx
!write(6,*) 'mjn(kk)',idx,vb1(idx),vb2(idx)
!write(6,*) 'bbs debug 3'

numdim = ndim
vb1(1:numdim) = dzero
vb2(1:numdim) = dzero

vb2(mjntm) = 1.d0
!rewind(nf7)
do i=1,mth_eigen-1
  !read(nf7) vb1(1:ndim)
  call read_ml(lucivec,1,vb1,ndim,i)
  call orth_ab(ndim,vb2,vb1)
end do
call read_bv(lucitv1,1,vb1,ndim)
call orth_ab(ndim,vb2,vb1)
call norm_a(ndim,vb2)
vb1(1:ndim) = vb2(1:ndim)

!write(6,*) 'bbs debug 4'
call read_ml(lucidia,1,vb2,ndim,1)
!call read_bv(nf8,1,vb2,ndim)
do i=1,ndim
  vb2(i) = vb1(i)*vb2(i)
end do
!write(6,*) 'bbs debug 5'

call matrix_vector_multi_parallel_drt(sechc)
call write_bv(lucitv1,2,vb1,ndim)
call write_bv(lucitv2,2,vb2,ndim)
!write(6,*) ' initial vector 5'
!call abend()
!write(6,*) 'bbs debug 6'

return
!...end of mrcibasis_rest

end subroutine mrcibasis_rest

subroutine mrcibasis(ndim,mroot,mjn,indx,vb1,vb2,vcien,mth_eigen,ncivec)
!***********************************************************************
! this subroutine is revised by suo bing. the initial trial
! vectors are calculated.
! on entry:
!-----------------------------------------------------------------------
!     ndim  - dimension of ci space
!     mroot - number of roots are calculated
!     mjn
!     indx  - index of mth vector in vb1 vector
!     vb1   - vector1
!     vb2   - vector2
!     mth_eigen - 0 or 1
!  on out:
!-----------------------------------------------------------------------
!     vb1   - trial vectors

implicit real*8(a-h,o-z)
#include "ci_parameter.fh"
#include "files_gugaci.fh"
!common /file_descript/nf1,nf2,nf3,nf4,nf7,nf8,nf9,nf10,nf11,nf13,nf15,nf20,nf21,nf22,nf23
!logical logic_tdav,logic_inivec_read
!#include "scratch.fh"
#include "program_control.fh"
data dzero/0.d0/
dimension vb1(ncivec*ndim), vb2(ncivec*ndim), vdia(2*mroot)
dimension indx(max_kspace), mjn(2*max_root), vcien(mroot)
dimension mjntmp(mroot*2), vdiatmp(2*mroot)

call read_ml(lucidia,1,vb2,ndim,1)

!call read_bv(nf8,1,vb2,ndim)
indx(1:max_kspace) = 0
indx0 = 0
do i=1,mroot
  indx(i) = indx0
  indx0 = indx0+ndim
  vdia(i) = vb2(mjn(i))
end do

do i=1,mroot
  write(6,'(2x,2i8,f18.8)') i,mjn(i),vb2(mjn(i))
end do
mjntmp(1:mroot) = mjn(1:mroot)

! initialize vb1-vector1 and th-vector2 to zero
if (mth_eigen == 0) then
  if (logic_inivec_read) then
    numdim = ndim*mroot
    do m=1,numdim
      vb1(m) = dzero
      vb2(m) = dzero
    end do
    call read_bv(nf23,1,vb1,numdim)
  else
    numdim = ndim*mroot
    do m=1,numdim
      vb1(m) = dzero
      vb2(m) = dzero
    end do
    do j=1,mroot
      ij = indx(j)
      vb1(ij+mjn(j)) = 1.0d0
      vb2(ij+mjn(j)) = vdia(j)
      vcien(j) = vdia(j)
    end do
  end if
else
  if (logic_inivec_read) then
    vb1(1:ndim) = dzero
    vb2(1:ndim) = dzero
    call read_bv(nf23,mth_eigen,vb1(1),ndim)
  else
    vb1(1:ndim) = dzero
    vb2(1:ndim) = dzero
    j = mth_eigen
    ij = indx(1)
    vb1(ij+mjn(j)) = 1.0d0
    vb2(ij+mjn(j)) = vdia(j)
    vcien(1) = vdia(j)
    mroot = 1
  end if
end if

!write(6,*) ' initial basis vector 0',nf23

if (logic_inivec_read) then
  call read_ml(lucidia,1,vb2,ndim,1)
  !call read_bv(nf8,1,vb2,ndim)
  if (mth_eigen == 0) then
    do j=2,mroot
      idx = indx(j)
      do k=idx+1,idx+ndim
        vb2(k) = vb1(k)*vb2(k-idx)
      end do
    end do
    do k=1,ndim
      vb2(k) = vb1(k)*vb2(k)
    end do
  else
    do k=1,ndim
      vb2(k) = vb1(k)*vb2(k)
    end do
  end if
  !do i=1,100
  !  write(6,'(2(f18.9,1x))') vb1(i),vb2(i)
  !end do
  !call abend()

  call matrix_vector_multi_parallel_drt(sechc)

  if (mth_eigen == 0) then
    do j=1,mroot
      vsum = 0.d0
      idx = indx(j)
      do k=idx+1,idx+ndim
        vsum = vsum+vb1(k)*vb2(k)
      end do
      vcien(j) = vsum
      write(6,*) 'vcien(j)',vsum
    end do
  else
    vsum = 0.d0
    do k=1,ndim
      vsum = vsum+vb1(k)*vb2(k)
    end do
    vcien(1) = vsum
  end if
else
  call matrix_vector_multi_v(sechc)
end if
! write vector1 and vector2 to fort3 and fort4
call write_bv(lucitv1,1,vb1,ndim*mroot)
call write_bv(lucitv2,1,vb2,ndim*mroot)

call read_ml(lucidia,1,vb1,ndim,1)
!call read_bv(nf8,1,vb1,ndim)
if (mth_eigen == 0) then
  do kk=mroot+1,2*mroot
    ij = indx(kk-mroot)
    vad = 0.d0
    idx = 0
    outer1: do m=1,ndim
      do l=1,kk-1
        if (m == mjntmp(l)) cycle outer1
      end do
      if (abs(vb2(m+ij)) > vad) then
        vad = abs(vb2(m+ij))
        idx = m
      end if
    end do outer1
    mjntmp(kk) = idx
    vdiatmp(kk) = vb1(idx)
    write(6,*) 'mjn(kk)',idx,vb1(idx+ij),vb2(idx+ij)
  end do
else
  kk = 2
  ij = indx(kk-mroot)
  vad = 0.d0
  idx = 0
  l = mth_eigen
  outer2: do m=1,ndim
    do k=1,mth_eigen
      if (m == mjn(k)) cycle outer2
    end do
    if (abs(vb2(m)) > vad) then
      vad = abs(vb2(m))
      idx = m
    end if
  end do outer2
  mjntmp(kk) = idx
  vdiatmp(kk) = vb1(idx)
  write(6,*) 'mjn(kk)',idx,vb1(idx),vb2(idx)
end if

numdim = ndim*mroot
do m=1,numdim
  vb1(m) = dzero
  vb2(m) = dzero
end do
!write(6,*) ' initial vector 3',vdia(2),mroot

if (mth_eigen == 0) then
  do m=mroot+1,2*mroot
    im = indx(m-mroot)
    vb1(im+mjntmp(m)) = 1.d0
    vb2(im+mjntmp(m)) = vdiatmp(m)
  end do

  if (logic_inivec_read) then
    call read_bv(lucitv1,1,vb2,numdim)
    do i=mroot+1,2*mroot
      idx = indx(i-mroot)+1
      do j=1,mroot
        jdx = indx(j)+1
        call orth_ab(ndim,vb1(idx),vb2(jdx))
      end do
      do j=mroot+1,i-1
        jdx = indx(j-mroot)+1
        call orth_ab(ndim,vb1(idx),vb1(jdx))
      end do
      call norm_a(ndim,vb1(idx))
    end do
  end if
else
  vb1(mjntmp(2)) = 1.d0
  vb2(mjntmp(2)) = vdiatmp(2)

  if (logic_inivec_read) then
    call read_bv(lucitv1,mth_eigen,vb2,ndim)
    call orth_ab(ndim,vb1,vb2)
    call norm_a(ndim,vb1)
  end if
end if

if (logic_inivec_read) then
  call read_ml(lucidia,1,vb2,ndim,1)

  !call read_bv(nf8,1,vb2,ndim)
  if (mth_eigen == 0) then
    do j=2,mroot
      idx = indx(j)
      do k=idx+1,idx+ndim
        vb2(k) = vb1(k)*vb2(k-idx)
      end do
    end do
    do k=1,ndim
      vb2(k) = vb1(k)*vb2(k)
    end do
  else
    do k=1,ndim
      vb2(k) = vb1(k)*vb2(k)
    end do
  end if

  !do i=1,100
  !  write(6,'(2(f18.9,1x))') vb1(i),vb2(i)
  !end do
  !call abend()

  call matrix_vector_multi_parallel_drt(sechc)
else
  call matrix_vector_multi_parallel_drt(sechc)
end if

call write_bv(lucitv1,2,vb1,ndim*mroot)
call write_bv(lucitv2,2,vb2,ndim*mroot)

return
!...end of mrcibasis

end subroutine mrcibasis

subroutine mrcibasis_init(ndim,mroot,mjn,indx,vb1,vb2,vcien,mth_eigen,ncivec)
!***********************************************************************
! this subroutine is revised by suo bing. the initial trial
! vectors are calculated.
! on entry:
!-----------------------------------------------------------------------
!     ndim  - dimension of ci space
!     mroot - number of roots are calculated
!     mjn
!     indx  - index of mth vector in vb1 vector
!     vb1   - vector1
!     vb2   - vector2
!     mth_eigen - 0 or 1
!  on out:
!-----------------------------------------------------------------------
!     vb1   - trial vectors

implicit real*8(a-h,o-z)
#include "ci_parameter.fh"
#include "files_gugaci.fh"
#include "program_control.fh"
data dzero/0.d0/,epc/5.0d-3/
dimension vb1(ncivec*ndim), vb2(ncivec*ndim), vdia(2*mroot)
dimension indx(max_kspace), mjn(2*max_root), vcien(mroot)
real*8, allocatable :: diagelement(:)

allocate(diagelement(ndim))

call read_ml(lucidia,1,diagelement,ndim,1)

!call read_bv(nf8,1,vb2,ndim)
indx(1:max_kspace) = 0
indx0 = 0
do i=1,mroot
  indx(i) = indx0
  indx0 = indx0+ndim
  vdia(i) = diagelement(mjn(i))
end do

do i=1,mroot
  write(6,'(2x,2i8,f18.8)') i,mjn(i),diagelement(mjn(i))
end do

! initialize vb1-vector1 and th-vector2 to zero
if (mth_eigen == 0) then
  numdim = ndim*mroot
  vb1(1:numdim) = dzero
  vb2(1:numdim) = dzero
  do j=1,mroot
    ij = indx(j)
    vb1(ij+mjn(j)) = 1.0d0
    vb2(ij+mjn(j)) = vdia(j)
    vcien(j) = vdia(j)
  end do
else
  vb1(1:ndim) = dzero
  vb2(1:ndim) = dzero
  j = mth_eigen
  vb1(mjn(j)) = 1.0d0
  vb2(mjn(j)) = vdia(j)
  vcien(1) = vdia(j)
  mroot = 1
end if

!write(6,*) ' initial basis vector 0',nf23
call matrix_vector_multi_parallel_drt(sechc)
! write vector1 and vector2 to fort3 and fort4
call write_bv(lucitv1,1,vb1,ndim*mroot)
call write_bv(lucitv2,1,vb2,ndim*mroot)

! init second vector
!call read_bv(nf8,1,vb1,ndim)
if (mth_eigen == 0) then
  vb1(1:ndim*mroot) = 0.d0
  do kk=1,mroot
    ioff = indx(kk)
    ij = mjn(kk)
    vadi = diagelement(ij)
    do m=1,ij-1
      fenmu = vadi-diagelement(m)
      if (abs(fenmu) < epc) fenmu = epc
      vb1(ioff+m) = vb2(ioff+m)/fenmu
      !write(6,'(i8,2f18.8)') m,vb2(ioff+m),fenmu
    end do
    do m=ij+1,ndim
      fenmu = vadi-diagelement(m)
      if (abs(fenmu) < epc) fenmu = epc
      vb1(ioff+m) = vb2(ioff+m)/fenmu
      !write(6,'(i8,2f18.8)') m,vb2(ioff+m),fenmu
    end do
  end do

  call read_bv(lucitv1,1,vb2,ndim*mroot)
  do m=1,mroot
    ! orth with vb1
    do n=1,mroot
      call orth_ab(ndim,vb1(indx(m)+1),vb2(indx(n)+1))
    end do
    ! orth with previous vector
    do n=1,m-1
      call orth_ab(ndim,vb1(indx(m)+1),vb1(indx(n)+1))
    end do
    call norm_a(ndim,vb1(indx(m)+1))
  end do

  !write(6,*)
  !do m=1,ndim
  !  write(6,'(2x,i5,f18.8)') m,vb1(m)
  !end do
  !call abend()
else
  call abend()
end if

numdim = ndim*mroot
vb2(1:numdim) = dzero
do m=1,mroot
  kk = indx(m)
  do i=1,ndim
    vb2(i+kk) = vb1(i+kk)*diagelement(i)
  end do
end do

call matrix_vector_multi_parallel_drt(sechc)
call write_bv(lucitv1,2,vb1,ndim*mroot)
call write_bv(lucitv2,2,vb2,ndim*mroot)

deallocate(diagelement)

return
!...end of mrcibasis_init

end subroutine mrcibasis_init

subroutine get_eigvector(mtsta0,vcien,valpha,diffci,vresid,mtheigen,log_muliter)

#ifdef _XIANEST_
use control, only: toptask
#endif
#include "drt_h.fh"
#include "files_gugaci.fh"
#include "scratch.fh"
dimension valpha(max_root), vcien(max_root), vresid(max_root)
!dimension dav2(max_root),dav3(max_root),remei(max_root)
dimension dav1(max_root), vcml(max_root)
logical log_muliter

vector1 = 0.d0
nd = nci_dim*(mroot-mtsta0+1)
call read_ml(lucidia,1,vector2,nd,2)
if (log_muliter) then
  nc = 1
  do i=1,mtsta0-1
    !mtidx = indx(i)
    call read_ml(lucivec,1,vector1(nc:nc+nci_dim-1),nci_dim,i)
    nc = nc+nci_dim
  end do
  vector1(nc:nc+nd-1) = vector2(1:nd)
  eci(1:mroot) = vcien(1:mroot)

  ! write eigenvector into fort7
  nc = 1
  do i=1,mroot
    call write_ml(lucivec,1,vector1(nc:nc+nci_dim-1),nci_dim,i)
    nc = nc+nci_dim
  end do
  !write(6,*) 'xxx cof'
  !write(6,*) vector1(nci_dim+1:nci_dim+5)
  !call abend()
else
  call write_ml(lucivec,1,vector2,nci_dim,mth_eigen)
end if
!rewind(nf7)

call memcidiag_alloc()
do mt=1,mroot
  !*********************************************************************
  ! do davidson correction

  !if (.not. log_muliter) then
  call read_ml(lucivec,1,vector1,nci_dim,mt)
  !end if
  !mtidx = indx(mt)
  vcml(mt) = 0.d0
  do j=1,ndim_h0
    jr = j
    if (logic_mr) jr = irfno(j)
    vcml(mt) = vcml(mt)+vector1(jr)*vector1(jr)
  end do
  if (log_muliter) then
    de = ecih0(mt)-vcien(mt)
  else
    de = ecih0(mth_eigen)-vcien(mt)
  end if

  dedav1 = (1.d0-vcml(mt))*de
  !dedav2=dedav1/(vcml(mt))
  !dedav3=dedav1/(2*vcml(mt)-1)
  !demei=dedav2*(n_electron*(n_electron-5)+6)/(n_electron*(n_electron-1))
  !demei=dedav2*(n_electron*(n_electron-5)+6)/(n_electron*(n_electron-1))
  dav1(mt) = vcien(mt)-dedav1
  !dav2(mt)=vcien(mt)-dedav2
  !dav3(mt)=vcien(mt)-dedav3
  !remei(mt)=vcien(mt)-demei

  write(6,900)
  !write(6,901) mt,vcien(mt),dav1(mt),dav2(mt),dav3(mt),remei(mt)
  !write(nf2,901) mt,vcien(mt),dav1(mt),dav2(mt),dav3(mt),remei(mt)
  write(6,901) mt,vcien(mt),dav1(mt),vcml(mt)
  !read(nf7) vector2(1:nci_dim)
  write(6,810) mt
  !mtidx = indx(mt)
  do i=1,nci_dim
    vcof = vector1(i)
    if (abs(vcof) > cm_cri) then
      call found_a_config(i,vcof,2)
    end if
  end do
end do
call memcidiag_dealloc()

if (log_muliter) then
  write(6,800)
  write(6,900)
  do mt=1,mroot
    write(6,901) mt,vcien(mt),dav1(mt),vcml(mt)
  end do
end if

#ifdef MOLPRO
if (toptask%task(8) == 1) then
  write(6,*)
  write(6,*) '++++++++  DATA CHECK +++++++++++++++++++++++++++++'
  call checkdata('r',mroot,8,0,vcien,'MRCI','ECI')
  call checkdata('r',mroot,8,0,dav1,'MRCI','ECI_DAV')
  write(6,*) '++++++++++ END DATA CHECK ++++++++++++++++++++++++'
  write(6,*)
end if
#else
call add_info('ECI',vcien,mroot,8)
call add_info('ECI_DAV',dav1,mroot,8)
#endif

return

800 format(/,1x,'mrsdci calculation converged')
810 format(/,1x,'the main references for root ',i3,/)
900 format(/1x,'nroot',6x,'ci energy',10x,'dav energy',8x,'coef')
901 format(1x,i3,2x,2(1x,f18.9),3x,f8.6)
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(valpha)
  call Unused_real(diffci)
  call Unused_real_array(vresid)
  call Unused_integer(mtheigen)
end if
!...end of get_eigvector

end subroutine get_eigvector

subroutine matrix_vector_multi_v(sechc)

#include "drt_h.fh"
#include "scratch.fh"
external inn_ext_ss_loop_unpack
external inn_ext_st_loop_unpack
external inn_ext_ts_loop_unpack
external inn_ext_tt_loop_unpack
external inn_ext_st_drl_loop_unpack
external inn_ext_ts_drl_loop_unpack

sc1 = c_time()
log_prod = 1
call readint(1,vint_ci)
call inner_space_loop()

call readint(2,vint_ci)
call dv_drt_ci_new()
call vd_drt_ci_new()

call readint(3,vint_ci)
call sv_drt_ci_new()
call tv_drt_ci_new()
sc2 = c_time()
sechc = sc2-sc1

return

end subroutine matrix_vector_multi_v

subroutine matrix_vector_multi_d(sechc)

#include "drt_h.fh"
#include "scratch.fh"
external inn_ext_ss_loop_unpack
external inn_ext_st_loop_unpack
external inn_ext_ts_loop_unpack
external inn_ext_tt_loop_unpack
external inn_ext_st_drl_loop_unpack
external inn_ext_ts_drl_loop_unpack

log_prod = 1
sc1 = c_time()

call readint(1,vint_ci)
call inner_space_loop()

call readint(2,vint_ci)
call sd_drt_ci_new()
call td_drt_ci_new()
call ds_drt_ci_new()
call dt_drt_ci_new()
call dv_drt_ci_new()
call vd_drt_ci_new()

call readint(3,vint_ci)
call dd_drt_ci_new()

call readint(4,vint_ci)
call ext_space_loop()

sc2 = c_time()
sechc = sc2-sc1

return

end subroutine matrix_vector_multi_d

subroutine matrix_vector_multi_parallel_drt(sechc)

#include "drt_h.fh"
#include "scratch.fh"
external inn_ext_ss_loop_unpack
external inn_ext_st_loop_unpack
external inn_ext_ts_loop_unpack
external inn_ext_tt_loop_unpack
external inn_ext_st_drl_loop_unpack
external inn_ext_ts_drl_loop_unpack

log_prod = 1
sc1 = c_time()

call readint(1,vint_ci)
call inner_space_loop()

call readint(2,vint_ci)
call sd_drt_ci_new()
call td_drt_ci_new()
call ds_drt_ci_new()
call dt_drt_ci_new()
call vd_drt_ci_new()
call dv_drt_ci_new()

call readint(3,vint_ci)
call dd_drt_ci_new()
call sv_drt_ci_new()
call tv_drt_ci_new()
call ss_drt_ci_new()
call st_drt_ci_new()
call tt_drt_ci_new()
call ts_drt_ci_new()

call readint(4,vint_ci)
call ext_space_loop()

sc2 = c_time()
sechc = sc2-sc1

return

end subroutine matrix_vector_multi_parallel_drt

subroutine matrix_vector_multi_parallel_prt(sc3)

#include "drt_h.fh"
#include "scratch.fh"
external inn_ext_ss_loop_unpack
external inn_ext_st_loop_unpack
external inn_ext_ts_loop_unpack
external inn_ext_tt_loop_unpack
external inn_ext_st_drl_loop_unpack
external inn_ext_ts_drl_loop_unpack

!write(6,*) 'error stop'
!call abend()
sc1 = c_time()

mtest = 1482
log_prod = 1

call readint(1,vint_ci)
call inner_space_loop()
write(6,*) ' inner',vector2(mtest)

call readint(2,vint_ci)
call sd_drt_ci_new()
write(6,*) ' sd',vector2(mtest)
call td_drt_ci_new()
write(6,*) ' td',vector2(mtest)
call ds_drt_ci_new()
write(6,*) ' ds',vector2(mtest)
call dt_drt_ci_new()
write(6,*) ' dt',vector2(mtest)
call vd_drt_ci_new()
write(6,*) ' vd',vector2(mtest)
call dv_drt_ci_new()
write(6,*) ' dv',vector2(mtest)

call readint(3,vint_ci)
call dd_drt_ci_new()
write(6,*) ' dd',vector2(mtest)
call sv_drt_ci_new()
write(6,*) ' sv',vector2(mtest)
call tv_drt_ci_new()
write(6,*) ' tv',vector2(mtest)
call ss_drt_ci_new()
write(6,*) ' ss',vector2(mtest)
call st_drt_ci_new()
write(6,*) ' st',vector2(mtest)
call tt_drt_ci_new()
write(6,*) ' tt',vector2(mtest)
call ts_drt_ci_new()
write(6,*) ' ts',vector2(mtest)

call readint(4,vint_ci)
call ext_space_loop()
write(6,*) ' exter',vector2(mtest)

sc2 = c_time()
write(6,*) '  end this matrix_vector_multi_parallel_drt, takes',sc2-sc1,'s'
sc3 = sc2-sc1

end subroutine matrix_vector_multi_parallel_prt

subroutine orthog(kval,iiter,msta,idxvec)
!***********************************************************************
! on entry:
!   kval - dimension of current k space
! on out:
!   kval - dimension of new k space
!   iiter - iiter+1

#include "drt_h.fh"
#include "files_gugaci.fh"
dimension idxvec(max_iter)
data dzero/0.d0/

do jiter=1,iiter
  if (jiter == 1) then
    jrst = 1
  else
    jrst = idxvec(jiter-1)+1
  end if
  jren = idxvec(jiter)
  nd = (jren-jrst+1)*nci_dim
  call read_bv(lucitv1,jiter,vector2,nd)
  do mt=msta,mroot
    mtidx = indx(mt-msta+1)
    do nt=jrst,jren
      ntidx = indx(nt-jrst+1)
      vsum = dzero
      do l=1,nci_dim
        vsum = vsum+vector1(mtidx+l)*vector2(ntidx+l)
      end do
      do l=1,nci_dim
        vector1(mtidx+l) = vector1(mtidx+l)-vsum*vector2(ntidx+l)
      end do
    end do
  end do
end do

! normalization of vector msta
vsum = dzero
call norm_a(nci_dim,vector1)
kval = kval+1

do mt=msta+1,mroot
  mtidx = indx(mt-msta+1)
  do nt=msta,mt
    ntidx = indx(nt-msta+1)
    vsum = dzero
    do l=1,nci_dim
      vsum = vsum+vector1(mtidx+l)*vector1(ntidx+l)
    end do
    do l=1,nci_dim
      vector1(mtidx+l) = vector1(mtidx+l)-vsum*vector1(ntidx+l)
    end do
  end do
  call norm_a(nci_dim,vector1(mtidx+1:mtidx+nci_dim))
  kval = kval+1
end do

iiter = iiter+1

return
!...end of orthog

end subroutine orthog

subroutine orthogonalization(j,n,m,ir)

#include "drt_h.fh"
#include "files_gugaci.fh"

ir = 1
vsmax2 = 1.d10
do while (j > 1)
  vsmax1 = 0.d0
  do i=1,j-1
    vsumtmp = 0.d0
    call read_bv(lucitv1,i,vector2,n)

    do l=1,n
      vsumtmp = vsumtmp+vector1(l)*vector2(l)
    end do
    do l=1,n
      vector1(l) = vector1(l)-vsumtmp*vector2(l)
    end do
    vsmax1 = max(vsmax1,abs(vsumtmp))
  end do
  if (vsmax1 < vortho_criterion) exit
  if (vsmax1 > vsmax2) then
    ir = -1
    return
  end if
  vsmax2 = vsmax1
end do

vsumtmp = 0.d0
do l=1,n
  vsumtmp = vsumtmp+vector1(l)*vector1(l)
end do

do m=1,mth_eigen-1
  !call read_bv(nf7,m,vector2,nci_dim)
  call read_ml(lucivec,1,vector2,nci_dim,m)
  call orth_ab(nci_dim,vector1,vector2)
end do
call norm_a(nci_dim,vector1)
!vsumtmp = 1.d0/sqrt(vsumtmp)
!do l=1,n
!  vector1(l) = vector1(l)*vsumtmp
!end do

call write_bv(lucitv1,j,vector1,n)

end subroutine orthogonalization

subroutine compute_vp_matrix(j,n) !nf3

#include "drt_h.fh"
#include "files_gugaci.fh"
ij = j*(j-1)/2
vsumtmp = 0.d0
do l=1,n
  vsumtmp = vsumtmp+vector1(l)*vector2(l)
end do
vp(ij+j) = vsumtmp
!rewind(nf3)
do i=1,j-1
  call read_bv(lucitv1,i,vector1,n)

  vsumtmp = 0.d0
  do l=1,n
    vsumtmp = vsumtmp+vector1(l)*vector2(l)
  end do
  vp(ij+i) = vsumtmp
end do
call read_bv(lucitv1,i,vector1,n)

end subroutine compute_vp_matrix

subroutine compute_residual_vector_mroot(kval,mtsta,iiter,idxvec,vresid,vcien)
!***********************************************************************
! 2 mar 2007 - revised
! calculate residual vector
!      e|ck>-h|ck>
!      |ck>=\sigma(vu(i,m)*vector_i),i=1,j

#include "drt_h.fh"
#include "files_gugaci.fh"
#include "scratch.fh"
dimension idxvec(max_kspace)
dimension vresid(max_root), vcien(max_root)
data dzero/0.d0/
real*8, allocatable :: diagelement(:)

depc = 1.d-3
nd = nci_dim*(mroot-mtsta+1)
vector1(1:nd) = dzero

do jiter=1,iiter
  if (jiter == 1) then
    irts = 1
  else
    irts = idxvec(jiter-1)+1
  end if
  irte = idxvec(jiter)
  nd = (irte-irts+1)*nci_dim
  call read_bv(lucitv1,jiter,vector2,nd)
  do mt=mtsta,mroot
    mtidx = indx(mt-mtsta+1)
    do irot=irts,irte
      vuim = vu(irot,mt)
      itidx = indx(irot-irts+1)
      do l=1,nci_dim
        vector1(mtidx+l) = vector1(mtidx+l)+vuim*vector2(itidx+l)
      end do
    end do
  end do
end do

! present ci vector, second record in lucidia
nd = nci_dim*(mroot-mtsta+1)
call write_ml(lucidia,1,vector1,nd,2)

do mt=mtsta,mroot
  mtidx = indx(mt-mtsta+1)
  !venergy = vcien(mt)
  do l=1,nci_dim
    !vector1(mtidx+l) = venergy*vector1(mtidx+l)
    vector1(mtidx+l) = 0.d0
  end do
end do

do jiter=1,iiter
  if (jiter == 1) then
    irts = 1
  else
    irts = idxvec(jiter-1)+1
  end if
  irte = idxvec(jiter)
  nd = (irte-irts+1)*nci_dim
  call read_bv(lucitv2,jiter,vector2,nd)
  do mt=mtsta,mroot
    mtidx = indx(mt-mtsta+1)
    do irot=irts,irte
      itidx = indx(irot-irts+1)
      vuim = vu(irot,mt)
      do l=1,nci_dim
        vector1(mtidx+l) = vector1(mtidx+l)+vuim*vector2(itidx+l)
      end do
    end do
  end do
end do

allocate(diagelement(nci_dim))
call read_ml(lucidia,2,diagelement,nci_dim,1)
call read_ml(lucidia,2,vector2,nd,2)

! new approximation vector
ij = 0
do mt=mtsta,mroot
  mtidx = indx(mt-mtsta+1)
  vtmp = dzero
  do i=1,nci_dim
    depcc = vcien(mt)-diagelement(i)
    if (abs(depcc) < depc) depcc = depc
    vector1(mtidx+i) = (vector1(mtidx+i)-vector2(ij+i)*vcien(mt))/depcc
    vtmp = vtmp+vector1(ij+i)*vector1(ij+i)
  end do
  vresid(mt) = vtmp
  ij = ij+nci_dim
end do

deallocate(diagelement)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(kval)
!...end of compute_residule_vector

end subroutine compute_residual_vector_mroot

subroutine compute_residual_vector(kval,mtsta,iiter,idxvec,vresid,vcien)
!***********************************************************************
! 2 mar 2007 - revised
! calculate residual vector
!      e|ck>-h|ck>
!      |ck>=\sigma(vu(i,m)*vector_i),i=1,j

#include "drt_h.fh"
#include "files_gugaci.fh"
#include "scratch.fh"
dimension idxvec(max_kspace)
dimension vresid(max_root), vcien(max_root)
data dzero/0.d0/

nd = nci_dim*(mroot-mtsta+1)
vector1(1:nd) = dzero

do jiter=1,iiter
  if (jiter == 1) then
    irts = 1
  else
    irts = idxvec(jiter-1)+1
  end if
  irte = idxvec(jiter)
  nd = (irte-irts+1)*nci_dim
  call read_bv(lucitv1,jiter,vector2,nd)
  do mt=mtsta,mroot
    mtidx = indx(mt-mtsta+1)
    do irot=irts,irte
      vuim = vu(irot,mt)
      itidx = indx(irot-irts+1)
      do l=1,nci_dim
        vector1(mtidx+l) = vector1(mtidx+l)+vuim*vector2(itidx+l)
      end do
    end do
  end do
end do

nd = nci_dim*(mroot-mtsta+1)
call write_ml(lucidia,1,vector1,nd,2)

do mt=mtsta,mroot
  mtidx = indx(mt-mtsta+1)
  venergy = vcien(mt)
  do l=1,nci_dim
    vector1(mtidx+l) = venergy*vector1(mtidx+l)
  end do
end do

do jiter=1,iiter
  if (jiter == 1) then
    irts = 1
  else
    irts = idxvec(jiter-1)+1
  end if
  irte = idxvec(jiter)
  nd = (irte-irts+1)*nci_dim
  call read_bv(lucitv2,jiter,vector2,nd)
  do mt=mtsta,mroot
    mtidx = indx(mt-mtsta+1)
    do irot=irts,irte
      itidx = indx(irot-irts+1)
      vuim = vu(irot,mt)
      do l=1,nci_dim
        vector1(mtidx+l) = vector1(mtidx+l)-vuim*vector2(itidx+l)
      end do
    end do
  end do
end do

do mt=mtsta,mroot
  mtidx = indx(mt-mtsta+1)
  vtmp = dzero
  do l=1,nci_dim
    li = l+mtidx
    vtmp = vtmp+vector1(li)*vector1(li)
  end do
  vresid(mt) = vtmp
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(kval)
!...end of compute_residule_vector

end subroutine compute_residual_vector

subroutine read_ref_state(nf)

#include "drt_h.fh"
character*32 refstring

do i=1,n_ref
  do j=1,norb_dz
    iref_occ(j,i) = 2
  end do
  read(nf,*) refstring(1:norb_act)
  do j=1,norb_act
    if (refstring(j:j) == '0') iref_occ(j+norb_dz,i) = 0
    if (refstring(j:j) == '1') iref_occ(j+norb_dz,i) = 1
    if (refstring(j:j) == '2') iref_occ(j+norb_dz,i) = 2
  end do
end do

write(6,*) '    multireference mode'
write(6,*) '    the reference states are:'
do i=1,n_ref
  write(6,'(12x,60(i1))') (iref_occ(j,i),j=1,norb_inn)
end do

end subroutine read_ref_state

function maxind(m,vector)

implicit real*8(a-h,o-z)
dimension vector(m)
l = 1
am = abs(vector(l))
do j=1,m
  if (abs(vector(j)) > am) then
    l = j
    am = abs(vector(j))
  end if
end do

maxind = l

end function maxind

subroutine matrmkmul(kval,mtsta,iiter,idxvec,irset)
!***********************************************************************
! 27 feb 2007 - written by suo bing
! construct p matrix in b and h*b space
! on entry:
!   msta - start index of b space
!   mend - end index of b space
!   iiter - ith iteration

#include "drt_h.fh"
#include "files_gugaci.fh"
dimension idxvec(max_iter)
data dzero/0.d0/

if ((iiter == 2) .and. (irset == 0)) then
  nroot = mroot-mtsta+1
  ! the first iteration
  call read_bv(lucitv1,1,vector1,nci_dim*mroot)
  call read_bv(lucitv2,1,vector2,nci_dim*mroot)
  do irot=1,mroot+nroot
    if (irot == mroot+1) then
      call read_bv(lucitv1,2,vector1,nci_dim*mroot)
    end if
    ij = irot*(irot-1)/2
    if (irot <= mroot) then
      iidx = indx(irot)
    else
      iidx = indx(irot-mroot)
    end if
    jend = min(irot,mroot)
    do jrot=1,jend
      jidx = indx(jrot)
      valsum = dzero
      do i=1,nci_dim
        valsum = valsum+vector1(iidx+i)*vector2(jidx+i)
      end do
      vp(ij+jrot) = valsum
    end do
  end do
  call read_bv(lucitv2,2,vector2,nci_dim*mroot)
  do irot=mroot+1,mroot+nroot
    ij = irot*(irot-1)/2
    iidx = indx(irot-mroot)
    do jrot=mroot+1,irot
      jidx = indx(jrot-mroot)
      valsum = dzero
      do i=1,nci_dim
        valsum = valsum+vector1(iidx+i)*vector2(jidx+i)
      end do
      vp(ij+jrot) = valsum
    end do
  end do
else
  ! other iteration
  msta = idxvec(iiter-1)+1
  mend = idxvec(iiter)
  do iit=1,iiter
    if (iit == 1) then
      jrs = 1
    else
      jrs = idxvec(iit-1)+1
    end if
    jre = idxvec(iit)
    nd = (jre-jrs+1)*nci_dim
    call read_bv(lucitv2,iit,vector2,nd)
    do jrot=jrs,jre
      jidx = indx(jrot-jrs+1)
      do irot=msta,mend
        iidx = indx(irot-msta+1)
        vsum = dzero
        do l=1,nci_dim
          vsum = vsum+vector1(iidx+l)*vector2(jidx+l)
        end do
        if (irot > jrot) then
          ij = irot*(irot-1)/2
          vp(ij+jrot) = vsum
        else
          ij = jrot*(jrot-1)/2
          vp(ij+irot) = vsum
        end if
      end do
    end do
  end do
end if

!write(6,*) 'vp'
!do i=1,kval
!  nd = i*(i-1)/2
!  write(6,'(5(1x,f15.8))') vp(nd+1:nd+i)
!end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(kval)
!...end of matrmkmul

end subroutine matrmkmul

subroutine orthogwconvec()

#include "drt_h.fh"
#include "files_gugaci.fh"

!rewind(nf7)
do i=1,mth_eigen-1
  !read(nf7) vector2(1:nci_dim)
  call read_ml(lucivec,1,vector1,nci_dim,i)
  call orth_ab(nci_dim,vector1,vector2)
end do
call norm_a(nci_dim,vector1)

return
!...end of orthogwconvec

end subroutine orthogwconvec

subroutine cielement()
! print a Row of CI matrix

#include "drt_h.fh"

write(6,*) 'mcroot',mcroot
indx(1) = 0
indx(2) = nci_dim
mcroot = 1

icolum = 807
vector1(1:nci_dim) = 0.d0
vector2(1:nci_dim) = 0.d0
vector1(icolum) = 1.d0
call matrix_vector_multi_parallel_drt(sc1)

!do i=1,nci_dim
!  write(21,'(i8,1x,f18.8)') i,vector2(i)
!end do

end subroutine cielement
