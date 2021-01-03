************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine geth0()
#include "drt_h.fh"
#include "pl_structure_h.fh"
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
c     character*256 filename
      dimension vad0(max_h0)
      real*8, allocatable :: vb1(:), vb2(:)
      dimension iselcsf_occ(max_innorb,max_ref)
      data dcrita/0.5d-5/

      allocate(vb1(max_h0*max_kspace),vb2(max_h0*max_kspace))
      if(.not.logic_mr) then
        call minevalue(iselcsf_occ)
      endif
!========================================================
! calculate ndim_h0
      if(norb_act.eq.0) then
        ndim_h0=1
        vector2(1)=vector1(1)
        irfno(1)=1
        goto 100
!        return
      endif
      ipae=1
      jpae=nu_ae(ipae)
      if(jpae.eq.0) return
      jpadl=1
      if(nu_ad(jpadl).eq.0) return
      jpad=jpadl
      call seg_drt()
      if(ndim .eq. 0) return
      ndim_h0=ndim

      if(mroot.gt.ndim_h0) then
        write(6,*) '    mroot> ndim_h0, mroot,ndim_h0=',mroot,ndim_h0
        mroot=min(ndim_h0,mroot)
        write(6,*) "   ",mroot,"roots are calculated"
      endif

      call copy_to_drtl()

      if ( logic_mr ) then
        call irfrst(iselcsf_occ)
        if(mroot.gt.ndim_h0) then
          write(6,*) '    mroot> ndim_h0, mroot,ndim_h0=',mroot,ndim_h0
          mroot=min(ndim_h0,mroot)
          write(6,*) "   ",mroot,"roots are calculated"
        endif
        call  minevalue(iselcsf_occ)
      endif

!=====================================================================
100   continue
c     if(logic_mr) ndim_h0=irf
      if(.not.logic_mr) then
        ndim0=ndim_h0
      else
        ndim0=irf
      endif
      write(6,*)'     ================================'
      write(6,*)'         step 1: diagnalization h0   '
      write(6,*)'            ndim_h0=', ndim0
      write(6,*)'     ================================'

      if(ndim_h0.eq.1) then
        ecih0(1)=escf(1)
        vcm(1)=1.d0
        goto 400
      endif

      call formh0()   ! for log_mr, ndim_h0 changed to irf in this subro
!=====================================================================
      if (associated(vcm)) deallocate(vcm)
      allocate(vcm(ndim_h0*mroot))

      if(ndim_h0.le.30) then
        call hotred(max_kspace,ndim_h0,vector2,vd,ve,vu)
        call qlcm(max_kspace,ndim_h0,vd,ve,vu)
        ijm=0
        do m=1,mroot
          ecih0(m)=vd(m)
          do l=1,ndim_h0
            vcm(ijm+l)=vu(m,l)
          enddo
          ijm=ijm+ndim_h0
        enddo
        goto 400
      else
        mmspace=mroot*3+10
        do i=1,mmspace
          indx(i)=(i-1)*ndim_h0
        enddo

        do i=1,ndim_h0
          mn=i*(i+1)/2
          vad0(i)=vector2(mn)
        enddo

        nxh=ndim_h0*(ndim_h0+1)/2
        nxb=ndim_h0*max_kspace

        call basis_2(ndim_h0,vb1,nxb,vad0,vector2,nxh)

        do m=1,mroot
          ecih0(m)=escf(m)
        enddo

        kval=mroot*2
        call hymat_2(max_root,max_kspace,ndim_h0,kval,mroot,dcrita,
     *               ecih0,vcm,indx,vector2,nxh,vb1,vb2,nxb,vad0)
c        vcm(1:mroot*ndim_h0)=vb1(1:mroot*ndim_h0)
c save ci vector in h0 into vb2
c        numh0=nci_h0 !iw_sta(2,1)
c        vb2(1:numh0*mroot)=0.d0
c        if(logic_mr) then
c          idx1=0
c          idx2=0
c          do i=1,mroot
c            do j=1,ndim_h0
c              m=irfno(j)
c              vb2(idx1+m)=vb1(idx2+j)
c            enddo
c            idx1=idx1+numh0
c            idx2=idx2+ndim_h0
c          enddo
c        else
c          vb2(1:numh0*mroot)=vb1(1:numh0*mroot)
c        endif
      endif
      deallocate(vcm)

400   continue
      write(6,*)
      do m=1,mroot
        write(6,'(5x,a7,i5,f18.8)')   ' root,',m,ecih0(m)
      enddo
      write(6,*)

c      filename=tmpdir(1:len_str)//"/fort7"
c      len=len_str+6
c      open(nf7,file=filename(1:len),form="unformatted")
c      write(nf7) vb2(1:numh0*mroot)
c      close(nf7)

c      open(100,file="tmp.dat")
c      do i=1,ndim_h0
c        write(100,"(1x,f18.9,2i8)") vb1(i),i,irfno(i)
c      enddo
c      write(100,*) "ndim_h0=",ndim_h0,"v0=",numh0
c      idx1=0
c      idx2=0
c      do i=1,mroot
c        write(100,*) "mroot=",mroot
c        do j=1,numh0
c          write(100,"(1x,i2,1x,i5,2x,f18.9)") 1,j,vb2(j)
c          m=ifrno(j)
c          write(100,"(2(1x,i8,1x,f18.9))") j,vb2(j+idx1),m,vb1(m+idx2)
c        enddo
c        idx1=idx1+numh0
c        idx2=idx2+ndim_h0
c      enddo
c      close(100)
c      stop 888
      return
      end

      subroutine formh0()
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
      common /scratch/ tmpdir,len_str
      character*256 tmpdir
      real*8, allocatable :: buff(:)

      num=ndim_h0*(ndim_h0+1)/2
      if(num.gt.max_vector) then
        write(6,*) " no enough space to store h0 matrix",num
#ifdef MOLPRO
#else
      call abend()
#endif
#ifdef _XIANEST_
#endif
!       call abend
!        stop 888
      endif

      vector2(1:lenvec)=0.d0
      log_prod=2
      call readint(1,vint_ci)
! act complete loop
      call cloop_in_act()
! dbl- act loop
      call ploop_in_act()

c save the h0 matrix into sracth file 23 to use later
      if ( logic_mr ) then ! rst
        allocate(buff(ndim_h0))
        buff(1:ndim_h0)=vector1(1:ndim_h0)
        vector1(1:lenvec)=0.d0
        mnh0=ndim_h0*(ndim_h0+1)/2
        vector1(1:mnh0)=vector2(1:mnh0)
        mnrh0=irf*(irf+1)/2
        vector2(1:mnrh0)=0.d0
        do ir1=1,irf
          iconf1=irfno(ir1)
          ii=ir1*(ir1+1)/2
          vector2(ii)=buff(iconf1)
          do ir2=1,ir1-1
            iconf2=irfno(ir2)
            ii=ir1*(ir1-1)/2+ir2
            iconfmax=max(iconf1,iconf2)
            iconfmin=min(iconf1,iconf2)
            iicc=iconfmax*(iconfmax-1)/2+iconfmin
            vector2(ii)=vector1(iicc)
          enddo
        enddo
        ndim_h0=irf
        deallocate(buff)
      else
        do i=1,ndim_h0   ! rcas
          ii=i*(i-1)/2+i
          vector2(ii)=vector2(ii)+vector1(i)
        enddo
      endif

      !do i=1,ndim_h0
      !  write(6,"(i4,1x,f12.6)") i,vector1(i)
      !enddo
      !stop 888
      write(6,"(a24,i5)") " dimension of h0 space= ",ndim_h0
      log_prod=1

c      open(100,file="h0_new")
c      do i=1,num
c        write(100,*) i,vector2(i)
c      enddo
c      close(100)
c      stop 777
      return
      end

c mroot=2 min_space=4
c    b1=(0,0,1.0,0,0,0...) b2=(0,0,0,0,1.0,0...)
      subroutine basis_2(ndim,vb1,nxb,vad,th,nxh)
#include "drt_h.fh"
      data dzero/0.d0/dcrita/1.0d-6/epc/5.0d-3/
      dimension vb1(max_kspace*ndim),vad(ndim)
      dimension th(nxh),ijb1(mroot)
      vb1=0.d0

      do j=1,mroot
        ij=indx(j)
        mief=mjn(j)
        if(logic_mr)then
          mjnj=mjn(j)
          mief=ifrno(mjnj)
        endif
        do 60 l=1,ndim
        vb1(ij+l)=dzero
60      continue
        vb1(ij+mief)=1.0d0
      enddo

c=================================================================
      j=mroot
      do m=1,mroot
        i=mjn(m)
        if(logic_mr)then
          mjnj=mjn(m)
          i=ifrno(mjnj)
        endif
        vadi=vad(i)
        ijh=i*(i-1)/2
        j=j+1
        ijb1(m)=indx(j)
        do l=1,i-1
          fenmu=vadi-vad(l)
          if(abs(fenmu).lt.epc) fenmu=epc
          vb1(ijb1(m)+l)=th(ijh+l)/fenmu
        enddo
        do l=i+1,ndim
          fenmu=vadi-vad(l)
          if(abs(fenmu).lt.epc) fenmu=epc
          ijh=l*(l-1)/2
          vb1(ijb1(m)+l)=th(ijh+i)/fenmu
        enddo
      enddo

c--------------------------------------------------------------------
c write out basis
c--------------------------------------------------------------------
c500    write(nf2,*)' l    vb5       vb6      vb7       vb8     '
c
c       do l=1,ndim
c        write(nf2,'(2x,i5,4f10.4)')l,vb1(indx(5)+l),vb1(indx(6)+l)
c     :                             ,vb1(indx(7)+l),vb1(indx(8)+l)
c     enddo
c      stop 777   !wyb_tmp
c--------------------------------------------------------------------
      do m0=1,mroot
        ib=m0+mroot
        call orthnor(ndim,ib,dcrita,vb1,nxb)
      enddo

      return
      end


      subroutine minevalue(iselcsf_occ)
#include "drt_h.fh"
#include "files_gugaci.fh"
      common/config/ndr,nwalk(0:max_orb)
      dimension iselcsf_occ(max_innorb,max_ref)
      dimension iwalktmp(0:max_orb)
      data dzero/0.d0/dcrita/1.0d-6/epc/5.0d-3/

      call read_ml(lucidia,1,vector1,nci_dim,1)

      vector2(1:nci_dim)=vector1(1:nci_dim)
      ndimh0=nci_h0 !iw_sta(2,1)

      if(.not.logic_mr) then
        do i=1,mroot
          l=1
          am=vector2(l)
          do j=1,ndimh0
            if(vector2(j).ne.dzero.and.vector2(j).lt.am) then
              l=j
              am=vector2(j)
            endif
          enddo
          mjn(i)=l
          vector2(l)=dzero
        enddo
      else
        do i=1,mroot
          l=1
          am=vector2(l)
          do j=1,irf
            jm=irfno(j)
            if(vector2(jm).ne.dzero.and.vector2(jm).lt.am) then
              l=jm
              am=vector2(jm)
            endif
          enddo
          mjn(i)=l
          vector2(l)=dzero
        enddo
      endif

      do m=1,mroot
        escf(m)=vector1(mjn(m))
      enddo

      write(6,*) '   mjn(k) :',mroot
      do m=1,mroot
        call found_a_config(mjn(m),escf(m),1)
        do i=1,norb_all
          iwalktmp(i)=nwalk(norb_all-i+1)
        enddo
        ij=norb_dz
        if(m.le.2*mroot) then
          do io=1,norb_act
            ij=ij+1
            if(iwalktmp(ij).eq.3) iselcsf_occ(io,m)=3
            if(iwalktmp(ij).eq.2) iselcsf_occ(io,m)=2
            if(iwalktmp(ij).eq.1) iselcsf_occ(io,m)=1
            if(iwalktmp(ij).eq.0) iselcsf_occ(io,m)=0
          enddo
c          write(6,"(16(i1))") iwalktmp(norb_dz+1:norb_dz+norb_act)
        endif
!        write(6 ,'(2x,2i8,f18.8)') m,mjn(m),escf(m)
!       write(nf2,'(2x,2i8,f18.8)') m,mjn(m),escf(m)
      enddo
      write(6,*)
      end

      subroutine orthnor_ab(n,av,bv,id)  !bv:basis, av:vector for orth a
      real*8 av(n),bv(n),s,ddot_,dcrita
      dcrita=1.0e-10
      if(id.ne.0) goto 150
c     orthogonization av,bv
      s=ddot_(n,av,1,bv,1)
      do i=1,n
        av(i)=av(i)-s*bv(i)
      enddo
c     normalization of av_eigenvector.
150   s=0.0d0
      s=ddot_(n,av,1,av,1)
      s=sqrt(s)
      s=max(s,dcrita)
      do i=1,n
        av(i)=av(i)/s
      enddo

      return
      end

      real*8 function ddot_bak(n,dx,dy)
      real*8 dx(n),dy(n),s
      s=0.0d0
      do l=1,n
        s=s+dx(l)*dy(l)
      enddo
      ddot_bak=s
      return
      end

      subroutine matrmk_1(k)
#include "drt_h.fh"
#include "files_gugaci.fh"
      do ibas=1,k
        call read_bv(lucitv1,ibas,vector1,nci_dim)
        ij=ibas*(ibas-1)/2
        do jbas=1,ibas
          call read_bv(lucitv2,jbas,vector2,nci_dim)
          vsumtmp=0.d0
          do l=1,nci_dim
            vsumtmp=vsumtmp+vector1(l)*vector2(l)
          enddo
          vp(ij+jbas)=vsumtmp
        enddo
      enddo
      write(6,*)
      il=0
      do l=1,k
        write(6,1112) (vp(i),i=il+1,il+l)
        il=il+l
      enddo
      write(6,*)
1112  format(2x,20f14.8)
      return
      end

      subroutine matrix_vector_mul_h0()


      end
