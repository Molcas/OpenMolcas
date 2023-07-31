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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine klvaa_oovo(x,g,vblock,N,nug,LU,last,ias)
!mp !subroutine klvaa_oovo(G,ix,it,ig,iscr,vblock,N,nug,LU,last,ias)
!mpn subroutine klvaa_oovo(ix,it,ig,iscr,vblock,N,nug,LU,last,ias)
!
! creates L(alpha>alpha,alpha-alpha)
! DA files LMATICH(ISP)ICH(ISP)
! max G at this place
!
! parallelization (seems to be) irrelevant at the moment
! implemented integer offsets, PV, 14 may 2004.

use ChT3_global, only: DimGrpaR, maxdim, nblock, nc, NNOAB, NNUAB, no, NOAB, NUAB, printkey
use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), allocatable, intent(_OUT_) :: x(:)
real(kind=wp), allocatable, intent(inout) :: g(:)
integer(kind=iwp), intent(in) :: vblock, N, nug, LU
integer(kind=iwp), intent(out) :: last, ias
integer(kind=iwp) :: A, A1, A2, a_tmp, AADT, ADIM, B1, b1_tmp, B2, b2_chk, i, i_blk, iasblock, IJS, IS2, isp, j, j_blk, j_tmp, k, &
                     KADT, KI, length, length1, length2, MAXDIM2, NGA, ngaa, ngaf, ngal, NGB, ngbf, ngbl, nind_ngaf, nind_ngal, &
                     nind_ngbf, nind_ngbl, nno, NSTEP, R, RAD
real(kind=wp), allocatable :: l1(:), l2(:), t2_exp(:), t2_tmp(:), tmp(:), tmp2(:), x2(:)

!mp write(u6,*) 'Entering klvaa_oovo'

ISP = 1
iasblock = NNOAB(ISP)*vblock*N/nblock
if ((iasblock*nblock) < (NNOAB(ISP)*vblock*N)) iasblock = iasblock+1
nno = nTri_Elem(noab(isp))
!!call w_alloc(ix,nnoab(ISP)*vblock*n,'IX klvaao')
!!call w_alloc(ig,noab(isp)*nuab(isp)*nnoab(isp),'IG klvaao')
!!call w_alloc(iscr,nuab(isp)*nuab(isp)*nnoab(isp),'ISCR klvaao')
!FN = 'OOVOI'//ICH(ISP)
!!call GET3DM(FN,G(ig),NUAB(ISP)*NOAB(ISP),NNOAB(ISP),0)

!call EXPA1_UHF(G(IT),nnoab(isp),NUAB(ISP),-1,G(ISCR))
! expa done here. remains in it address
!mpn IJS = IT
!mpn do i=2,noab(isp)
!mpn   do j=1,i-1
!mpn     KADT = IT+(i-1)*noab(isp)*nnuab(3)+(j-1)*nnuab(3)
!mpn     !!call dcopy_(NNUAB(3),G(KADT),1,G(IJS),1)
!mpncmp  !call vneg_cht3(G(KADT),1,G(IJS),1,NNUAB(3))
!mpn     Work(IJS:IJS+NNUAB(3)-1) = -Work(KADT:KADT+NNUAB(3)-1)
!mpncmp  !call TRANSM_A(G(KADT),G(IJS),NUAB(ISP),NUAB(ISP))
!mpn     call TRANSM_A(Work(KADT),Work(IJS),NUAB(ISP),NUAB(ISP))
!mpn     IJS = IJS+NNUAB(3)
!mpn   end do ! j
!mpn end do ! i

IAS = 1
ngaa = 0
!mpn
do A1=1,NUAB(ISP),vblock
  ! not needed  call zeroma(g(ix),1,nnoab(isp)*vblock*n)
  adim = min(vblock,nuab(isp)-A1+1)
  A2 = A1+adim-1
  !mpn
  ngaa = ngaa+1

  !mp write(u6,*)
  !mp write(u6,*) '================================='
  !mp write(u6,*) ' nga ',ngaa
  !mp write(u6,*) '================================='
  !mp write(u6,*)

  !mp write(u6,'(A,2(i5,2x))') 'b1,b2 = ',a1,a2

  call block_interf(1,1,a1,a2,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

  !mp write(u6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',ngbf,ngbl,nind_ngbf,nind_ngbl

  ! - read amplitudes T2(nv,vblock,j<i)

  ! - calculate memory requirements

  length1 = nuab(isp)

  length2 = 0
  do i_blk=ngbf,ngbl
    length2 = length2+DimGrpaR(i_blk)
  end do

  !mp write(u6,*) 'length1, NUAB   = ',length1,NUAB(1)
  !mp write(u6,*) 'length2, vblock = ',length2,vblock

  length = length1*length2*nTri_Elem(no-1)
  !mp write(u6,*) 'length for blocked T2 amplitudes = ',length

  ! - setup memory

  !mp write(u6,*) 'allocating t2_exp = ',length
  call mma_allocate(t2_exp,length,label='t2_length')

  ! - read pertinent files and store them in the new blocked structure

  call mma_allocate(t2_tmp,maxdim*maxdim*no*no,label='t2_tmp')
  call mma_allocate(tmp2,maxdim*maxdim*no*no,label='tmp2')

  call gather_t2anti_blocked(length1,length2,ngbf,ngbl,t2_exp,t2_tmp,tmp2)

  call mma_deallocate(tmp2)
  call mma_deallocate(t2_tmp)

  !mpn
  ! copies oovo
  ! oovo (k>l,ai)
  K = 0
  do I=2,NOAB(ISP)
    do J=1,I-1
      K = K+1
      !!do K = 1,NNOAB(ISP)
      IJS = (K-1)*adim*N
      !!KADT = (a1-1)*noab(isp)+(K-1)*noab(isp)*nuab(isp)+1
      do a=A1,A2
        KADT = (J-1)*nno*nuab(isp)+(A-1)*nno
        do r=1,noab(isp)
          ! <ra|ij>
          !mp G(IJS+r) = G(KADT+iTri(r,i))
          x(IJS+r) = g(KADT+iTri(r,i))
        end do
        KADT = (I-1)*nno*nuab(isp)+(A-1)*nno
        do r=1,noab(isp)
          ! -<ra|ji>
          !mp G(IJS+r) = G(IJS+r)-G(KADT+iTri(r,j))
          x(IJS+r) = x(IJS+r)-g(KADT+iTri(r,j))
        end do
        !!call dcopy_(noab(isp),G(KADT),1,G(IJS),1)
        IJS = IJS+N
        !!KADT = KADT+noab(isp)
      end do
      ! now the T2
      !mpn
      !!KADT = ISCR+(K-1)*NUAB(ISP)*NUAB(ISP)+(A1-1)*NUAB(ISP)
      !mp KADT = IT+(K-1)*NUAB(ISP)*NUAB(ISP)+(A1-1)*NUAB(ISP)
      IJS = 1+NOAB(ISP)+(K-1)*adim*N

      a_tmp = a1-nind_ngbf
      !mp write(u6,*) 'K, a1, a_tmp',k,a1,a_tmp

      KADT = 1+(K-1)*NUAB(ISP)*length2+(a_tmp-1)*NUAB(ISP)

      do a=A1,A2
        !mp call dcopy_(NUAB(isp),G(KADT),1,G(IJS),1)
        x(IJS:IJS+NUAB(isp)-1) = t2_exp(KADT:KADT+NUAB(isp)-1)
        !mp write(u6,*) (t2_exp(KADT+a_tmp),a_tmp=0,NUAB(isp)-1)
        !!write(u6,'(A,2I3,11D10.4)') 'OT',K,a,(G(r),r=IJS-noab(isp),IJS+nuab(isp)-1)

        KADT = KADT+NUAB(ISP)
        IJS = IJS+N
      end do     !A
    end do       !J
  end do         !I
  !!end do         !K
  !!write(u6,'(A,2I5,4x,5D15.10)') 'block-m:a1,IAS,ddot',a1,ias,ddot_(N*adim*nnoab(ISP),G(IX),1,G(IX),1),(G(I),I=IX,IX+3)
  !mp call multi_wridir(G(IX),n*adim*nnoab(isp),LU,IAS,last)
  !mp
  !mpn do jjj=1,n*adim*nnoab(isp)
  !mpn   if (abs(x(jjj)) > 1.0e5_wp) then
  !mpn     write(u6,*) 'fucko 2'
  !mpn     write(u6,*) jjj,x(jjj)
  !mpn     stop
  !mpn   end if
  !mpn end do
  !mp
  call multi_wridir(x,N*adim*nnoab(isp),LU,IAS,last)
  IAS = IAS+iasblock

  !mp
  call mma_deallocate(t2_exp)
  !mp
end do     ! A1
  !mp write(u6,*) 'cast 1 ok'
  !mpn

IS2 = 3-ISP
! lmat
!mp ! Mozes odjbt ix, iscr, ig
!mp write(u6,*) 'test 2 na iscr ',vblock*noab(1)*noab(1)
!mp write(u6,*) 'test 2 na ig ',noab(1)*nuab(1)*nno
!mp write(u6,*) 'test 2 na ix ',noab(1)*noab(1)*vblock*n
!mp write(u6,*) 'nno (2) = ',nno
!call GetMem('c2_iscr','Free','Real',iscr,vblock*noab(1)*noab(1))
call mma_deallocate(g)
call mma_deallocate(x)
!mp write(u6,*) 'n (2) = ',n
!mp call w_memchk('IX klvaa ')
!mp call w_free(g(ix),0,'klvaaix')
!mp call w_alloc(ix,nnoab(isp)*vblock*vblock,'Ix klvaa-v')
!mp
call mma_allocate(x2,nnoab(isp)*vblock*vblock,label='klv_oo_ix')
!mp
! starts <aa||oo> integrals
!FN = 'VVOOI'//ich(3)
!!FN = 'VVOOI'//ich(isp)
!mp !call GET3DM(FN,G(it),NNUAB(3),NNOAB(3),0)
!mp
!mp !call w_alloc(il1,nc*no*maxdim,'IL1 klvaa-v')
!mpn call mma_allocate(l1,nc*no*maxdim,label='klv_oo_il1')
!mp !call w_alloc(itmp,max(nc*no*maxdim,maxdim*maxdim*no*no),'IL2 klvaa-v')
!mpn call mma_allocate(tmp,max(nc*no*maxdim,maxdim*maxdim*no*no),label='klv_oo_itmp')
!mp !call w_alloc(il2,max(nc*no*maxdim,maxdim*maxdim*no*no),'IL2 klvaa-v')
!mpn call mma_allocate(l2,max(nc*no*maxdim,maxdim*maxdim*no*no),label='klv_oo_il2')

!mp !call gen_vvoo(G(it),G(il1),G(itmp),G(il2))
!mpn call gen_vvoo(Work(it),l1,tmp,l2)

!mp !open(unit=36,file='vvoo_moje')
!mp !do i=0,NNUAB(3)*NNOAB(3)-1
!mp !  if (abs(G(it+i)) < 1.0e-7_wp) G(it+i) = Zero
!mp !  write(36,*) i,G(it+i)
!mp !end do
!mp !close (36)

!mp !call w_free(G(il1),0,'IL1 klvaa-v')
!mpn call mma_deallocate(l2)
!mpn call mma_deallocate(tmp)
!mpn call mma_deallocate(l1)

!mp write(u6,*) ddot_(nnoab(3)*nnuab(3),G(it),1,G(it),1)

!mp
!!call dscal_(NNUAB(3)*NNOAB(3),-One,G(it),1)
!!call dscal_(NNUAB(isp)*NNOAB(ISP),-One,G(it),1)

! number of blocks written in a single multiwrite

iasblock = vblock*vblock*nnoab(isp)/nblock
if (iasblock*nblock < vblock*vblock*nnoab(isp)) iasblock = iasblock+1
!!write(u6,*) 'create_aa vvoo   iasblock',iasblock

!!FN = 'VMAT'//ich(isp)//ich(isp)
!!call multi_opendir(FN,LU)
! currently using 3-dim (big field) - will be replaced after changing
! stepiv and the rest
do nga=1,nug
  A1 = (nga-1)*vblock+1
  adim = min(vblock,nuab(isp)-A1+1)
  A2 = A1+adim-1
  do ngb=1,nga
    if (nga == ngb) then
      maxdim2 = nTri_Elem(adim-1)
    else
      maxdim2 = adim*vblock
    end if
    B1 = (ngb-1)*vblock+1
    !mpn
    !mp write(u6,*)
    !mp write(u6,*) '================================='
    !mp write(u6,*) ' nga, ngb',nga,ngb
    !mp write(u6,*) '================================='
    !mp write(u6,*)

    !mpn
    ! - check the largest b2

    !mp b2_chk = b1-1+min(vblock,a2-b1)
    !mpn b2_chk = b1-1+vblock
    b2_chk = b1+min(vblock,nuab(isp)-b1+1)-1

    ! - find out which T2 blocked files will be needed for particular nga, ngb

    !mp write(u6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk

    call block_interf(a1,a2,b1,b2_chk,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

    !mp write(u6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',ngaf,ngal,nind_ngaf,nind_ngal
    !mp write(u6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',ngbf,ngbl,nind_ngbf,nind_ngbl

    ! - read amplitudes from T2_ngaf_ngbf ...  T2_ngaf_ngbl, nga>=ngb
    !                           ....               ....
    !                        T2ngal_ngbf  ...  T2_ngal_ngbl, nga>=ngb

    ! - calculate memory requirements (consider squared T2(a',a'))

    length1 = 0
    do i_blk=ngaf,ngal
      length1 = length1+DimGrpaR(i_blk)
    end do

    length2 = 0
    do j_blk=ngbf,ngbl
      length2 = length2+DimGrpaR(j_blk)
    end do

    !mp write(u6,*) 'length1, vblock = ',length1,vblock
    !mp write(u6,*) 'length2, vblock = ',length2,vblock

    length = length1*length2*no*no
    !mp write(u6,*) 'length for blocked VVOO integrals = ',length

    ! - setup memory

    !mp write(u6,*) 'allocating t2_exp = ',length
    call mma_allocate(t2_exp,length,label='t2_exp')

    ! - read pertinent files and generate block of vvoo integrals

    call mma_allocate(l1,nc*no*maxdim,label='vvooil1')
    call mma_allocate(tmp,max(nc*no*maxdim,maxdim*maxdim*no*no),label='vvooitmp')
    call mma_allocate(l2,max(nc*no*maxdim,maxdim*maxdim*no*no),label='vvooil2')
    !mp
    call gen_vvoo_blocked(t2_exp,l1,tmp,l2,length1,length2,ngaf,ngal,ngbf,ngbl)
    !mp
    call mma_deallocate(l2)
    call mma_deallocate(tmp)
    call mma_deallocate(l1)
    !mpn
    do a=a1,a2
      B2 = B1-1+min(vblock,A-B1)
      NSTEP = B2-B1+1
      if (nstep /= 0) then
        if (nga == ngb) then
          IJS = nTri_Elem(a-a1-1)+1
        else
          IJS = (a-a1)*vblock+1
        end if
        R = 0
        do I=2,NOAB(ISP)
          do J=1,I-1
            !!R = R+1
            !!do R=1,NNOAB(ISP)
            !!  KADT = (R-1)*NNUAB(ISP)
            !mpn R = (J-1)*noab(isp)+I
            !mpn KADT = (R-1)*NNUAB(3)
            !mpn KADT = KADT+(a-1)*NUAB(ISP)+B1+IT-1

            a_tmp = a-nind_ngaf
            b1_tmp = b1-nind_ngbf

            !mp write(u6,'(A,5(i4,x))') 'a,b1,I,J,nstep ',a,b1,I,J,nstep
            !mp write(u6,'(A,5(i4,x))') 'a_tmp,b1_tmp   ',a_tmp,b1_tmp
            do j_tmp=0,nstep-1

              KADT = (i-1)*noab(isp)*length1*length2+(j-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+a_tmp

              x2(IJS+j_tmp) = t2_exp(KADT)

            end do

            ! address and block for the A1-A2x B1-B2
            !!KADT = KADT+nTri_Elem(a-2)+B1+IT-1
            ! VO >>> G(IX)
            !mpn x2(IJS:IJS+NSTEP-1) = t2_exp(KADT:KADT+NSTEP-1)
            !mpn R = (I-1)*noab(isp)+J
            !mpn KADT = (R-1)*NNUAB(3)
            !mpn KADT = KADT+(a-1)*NUAB(ISP)+B1+IT-1
            do j_tmp=0,nstep-1

              KADT = (j-1)*noab(isp)*length1*length2+(i-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+a_tmp

              x2(IJS+j_tmp) = x2(IJS+j_tmp)-t2_exp(KADT)

            end do
            !mpn
            !mpn x2(IJS:IJS+NSTEP-1) = x2(IJS:IJS+NSTEP-1)-t2_exp(KADT:KADT+NSTEP-1)
            !!write(u6,'(A,2I3,8D15.8)')'T',nTri_Elem(I-2)+j,a,(G(r),r=IJS,IJS+NSTEP-1)
            IJS = IJS+maxdim2
          end do       ! RI
        end do         ! RJ
      end if
    end do             ! a1
    !!write(u6,'(A,4I5,4x,5D15.10)') 'block-m: a1,b1,IAS,ddot',a1,b1,ias,maxdim2,ddot_(nnoab(isp)*maxdim2,G(IX),1,G(IX),1), &
    !!                               (G(I),I=IX,IX+3)
    if (maxdim2 == 0) then
      maxdim2 = 1
    end if
    !mp call multi_wridir(G(IX),nnoab(isp)*maxdim2,LU,IAS,last)
    !mp
    !mpn do jjj=1,nnoab(isp)*maxdim2
    !mpn   if (abs(x2(jjj)) > 1.0e5_wp) then
    !mpn     write(u6,*) 'fucko 3'
    !mpn     write(u6,*) jjj,x2(jjj)
    !mpn     stop
    !mpn   end if
    !mpn end do
    !mp
    call multi_wridir(x2,nnoab(isp)*maxdim2,LU,IAS,last)
    ias = ias+iasblock

    !mp write(u6,*) 'deallocating t2_exp = ',length
    call mma_deallocate(t2_exp)

  end do    ! ngb
end do      ! nga
!mp write(u6,*) 'cast 2 ok'
!mp call w_memchk('all klvaa ')
!mp call w_free(g(ix),0,'klvaa ')
!mp
call mma_deallocate(x2)
!mp
!!call w_free(g(it),0,'klvaa ')
!mp call w_alloc(ix,nnoab(3)*vblock*vblock,'ix-vvoo')
!mp
call mma_allocate(x2,nnoab(3)*vblock*vblock,label='klv_oo_ix')
!mp
!!call w_alloc(ig,nnoab(3)*nnuab(3),'ig-vvoo')
!mpn ig = it
! from now on as for uhf
iasblock = nnoab(3)*vblock*vblock/nblock
if (iasblock*nblock < nnoab(3)*vblock**2) iasblock = iasblock+1
!!FN = 'VVOOI'//ICH(3)
!!CALL GET3DM(FN,G(IG),NNUAB(3),NNOAB(3),0)
! (c>d|AK)

!mpn
nga = 0
!mpn
do A1=1,NUAB(ISP),vblock
  A2 = A1+min(vblock,nuab(isp)-A1+1)-1
  adim = A2-A1+1
  !mpn
  nga = nga+1
  ngb = 0
  !mpn
  do B1=1,NUAB(IS2),vblock
    !mpn
    ngb = ngb+1
    !mpn
    NSTEP = min(vblock,nuab(is2)-B1+1)
    !mpn
    !mpn
    !mp write(u6,*)
    !mp write(u6,*) '================================='
    !mp write(u6,*) ' nga, ngb',nga,ngb
    !mp write(u6,*) '================================='
    !mp write(u6,*)

    ! - check the largest b2

    !mp b2_chk = b1-1+min(vblock,nuab(is2)-B1+1)
    !mpn b2_chk = b1-1+vblock
    b2_chk = b1+min(vblock,nuab(is2)-B1+1)-1

    ! - find out which T2 blocked files will be needed for particular nga, ngb

    !mp write(u6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk

    if (nga < ngb) then
      !mp write(u6,*) 'switching nga, ngb',ngb,nga

      call block_interf(b1,b2_chk,a1,a2,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

    else

      call block_interf(a1,a2,b1,b2_chk,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

    end if

    !mp write(u6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',ngaf,ngal,nind_ngaf,nind_ngal
    !mp write(u6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',ngbf,ngbl,nind_ngbf,nind_ngbl

    ! - read amplitudes from T2_ngaf_ngbf ...  T2_ngaf_ngbl, nga>=ngb
    !                           ....               ....
    !                        T2ngal_ngbf  ...  T2_ngal_ngbl, nga>=ngb

    ! - calculate memory requirements (consider squared T2(a',a'))

    length1 = 0
    do i_blk=ngaf,ngal
      length1 = length1+DimGrpaR(i_blk)
    end do

    length2 = 0
    do j_blk=ngbf,ngbl
      length2 = length2+DimGrpaR(j_blk)
    end do

    !mp write(u6,*) 'length1, vblock = ',length1,vblock
    !mp write(u6,*) 'length2, vblock = ',length2,vblock

    length = length1*length2*no*no
    !mp write(u6,*) 'length for blocked VVOO integrals = ',length

    ! - setup memory

    !mp write(u6,*) 'allocating t2_exp = ',length
    call mma_allocate(t2_exp,length,label='t2_exp')

    ! - read pertinent files and generate block of vvoo integrals

    call mma_allocate(l1,nc*no*maxdim,label='vvooil1')
    call mma_allocate(tmp,max(nc*no*maxdim,maxdim*maxdim*no*no),label='tmp')
    call mma_allocate(l2,max(nc*no*maxdim,maxdim*maxdim*no*no),label='vvooil2')
    !mp
    call gen_vvoo_blocked(t2_exp,l1,tmp,l2,length1,length2,ngaf,ngal,ngbf,ngbl)
    !mp
    call mma_deallocate(l2)
    call mma_deallocate(tmp)
    call mma_deallocate(l1)
    !mpn
    !!bdim = NSTEP
    ! mv T2(B,A,I,K) >> G(ix)
    KI = 0
    do K=1,noab(isp)
      KADT = (K-1)*NNUAB(3)*(NOAB(2)*(2-ISP)+ISP-1)
      do I=1,NOAB(IS2)
        KI = KI+1
        !mpn IADT = (I-1)*NNUAB(3)*(NOAB(2)*(2-IS2)+IS2-1)
        !mpn BADT = (2-ISP)*B1+(ISP-1)*(B1-1)*NUAB(2)
        do A=A1,A2
          !mpn
          if (nga < ngb) then
            a_tmp = a-nind_ngbf
            b1_tmp = b1-nind_ngaf
          else
            a_tmp = a-nind_ngaf
            b1_tmp = b1-nind_ngbf
          end if
          !mpn
          !mpn AADT = IADT+KADT+BADT+A*(ISP-1)+(2-ISP)*(A-1)*NUAB(2)+IG-1
          ! T2 <A,B,J,K> for isp=1 T2 <B,A,K,J> for isp=2

          RAD = (KI-1)*adim*nstep+(A-A1)*nstep+1
          !mpn ISTEP = (ISP-1)*NUAB(2)+2-ISP
          !mpn
          if (nga >= ngb) then ! nga> ngb

            !mp write(u6,'(A,4(i5,2x),3x,i5)') '(I)  a_tmp,b1_tmp,k,i  nstep = ',a_tmp,b1_tmp,k,i,nstep
            ! T2(B1,A,I,K) =? T2(A,B1,K,I)
            do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim

              AADT = (I-1)*length1*length2*NOAB(2)+(K-1)*length1*length2+(b1_tmp-1+j_tmp)*length1+a_tmp

              x2(RAD+j_tmp) = t2_exp(AADT)

            end do

          else ! nga < ngb

            !mp write(u6,'(A,4(i5,2x),3x,i5)') '(II) b1_tmp,a_tmp,k,i   nstep = ',b1_tmp,a_tmp,i,k,nstep
            ! T2(B1,A,I,K)
            do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim

              AADT = (k-1)*length1*length2*NOAB(2)+(i-1)*length1*length2+(a_tmp-1)*length1+b1_tmp+j_tmp

              x2(RAD+j_tmp) = t2_exp(AADT)

            end do
          end if
          !mpn
          !mpn x2(RAD:RAD+NSTEP-1) = t2_exp(AADT:AADT+NSTEP-1)

        end do    ! A
      end do      ! I
    end do        ! K
    !!write(u6,'(A,3I5,4x,5D15.10)') 'block-v: a1,b1,IAS,ddot',a1/vblock+1,b1/vblock+1,ias, &
    !!                               ddot_(NNOAB(3)*adim*nstep,G(IX),1,G(IX),1),(G(I),I=IX,IX+3)
    !mp call multi_wridir(G(IX),NNOAB(3)*adim*nstep,LU,IAS,last)
    !mp
    !mpn do jjj=1,NNOAB(3)*adim*nstep
    !mpn   if (abs(x2(jjj)) > 1.0e5_wp) then
    !mpn     write(u6,*) 'fucko 1'
    !mpn     write(u6,*) jjj,x2(jjj)
    !mpn     stop
    !mpn   end if
    !mpn end do
    !mp
    call multi_wridir(x2,NNOAB(3)*adim*nstep,LU,IAS,last)
    ias = ias+iasblock
    !mpn
    call mma_deallocate(t2_exp)
    !mpn
  end do    ! B1
end do      ! A1
!mp write(u6,*) 'cast 3 ok'
!!close(LU)   in calling routine
!write(u6,*) FN,isp,IAS
!!dupblk(ndup) = last   in calling routine
!mp call w_memchk('all klvaa ')
!!call w_free(g(ix),0,'IT klvaa ')   in calling routine
!mp
call mma_deallocate(x2)
!mp
!mpn
if (printkey > 1) write(u6,*) 'VVOO integrals regenerated from MOLCAS'

!mpn
call xflush(u6)
!!stop

return

end subroutine klvaa_oovo
