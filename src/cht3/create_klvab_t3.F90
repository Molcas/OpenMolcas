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

subroutine create_klvab_t3(vblock)
!mp subroutine create_klvab_t3(G,vblock)
!
!  creates K(alpha-beta,alpha-beta),K(beta-alpha,alpha-beta)
!  DA files KMATBA and KMATAB, LMATBA and LMATAB
!  creates L(alpha-beta,alpha-beta),K(beta-alpha,alpha-beta)
!  max G at this place
!  both matrices have to be available at a time (unfortunately)
!
!  parallelization irrelevant at the moment
!  implemented integer offsets, PV, 14 may 2004.

use ChT3_global, only: DimGrpaR, ICH, IOPT, maxdim, nblock, nc, NNOAB, no, NOAB, NUAB, nv, NvGrp, printkey
use ChT3_procedures, only: klvaa_oovo
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: vblock
integer(kind=iwp) :: A, A1, A2, a_tmp, AADT, adim, B1, b1_tmp, b2_chk, i, i_blk, IADR, ias, ias_aa, iasblock, IJS, is2, isp, ISPA, &
                     ISTEP, IUHF, j, j_blk, j_tmp, k, last, last_aa, length, length1, length2, lu, n, nga, ngaf, ngal, ngb, ngbf, &
                     ngbl, nind_ngaf, nind_ngal, nind_ngbf, nind_ngbl, NNO, NNU, NSTEP, nug, RAD
character(len=6) :: FN
real(kind=wp), allocatable :: g(:), l0(:), l1(:), l1_1(:), l2_1(:), scr(:), t2_exp(:), t2_tmp(:), tmp(:), tmp2(:), x(:)

if (printkey >= 10) then
  write(u6,*)
  write(u6,*) '------ DimGrpaR ------'
  write(u6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
  write(u6,*)
end if

! - calculate overall memory requirements

call check_create_klvab_t3_mem(vblock)

!mp call w_rescope(G,'G create KL')
!mp call w_free(g,0,'G klvab ')
N = noab(1)+nuab(1)
!NNRED = nTri_Elem(NOAB(1))
IUHF = 0
if (IOPT(1) /= 0) then
  IUHF = 1
  !NNRED = NNOAB(3)
end if
LU = 98

!mp call w_alloc(it,NNOAB(3)*NNUAB(3),'IT klvab')
!mp call w_alloc(ix,vblock*vblock*n,'IX klvab')
!mpn call mma_allocate(t,NNOAB(3)*NNUAB(3),label='create_it')
call mma_allocate(x,vblock*vblock*n,label='c1_ix')
!mp
!!FN = 'T2OLDC'
!!call GET3DM(FN,G(it),NNUAB(3),NNRED,0)
!mp

!mp call w_alloc(it2_tmp,maxdim*maxdim*no*no,'ITMP klvab')
!mp call w_alloc(itmp,maxdim*maxdim*no*no,'ITMP klvab')
!mpn call mma_allocate(t2_tmp,maxdim*maxdim*no*no,label='create_it2')
!mpn call mma_allocate(tmp,maxdim*maxdim*no*no,label='create_tmp')
!
!mp call gather_t2(G(it),G(it2_tmp),G(itmp))
!mpn call gather_t2(t,t2_tmp,tmp)
!
!
!mp call w_free(G(it2_tmp),0,'ITMP klvab ')
!mpn call mma_deallocate(tmp)
!mpn call mma_deallocate(t2_tmp)

!mpn write(u6,*) 'T2 regenerated from MOLCAS'
!mpn write(u6,*)
!
!mp call dscal_(NNUAB(3)*NNRED,-One,G(it),1)
!mp call dscal_(NNUAB(3)*NNOAB(3),-One,G(it),1)
!mpn t(:) = -t
!mp if (IUHF == 0) call decomp2ind(G(it),NNUAB(3),noab(1),NUAB(1))
!!write(u6,*) ddot_(nnoab(3)*nnuab(3),G(it),1,G(it),1)

! number of blocks written in a single multiwrite

iasblock = vblock*vblock*N/nblock
if ((iasblock*nblock) < (vblock*vblock*N)) iasblock = iasblock+1

do isp=1,IUHF+1
  is2 = 3-isp
  FN = 'KMAT'//ich(3-isp)//ich(isp)
  write(u6,*) 'FN,LU=',FN,LU
  call multi_opendir(FN,LU)
  !ndup = ndup+1
  !if (ndup > ndupmx) then
  !  write(u6,*) 'create_klvab_t3 -- ndupmx exceeded'
  !  call abend()
  !end if
  !!write(u6,*) FN,isp,ndup
  !dupfil(ndup) = FN
  if (IUHF == 0) then
    FN(5:6) = ich(isp)//ich(isp)
    call multi_opendir(FN,LU+1)
    !!ndup = ndup+1
    !if (ndup > ndupmx) then
    !  write(u6,*) 'create_klvab_t3 -- ndupmx exceeded'
    !  call abend()
    !end if
    !!write(u6,*) FN,isp,ndup+1
    !dupfil(ndup+1) = FN
  end if
  ! currently using 3-dim (big field) - will be replaced after changing
  ! stepiv and the rest
  nnu = nTri_Elem(nuab(is2))
  !mp call w_alloc(ig,(nuab(isp)*nnu),'IG klvab')
  !mp call w_alloc(iscr,nuab(is2)*nuab(is2),'IG iscr')
  !mp
  call mma_allocate(g,nuab(isp)*nnu,label='c1_ig')
  call mma_allocate(scr,nuab(is2)*nuab(is2),label='c1_iscr')
  !mp
  IAS = 1
  IAS_AA = 1

  do K=1,noab(isp)
    !!if (IUHF == 1) then
    !!  KADT = (K-1)*NNUAB(3)*(NOAB(2)*(2-ISP)+ISP-1)
    !!else
    !!  KADT = 0
    !!end if

    !mp !FN = 'VVVAI'//ICH(ISP)
    !mp !call GET3DM(FN,G(IG),NNU,NUAB(ISP),K)
    if (printkey > 1) write(u6,*) 'Regenerating VVVo integrals for o = ',K
    !mp
    !mp call w_alloc(il1_1,nc*maxdim,'IL1_1 iscr')
    !mp call w_alloc(il2_1,nc*maxdim*maxdim,'IL2_1 iscr')
    !mp call w_alloc(itmp,max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim),'ITMP iscr')
    !mp
    call mma_allocate(l1_1,nc*maxdim,label='cc_il1_1')
    call mma_allocate(l2_1,nc*maxdim*maxdim,label='cc_il2_1')
    call mma_allocate(tmp,max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim),label='cc_itmp')
    !mp
    !mp call gen_vvvo(K,G(IG),G(il1_1),G(il2_1),G(itmp))
    call gen_vvvo(K,g,l1_1,l2_1,tmp)
    !mp write(u6,*) ddot_(nnu*nuab(isp),G(ig),1,G(ig),1)
    !mp call zeroma(G(ig),1,NNU*NUAB(ISP))
    !
    !mp call w_free(G(il1_1),0,'IL1_1 iscr')
    call mma_deallocate(tmp)
    call mma_deallocate(l2_1)
    call mma_deallocate(l1_1)
    !mp
    call delf(FN,K,K)
    if (iuhf == 0) call klvaa_vvv(x,g,vblock,N,nug,LU+1,last_aa,iasblock,K,ias_aa)
    !mpn if (iuhf == 0) call klvaa_vvv(ix,it,ig,iscr,vblock,N,nug,LU+1,last_aa,iasblock,K,ias_aa)
    !mp if (iuhf == 0) call klvaa_vvv(G,ix,it,ig,iscr,vblock,N,nug,LU+1,last_aa,iasblock,K,ias_aa)

    !!call xflush(u6)
    ! (c>d|AK)

    !mpn
    nga = 0
    !mpn
    do A1=1,NUAB(ISP),vblock
      A2 = A1+min(vblock,nuab(isp)-A1+1)-1
      adim = A2-A1+1
      nga = nga+1
      ngb = 0
      do B1=1,NUAB(IS2),vblock
        ngb = ngb+1
        NSTEP = min(vblock,nuab(is2)-B1+1)
        !!bdim = NSTEP
        IJS = (A1-1)*NNU+1
        !mpn
        !mpn write(u6,*)
        !mpn write(u6,*) '================================='
        !mpn write(u6,*) ' nga, ngb',nga,ngb
        !mpn write(u6,*) '================================='
        !mpn write(u6,*)

        ! - check the largest b2

        !mpn b2_chk = b1-1+vblock
        b2_chk = b1-1+min(vblock,nuab(is2)-B1+1)

        ! - find out which T2 blocked files will be needed
        !   for particular nga, ngb

        !mpn write(u6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk

        !switch = .false.
        if (nga < ngb) then
          !mpn write(u6,*) 'switching nga, ngb',ngb,nga
          !switch = .true.

          call block_interf(b1,b2_chk,a1,a2,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

        else

          call block_interf(a1,a2,b1,b2_chk,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

        end if

        !mpn write(u6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',ngaf,ngal,nind_ngaf,nind_ngal
        !mpn write(u6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',ngbf,ngbl,nind_ngbf,nind_ngbl

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

        !mpn write(u6,*) 'length1, vblock = ',length1,vblock
        !mpn write(u6,*) 'length2, vblock = ',length2,vblock

        length = length1*length2*no*no
        !mpn write(u6,*) 'length for blocked T2 amplitudes = ',length

        ! - setup memory

        !mpn write(u6,*) 'allocating t2_exp = ',length
        call mma_allocate(t2_exp,length,label='t2_exp')

        ! - read pertinent files and store them in the new blocked structure

        call mma_allocate(t2_tmp,maxdim*maxdim*no*no,label='cd_it2tmp')
        call mma_allocate(tmp,maxdim*maxdim*no*no,label='cd_itmp')

        call gather_t2_blocked(length1,length2,ngaf,ngal,ngbf,ngbl,t2_exp,t2_tmp,tmp)

        call mma_deallocate(tmp)
        call mma_deallocate(t2_tmp)

        t2_exp(:) = -t2_exp

        do A=A1,A2
          !mp call EXPA1_UHF(G(IJS),1,NUAB(IS2),1,G(ISCR))
          call EXPA1_UHF(g(IJS),1,NUAB(IS2),1,SCR)
          IJS = IJS+NNU
          ! Gix       <ng*ng|R,k>
          !! not needed B2 = B1+min(ng,nuab(is2)-B1+1))-1
          !mpn
          if (nga < ngb) then
            a_tmp = a-nind_ngbf
            b1_tmp = b1-nind_ngaf
          else
            a_tmp = a-nind_ngaf
            b1_tmp = b1-nind_ngbf
          end if
          !mp write(u6,'(A,2(i5,2x),A,i5,2x)') 'b1,   a     ,i,k = ',b1,a,'    i',k
          ! mv T2(B,A,I,K) >> G(ix)
          !mpn
          do I=1,NOAB(IS2)
            ISPA = ISP
            !!if (IUHF == 1) then
            !!  IADT = (I-1)*NNUAB(3)*(NOAB(2)*(2-IS2)+IS2-1)
            !!else
            !!  IADT = iTri(K,I)
            !!  IADT = (IADT-1)*NNUAB(3)
            !!  if (K < I) ISPA = IS2
            !!end if
            !mpn BADT = (2-ISPA)*B1+(ISPA-1)*(B1-1)*NUAB(2)
            !mpn!!!
            !mpn AADT_tmp = IADT+KADT+BADT+A*(ISPA-1)+(2-ISPA)*(A-1)*NUAB(2)
            !mpn!!!
            ! T2 <A,B,J,K> for isp=1 T2 <B,A,K,J> for isp=2
            !! RAD = (I-1)*vblock*vblock+(A-A1)*vblock+1
            RAD = (I-1)*adim*nstep+(A-A1)*nstep+1
            ISTEP = (ISPA-1)*NUAB(2)+2-ISPA
            !mpn
            if (nga >= ngb) then ! nga> ngb

              !mp write(u6,'(A,4(i5,2x),3x,i5)') '(I)  a_tmp,b1_tmp,k,i  nstep = ',a_tmp,b1_tmp,k,i,nstep
              ! T2(B,A,I,K) =? T2(A,B1,K,I)
              do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim

                AADT = (I-1)*length1*length2*NOAB(2)+(K-1)*length1*length2+(b1_tmp-1+j_tmp)*length1+a_tmp

                x(RAD+j_tmp) = t2_exp(AADT)
                !mp
                !mpn if (abs(t2_exp(AADT)-t(AADT_tmp+j_tmp)) > 1.0e-5_wp) then
                !mpn   write(u6,*) 'halohaha 1',AADT,AADT_tmp+j_tmp,t2_exp(AADT),t(AADT_tmp+j_tmp)
                !mpn   stop
                !mpn end if
                !mp
              end do

            else ! nga < ngb

              !mp write(u6,'(A,4(i5,2x),3x,i5)') '(II) b1_tmp,a_tmp,k,i   nstep = ',b1_tmp,a_tmp,i,k,nstep
              ! T2(B,A,I,K)
              do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim

                !mp AADT = (k-1)*length1*length2*NOAB(2)+(i-1)*length1*length2+(a_tmp-1)*length1+b1_tmp+j_tmp
                !mp !AADT = (i-1)*length1*length2*NOAB(2)+(k-1)*length1*length2+(a_tmp-1)*length1+b1_tmp+j_tmp
                AADT = (K-1)*length1*length2*NOAB(2)+(I-1)*length1*length2+(a_tmp-1)*length1+b1_tmp+j_tmp

                x(RAD+j_tmp) = t2_exp(AADT)
                !mp
                !mpn if (abs(t2_exp(AADT)-t(AADT_tmp+j_tmp)) > 1.0e-5_wp) then
                !mpn   write(u6,*) 'halohaha 2',AADT,AADT_tmp+j_tmp,t2_exp(AADT),t(AADT_tmp+j_tmp)
                !mpn   stop
                !mpn end if
                !mp

              end do
            end if

            !mpn call dcopy_(NSTEP,t2_exp(AADT),ISTEP,x(RAD),1)
            !mpn
          end do      ! I
          RAD = noab(is2)*adim*nstep+(A-A1)*nstep+1
          do IADR=B1,NUAB(IS2)*NUAB(IS2),NUAB(IS2)
            !mp call dcopy_(NSTEP,G(IADR),1,G(RAD),1)
            x(RAD:RAD+NSTEP-1) = scr(IADR:IADR+NSTEP-1)
            RAD = RAD+adim*nstep
          end do      ! IADR
        end do        ! A
        !!write(u6,'(A,4I5,4x,D15.10)') 'block-w: K,a1,b1,IAS,ddot',K,a1,b1,ias,ddot_(N*vblock*vblock,G(IX),1,G(IX),1)
        !mp call multi_wridir(G(IX),N*vblock*vblock,LU,IAS,last)

        !mp
        !mp !!do jjj = 1,N*vblock*vblock
        !mp !!if (abs(x(jjj)) > 1.0e5_wp) then
        !mp !!  write(u6,*) 'prasa 1 ',jjj,x(jjj)
        !mp !!  stop
        !mp !!end if
        !mp !!end do
        !mp
        call multi_wridir(x,N*vblock*vblock,LU,IAS,last)
        !!write(u6,*) 'N*vblock*vblock,LU,last ',N*vblock*vblock,LU,last
        ias = ias+iasblock
        !mp
        call mma_deallocate(t2_exp)
        !mp
      end do      ! B1
    end do        ! A1
  end do          ! K
  if (printkey > 1) then
    write(u6,*) 'VVVo integrals regenerated from MOLCAS'
    write(u6,*)
  end if
  !mp
  close(LU)
  !dupblk(ndup) = last
  if (IUHF == 0) then
    close(LU+1)
    !ndup = ndup+1
    !dupblk(ndup) = last_aa
  end if
  !write(u6,*) FN,isp,IAS
  !mp call w_memchk('IG klvab ')
  !mp call w_free(g(ig),0,'IG klvab ')
  !mp
  call mma_deallocate(scr)
  call mma_deallocate(g)
  !mp
end do     ! ISP

!mp call dscal_(NNUAB(3)*NNOAB(3),-One,G(it),1)
!mpn t(:) = -t
!mp??
call mma_deallocate(x)
!mp??
do isp=1,IUHF+1
  !mp call w_memchk('IX klvab ')
  !mp call w_free(g(ix),0,'IX klvab ')
  !mp
  is2 = 3-isp
  iasblock = nnoab(3)*vblock*N/nblock
  if ((iasblock*nblock) < (nnoab(3)*vblock*N)) iasblock = iasblock+1
  FN = 'LMAT'//ich(3-isp)//ich(isp)
  call multi_opendir(FN,LU)
  !ndup = ndup+1
  !if (ndup > ndupmx) then
  !  write(u6,*) 'create_klvab_t3 -- ndupmx exceeded'
  !  call abend()
  !end if
  !write(u6,*) FN,isp,ndup
  !dupfil(ndup) = FN

  FN = 'OOVAI'//ICH(ISP)
  !mp call w_alloc(ix,noab(isp)*noab(IS2)*vblock*n,'IX klvabo')
  call mma_allocate(x,noab(isp)*noab(IS2)*vblock*n,label='c2_ix')
  nno = nTri_Elem(noab(is2))
  !mp call w_alloc(ig,noab(isp)*nuab(isp)*nno,'IG klvabo')
  !mp call w_alloc(iscr,vblock*noab(IS2)*noab(IS2),'ISCRo klvabo')
  call mma_allocate(g,noab(isp)*nuab(isp)*nno,label='c2_ig')
  call mma_allocate(scr,vblock*noab(IS2)*noab(IS2),label='c2_iscr')
  !mp
  !mp !call GET3DM(FN,G(ig),NNO,NUAB(ISP)*NOAB(ISP),0)
  !mp
  !mp call w_alloc(il0,nc*nno,'IL0 klvabo')
  !mp call w_alloc(il1,nc*no*nv,'IL1 klvabo')
  !mp call w_alloc(itmp,max(nc*nno,nc*no*maxdim,nc*no*nv),'ITMP klvabo')
  !mp
  call mma_allocate(l0,nc*nno,label='cr_il0')
  call mma_allocate(l1,nc*no*nv,label='cr_il1')
  call mma_allocate(tmp,max(nc*nno,nc*no*maxdim,nc*no*nv),label='cr_itmp')
  !mp
  call gen_oovo(g,l0,l1,tmp)
  !mp call gen_oovo(G(ig),G(il0),G(il1),G(itmp))
  !mp write(u6,*) ddot_(NNO*NUAB(ISP)*NOAB(ISP),G(ig),1,G(ig),1)
  !mp call zeroma(G(ig),1,NNO*NUAB(ISP)*NOAB(ISP))

  !mp call w_free(G(il0),0,'IL0 klvab')
  call mma_deallocate(tmp)
  call mma_deallocate(l1)
  call mma_deallocate(l0)

  if (printkey > 1) write(u6,*) 'OOVO integrals regenerated from MOLCAS'
  !mp
  IAS = 1
  !mpn
  nga = 0
  !mpn
  do A1=1,NUAB(ISP),vblock
    A2 = A1+min(vblock,nuab(isp)-(A1-1))-1
    NSTEP = min(vblock,nuab(isp)-(A1-1))
    !mpn
    nga = nga+1

    !mp write(u6,*)
    !mp write(u6,*) '================================='
    !mp write(u6,*) ' nga ',nga
    !mp write(u6,*) '================================='
    !mp write(u6,*)

    !mp write(u6,'(A,2(i5,2x))') 'b1,b2 = ',a1,a2

    !mp call block_interf(1,1,a1,a2,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)
    call block_interf(1,nuab(1),a1,a2,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

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

    length = length1*length2*no*no
    !mp write(u6,*) 'length for blocked T2 amplitudes = ',length

    ! - setup memory

    !mp write(u6,*) 'allocating t2_exp = ',length
    call mma_allocate(t2_exp,length,label='t2_exp')

    ! - read pertinent files and store them in the new blocked structure

    call mma_allocate(t2_tmp,maxdim*maxdim*no*no,label='t2_tmp')
    call mma_allocate(tmp2,maxdim*maxdim*no*no,label='tmp2')

    call gather_t2_fblocked(length1,length2,ngbf,ngbl,t2_exp,t2_tmp,tmp2)

    call mma_deallocate(tmp2)
    call mma_deallocate(t2_tmp)

    !mpn
    do k=1,noab(isp)
      !!if (IUHF == 1) then
      !!  KADT = (K-1)*NNUAB(3)*(NOAB(2)*(2-ISP)+ISP-1)
      !!else
      !!  KADT = 0
      !!end if
      IJS = 1+(A1-1)*nno+(k-1)*nuab(isp)*nno
      !mp call EXPA1_UHF(G(IJS),nstep,NOAB(IS2),1,G(iscr))
      call EXPA1_UHF(g(IJS),nstep,NOAB(IS2),1,scr)
      do I=1,noab(is2)
        ISPA = ISP
        !!if (IUHF == 1) then
        !!  IADT = (I-1)*NNUAB(3)*(NOAB(2)*(2-IS2)+IS2-1)
        !!else
        !!  IADT = iTri(K,I)
        !!  IADT = (IADT-1)*NNUAB(3)
        !!  if(K < I) ISPA = IS2
        !!end if

        do A=A1,A2

          !mpn
          a_tmp = a-nind_ngbf

          AADT = (K-1)*noab(2)*NUAB(2)*length2+(i-1)*NUAB(2)*length2+(a_tmp-1)*NUAB(2)+1
          !mpn
          ! copies T
          !mpn BADT = 2-ISPA+(ISPA-1)*NUAB(2)
          !  T2 <A,B,J,K> for isp=1 T2 <B,A,K,J> for isp=2
          ISTEP = (ISPA-1)*NUAB(2)+2-ISPA
          RAD = (I-1)*nstep*n+(A-A1)*n+noab(is2)+(K-1)*noab(is2)*nstep*N+1
          !mp call dcopy_(nuab(is2),G(AADT),ISTEP,G(RAD),1)
          !mp
          call dcopy_(nuab(is2),t2_exp(AADT),ISTEP,x(RAD),1)
        end do     ! A
      end do       ! I
      ! copies OOVO
      do J=1,noab(is2)
        do A=1,NSTEP
          RAD = 1+(A-1)*n+(J-1)*nstep*N+(K-1)*noab(is2)*nstep*N
          IADR = 1+(A-1)*noab(is2)*noab(is2)+(j-1)*noab(is2)
          !mp call dcopy_(noab(is2),G(IADR),1,G(RAD),1)
          x(RAD:RAD+noab(is2)-1) = scr(IADR:IADR+noab(is2)-1)
        end do     ! A
      end do       ! J
    end do         ! K
    !!write(u6,'(A,2I5,4x,D15.10)') 'block-w:a1,IAS,ddot',a1,ias,ddot_(N*vblock*nnoab(3),g(ix),1,g(ix),1)
    !!write(u6,'(a,a,2I4,D16.8)') 'block-w',ich(isp),(A1/vblock)+1,IAS,ddot_(N*nnoab(3)*NSTEP,g(ix),1,g(ix),1)
    !mp call multi_wridir(G(IX),N*nstep*nnoab(3),LU,IAS,last)
    !mp
    !mp !!do jjj=1,N*nstep*nnoab(3)
    !mp !!  if (abs(x(jjj)) > 1.0e5_wp) then
    !mp !!    write(u6,*) 'prasa 2 ',jjj,x(jjj)
    !mp !!    stop
    !mp !!  end if
    !mp !!end do
    !mp
    call multi_wridir(x,N*nstep*nnoab(3),LU,IAS,last)
    !mp
    IAS = IAS+iasblock
    !mpn
    call mma_deallocate(t2_exp)
    !mpn
  end do      ! A1
  close(LU)
  call mma_deallocate(scr)
  !write(u6,*) FN,isp,IAS
  !dupblk(ndup) = last
  if (IUHF == 0) then
    FN = 'LMAT'//ich(isp)//ich(isp)
    call multi_opendir(FN,LU)
    !ndup = ndup+1
    !if (ndup > ndupmx) then
    !  write(u6,*) 'create_klvab_t3 -- ndupmx exceeded'
    !  call abend()
    !end if
    !write(u6,*) FN,isp,ndup
    !dupfil(ndup) = FN
    ! this is to ensure correct copy to slaves
    IAS_AA = 1
    !mp write(u6,*) 'test 1 na iscr ',vblock*noab(IS2)*noab(IS2)
    !mp write(u6,*) 'test 1 na ig ',noab(isp)*nuab(isp)*nno
    !mp write(u6,*) 'test 1 na ix ',noab(isp)*noab(IS2)*vblock*n
    !mp call klvaa_oovo(G,ix,it,ig,iscr,vblock,N,nug,LU,last_aa,ias_aa)
    !mpn call klvaa_oovo(ix,it,ig,iscr,vblock,N,nug,LU,last_aa,ias_aa)
    call klvaa_oovo(x,g,vblock,N,nug,LU,last_aa,ias_aa)
    !mp
    !mp write(u6,*) 'klvaa_oovo finished'
    !mp
    close(LU)
    !ndup = ndup
    !dupblk(ndup) = last_aa
  end if
end do     ! ISP
!mp call w_memchk('all klvab ')
!mp call w_free(g(it),0,'IT klvab ')
! ig sa odalokuju v klva_oovo
!mp @@@ call GetMem('create_ig','Free','Real',ig,noab(isp)*nuab(isp)*nno)
!mp ???
!mp?? call GetMem('c2_ix','Free','Real',ix,noab(isp)*noab(IS2)*vblock*n)
!mpn call mma_deallocate(t)
call xflush(u6)

return

end subroutine create_klvab_t3
