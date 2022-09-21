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

subroutine klvaa_vvv(ix,ig,vblock,N,nug,lu,last,iasblock,K,ias)
!mp subroutine klvaa_vvv(G,ix,it,ig,iscr,vblock,N,nug,lu,last,iasblock,K,ias)
!mpn subroutine klvaa_vvv(ix,it,ig,iscr,vblock,N,nug,lu,last,iasblock,K,ias)
!
!  creates K(alpha>alpha,alpha-alpha) or K(beta>beta,beta,beta)
!  DA files KMATICH(ISP)ICH(ISP)     ISP=A
!  generates antisymmetrized integrals from abab (T2) or bbaa (vvvo)
!  ix, it, ig  allocated in create_klvab_t3
!  K,ias defined in create_klvab_t3.f
!
!  parallelization (seems to be) irrelevant at the moment
!  implemented integer offsets, PV, 14 may 2004.

use Index_Functions, only: iTri
use Constants, only: One
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ix, ig, vblock, N, nug, lu, last, iasblock, K, ias
#include "WrkSpc.fh"
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "ccsd_t3compat.fh"
#include "uhf.fh"
integer(kind=iwp) :: A, A1, A2, a_tmp, AADT, ADIM, B, B1, b1_tmp, B2, b2_chk, i_blk, IJS, isp, it2_tmp, it_exp, itmp, j_blk, &
                     j_tmp, KADT, KI, length, length1, length2, MAXDIMM, NGA, ngaf, ngal, NGB, ngbf, ngbl, nind_ngaf, nind_ngal, &
                     nind_ngbf, nind_ngbl, NNU, NSTEP, R
logical(kind=iwp) :: switch

ISP = 1
nug = nuab(isp)/vblock
if ((vblock*nug) /= nuab(isp)) nug = nug+1
NNU = NUAB(1)*(NUAB(1)+1)/2
if (printkey >= 11) then
  write(u6,'(A,10I5)') 'entering klvaa_vvv',vblock,N,nug,lu,last,iasblock,K,ias
  !mpn write(u6,*) ddot_(nnoab(3)*nnuab(3),Work(it),1,Work(it),1)
  !mp write(u6,*) ddot_(nnoab(3)*nnuab(3),G(it),1,G(it),1)
end if
!!do K=1,noab(isp)
!!  FN = 'VVVOI'//ICH(ISP)
!!  call GET3DM(FN,G(ig),NNUAB(isp),NUAB(ISP),K)
!!  call delf(FN,K,K)  ! all done in create_klvab_t3.f
do nga=1,nug
  A1 = (nga-1)*vblock+1
  adim = min(vblock,nuab(isp)-A1+1)
  A2 = A1+adim-1
  do ngb=1,nga
    if (nga == ngb) then
      maxdimm = adim*(adim-1)/2
    else
      maxdimm = adim*vblock
    end if
    call zeroma(Work(IX),1,N*maxdimm)
    !mp call zeroma(G(IX),1,N*maxdim)
    B1 = (ngb-1)*vblock+1

    !mp !write(u6,*)
    !mp !write(u6,*) '================================='
    !mp !write(u6,*) 'nga, ngb',nga,ngb
    !mp !write(u6,*) '================================='
    !mp !write(u6,*)

    !mpn
    ! - check the largest b2

    !mp b2_chk = b1-1+min(vblock,a2-b1)
    !mpn b2_chk = b1-1+vblock
    b2_chk = b1+min(vblock,nuab(isp)-b1+1)-1

    ! - find out which T2 blocked files will be needed
    !   for particular nga, ngb

    !mp !write(u6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk

    call block_interf(a1,a2,b1,b2_chk,ngaf,ngal,nind_ngaf,nind_ngal,ngbf,ngbl,nind_ngbf,nind_ngbl)

    !mp !write(u6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',ngaf,ngal,nind_ngaf,nind_ngal
    !mp !write(u6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',ngbf,ngbl,nind_ngbf,nind_ngbl
    !mp !write(u6,*) 'maxdim = ',maxdim

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

    !mp !write(u6,*) 'length1, vblock = ',length1,vblock
    !mp !write(u6,*) 'length2, vblock = ',length2,vblock

    length = length1*length2*no*no
    !mp !write(u6,*) 'length for blocked T2 amplitudes = ',length

    ! - setup memory

    !mp !write(u6,*) 'allocating t2_exp = ',length
    call GetMem('it1_exp','Allo','Real',it_exp,length)

    ! - read pertinent files and store them in the new blocked structure

    call GetMem('c1_it2tmp','Allo','Real',it2_tmp,maxdim*maxdim*no*no)
    call GetMem('c1_itmp','Allo','Real',itmp,maxdim*maxdim*no*no)

    switch = .false.
    call gather_t2_blocked(length1,length2,ngaf,ngal,ngbf,ngbl,Work(it_exp),Work(it2_tmp),Work(itmp),switch)

    call GetMem('c1_itmp','Free','Real',itmp,maxdim*maxdim*no*no)
    call GetMem('c1_it2tmp','Free','Real',it2_tmp,maxdim*maxdim*no*no)

    call dscal_(length,-One,Work(it_exp),1)

    !mp
    do a=a1,a2
      B2 = B1-1+min(vblock,A-B1)
      NSTEP = B2-B1+1
      !mp !write(u6,*) 'NSTEP = ',NSTEP
      if (nstep /= 0) then
        if (nga == ngb) then
          IJS = (a-a1-1)*(a-a1)/2+IX
          !!IJS = (a-a1-2+1)*(a-a1-1+1)/2+1
        else
          IJS = (a-a1)*vblock+IX
        end if
        do R=1,K-1
          ! copies (b1..b2,a,r,k)-(b1..b2,a,k,r)
          !!KADT = (K-2)*(K-1)/2+R
          !!KADT = (KADT-1)*NNUAB(ISP)
          ! address and block for the A1-A2x B1-B2
          !!KADT = KADT+(a-1)*(a-2)/2 +B1 +IT -1
          ! T2 >>> K
          !!call dcopy_(NSTEP,G(KADT),1,G(IJS),1)
          !mpn KADT = (k-1)*noab(1)*NNUAB(3)+(r-1)*NNUAB(3)+(a-1)*nuab(1)+B1+IT-1
          !mpn AADT = (r-1)*noab(1)*NNUAB(3)+(k-1)*NNUAB(3)+(a-1)*nuab(1)+B1+IT-1

          a_tmp = a-nind_ngaf
          b1_tmp = b1-nind_ngbf
          !mp !write(u6,'(A,4(i5,2x))') 'b1,   a     ,r,k = ',b1,a,r,k
          !mp !write(u6,'(A,4(i5,2x))') 'a_tmp,b1_tmp,k,r = ',a_tmp,b1_tmp,k,r

          do j_tmp=0,NSTEP-1

            KADT = (r-1)*noab(1)*length1*length2+(k-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+a_tmp+it_exp-1

            AADT = (k-1)*noab(1)*length1*length2+(r-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+a_tmp+it_exp-1

            Work(IJS+j_tmp) = Work(AADT)-Work(KADT)
            !mp !write(u6,'(A,2(i5,2x),3(f18.10,2x))') 'XXXN = ',a_tmp,b1_tmp+j_tmp-1,Work(IJS+j_tmp),Work(AADT),Work(KADT)

          end do

          !mp call vsub(G(KADT),1,G(AADT),1,G(IJS),1,NSTEP)
          !mpn call vsub(Work(KADT),1,Work(AADT),1,Work(IJS),1,NSTEP)
          !!write(u6,'(A,3I3,8D15.8)') 'T',K,R,a,(G(I),I=IJS,IJS+NSTEP-1)
          IJS = IJS+maxdimm
        end do      ! R
          !mp call zeroma(G(IJS),1,NSTEP)
        call zeroma(Work(IJS),1,NSTEP)
        IJS = IJS+maxdimm
        ! T2 >>> K  (transposed)
        ! copies (b1..b2,a,r,k)-(b1..b2,a,k,r)
        !mpn

        ! (b1..b2,a,r,k) =? (a,b1...b2,k,r)
        ! (b1..b2,a,k,r) =? (a,b1...b2,r,k)

        !mpn
        do R=K+1,NOAB(ISP)
          !!KADT = (R-2)*(R-1)/2+K
          !!KADT = (KADT-1)*NNUAB(ISP)
          !!KADT = KADT+(a-1)*(a-2)/2+B1+IT-1
          !!call vneg_cht3(G(KADT),1,G(IJS),1,NSTEP)

          !mpn KADT = (k-1)*noab(1)*NNUAB(3)+(r-1)*nnUab(3)+(a-1)*nuab(1)+B1+IT-1
          !mpn AADT = (r-1)*noab(1)*NNUAB(3)+(k-1)*nnUab(3)+(a-1)*nuab(1)+B1+IT-1

          a_tmp = a-nind_ngaf
          b1_tmp = b1-nind_ngbf
          !mp !write(u6,'(A,4(i5,2x))') 'b1,   a     ,r,k = ',b1,a,r,k
          !mp !write(u6,'(A,4(i5,2x))') 'a_tmp,b1_tmp,k,r = ',a_tmp,b1_tmp,k,r

          do j_tmp=0,NSTEP-1

            KADT = (r-1)*noab(1)*length1*length2+(k-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+a_tmp+it_exp-1

            AADT = (k-1)*noab(1)*length1*length2+(r-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+a_tmp+it_exp-1

            Work(IJS+j_tmp) = Work(AADT)-Work(KADT)
            !mp write(u6,'(A,2(i5,2x),3(f18.10,2x))') 'XXXN = ',a_tmp,b1_tmp+j_tmp-1,Work(IJS+j_tmp),Work(AADT),Work(KADT)

          end do

          !mp call vsub(G(KADT),1,G(AADT),1,G(IJS),1,NSTEP)
          !mpn call vsub(Work(KADT),1,Work(AADT),1,Work(IJS),1,NSTEP)
          !!write(u6,'(A,3I3,8D15.8)') 'T',K,R,a,(G(I),I=IJS,IJS+NSTEP-1)
          !!call daxpy_(NSTEP,-One,G(KADT),1,G(IJS),1)
          IJS = IJS+maxdimm
        end do      ! R
        ! <VV|VO> >>> K original for aaaa
        ! (VV|VO( >>> K now          aabb needed (AR|BK)-(BR|AK)

        !!call EXPA1_UHF(G(IJS),1,NUAB(IS2),1,G(ISCR))
        !!KADT = IG-1+(a-1)*(a-2)/2 +B1
        KADT = IG-1+(B1-1)*NNU
        AADT = IG-1+(A-1)*NNU
        do R=1,NUAB(ISP)
          !!call dcopy_(NSTEP,G(KADT),1,G(IJS),1)
          ! first (AR|BK)
          !mp call dcopy_(NSTEP,G(KADT+iTri(a,r)),NNU,G(IJS),1)
          call dcopy_(NSTEP,Work(KADT+iTri(a,r)),NNU,Work(IJS),1)
          ! now  -(BR|AK)
          KI = IJS
          do B=B1,B2
            !mp G(KI) = G(KI)-G(AADT+iTri(B,R))
            Work(KI) = Work(KI)-Work(AADT+iTri(B,R))
            KI = KI+1
          end do
          !!write(u6,'(A,3I3,8D15.8)') 'V',K,R,a,(G(I),I=IJS,IJS+NSTEP-1)
          IJS = IJS+maxdimm
          !!KADT = KADT+NNUAB(ISP)
        end do   !R
      end if     ! nstep == 0
    end do       !A
    !  <VV|VO> >>> K
    !mp !write(u6,'(A,5I5,4x,5D15.10)') 'block-w: K,a1,b1,IAS,ddot',K,a1,b1,ias,maxdim,ddot_(N*maxdim,Work(IX),1,Work(IX),1), &
    !                                   (Work(I),I=IX,IX+3)
    !mp !write(u6,'(A,5I5,4x,5D15.10)') 'block-w: K,a1,b1,IAS,ddot',K,a1,b1,ias,maxdim,ddot_(N*maxdim,G(IX),1,G(IX),1), &
    !                                   (G(I),I=IX,IX+3)
    !mp  call multi_wridir(G(IX),N*maxdim,LU,IAS, last)
    !mp
    !mp !!do jjj=0,N*maxdimm-1
    !mp !!  if (abs(Work(ix+jjj)) > 1.0e5_wp) then
    !mp !!    write(u6,*) 'prasa 3 ',jjj,Work(ix+jjj)
    !mp !!    stop
    !mp !!  end if
    !mp !!end do
    !mp
    call multi_wridir(Work(IX),N*maxdimm,LU,IAS,last)
    ias = ias+iasblock
    !mpn
    !mp write(u6,*) 'deallocating t2_exp = ',length
    call GetMem('it1_exp','Free','Real',it_exp,length)
    !mpn
  end do       ! ngb
end do         ! nga
!!end do         ! K
!mp !write(u6,'(A,10I5)')'leaving klvaa_vvv',vblock,N,nug,lu,last,iasblock,K,ias
call xflush(u6)
!!stop

return

end subroutine klvaa_vvv
