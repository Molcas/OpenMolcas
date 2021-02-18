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
cmp!      SUBROUTINE klvaa_oovo(G,ix,it,ig,iscr,vblock,N,nug,
cmpn      SUBROUTINE klvaa_oovo(ix,it,ig,iscr,vblock,N,nug,
      SUBROUTINE klvaa_oovo(ix,ig,iscr,vblock,N,nug,
     $ LU,last,ias)
C
C  creates L(alpha>alpha,alpha-alpha)
C  DA files LMATICH(ISP)ICH(ISP)
C  max G at this place
C
C  parallelization (seems to be) irrelevant at the moment
C  implemented integer offsets, PV, 14 may 2004.
C
      IMPLICIT NONE
#include "ndisk.fh"
#include "dupfiles.fh"
cmp
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
#include "WrkSpc.fh"
cmpn
#include "cht3_reord.fh"
        integer i_blk,j_blk,b2_chk
        integer ngaf,ngal,ngbf,ngbl
        integer nind_ngbf,nind_ngbl,nind_ngaf,nind_ngal
        integer length,length1,length2
        integer it_exp,it2_tmp,itmp2
        integer ngaa,a_tmp,b1_tmp
        integer j_tmp,RAD_tmp
cmp
c       integer jjj
cmp
cmpn
cmp
cmp      real*8 G(*),ddot_,one
      real*8 one
c     real*8 ddot_
      parameter (one=1.d0)
      integer ix,ig,iscr, KADT, IJS, RAD, AADT
      integer isp,ias,vblock,n,i,j,k,lu,iasblock
      INTEGER A,A1,A2,B1,B2,NSTEP,ADIM,NUG,NGA,NGB,R,MAXDIM2
      INTEGER IS2,KI, last,indab
      CHARACTER ich*1
c      CHARACTER FN*6
      INTEGER IOPT,NOAB,NNOAB,NUAB,NNUAB,nno
c     INTEGER m
cmp
        integer il1,itmp,il2
cmp
      COMMON/UHF/NOAB(2),NNOAB(3),NUAB(2),NNUAB(3),ICH(3)
      COMMON/IOIND/IOPT(96)
      indab(i,j)=(max(i,j)-1)*max(i,j)/2+min(i,j)
c
C
cmp      write(6,*) 'Entering klvaa_oovo'
C
      ISP=1
      iasblock=NNOAB(ISP)*vblock*N/nblock
      if((iasblock*nblock).lt.(NNOAB(ISP)*vblock*N))iasblock=iasblock+1
      nno=(noab(isp)+1)*noab(isp)/2
      !!call w_alloc(ix,nnoab(ISP)*vblock*n,'IX klvaao ')
      !!call w_alloc(ig,noab(isp)*nuab(isp)*nnoab(isp),'IG klvaao')
      !!call w_alloc(iscr,nuab(isp)*nuab(isp)*nnoab(isp),'ISCR klvaao')
      !FN(1:5)='OOVOI'
      !FN(6:6)=ICH(ISP)
      !!CALL GET3DM(FN,G(ig),NUAB(ISP)*NOAB(ISP),NNOAB(ISP),0)
C
C      call EXPA1_UHF(G(IT),nnoab(isp),NUAB(ISP),-1,G(ISCR))
C       expa done here. remains in it address
cmpn      IJS=IT
cmpn       do i=2,noab(isp)
cmpn       do j=1,i-1
cmpn       KADT=IT+(i-1)*noab(isp)*nnuab(3)+(j-1)*nnuab(3)
cmpn       !!call dcopy_(NNUAB(3),G(KADT),1,G(IJS),1)
cmpncmp!       call vneg_cht3(G(KADT),1,G(IJS),1,NNUAB(3))
cmpn       call vneg_cht3(Work(KADT),1,Work(IJS),1,NNUAB(3))
cmpncmp!       CALL TRANSM_A(G(KADT),G(IJS),NUAB(ISP),NUAB(ISP))
cmpn       CALL TRANSM_A(Work(KADT),Work(IJS),NUAB(ISP),NUAB(ISP))
cmpn       IJS=IJS+NNUAB(3)
cmpn       enddo ! j
cmpn       enddo ! i
c
      IAS=1
        ngaa=0
cmpn
      do A1=1,NUAB(ISP),vblock
C  not needed    call zeroma(g(ix),1,nnoab(isp)*vblock*n)
         adim=min(vblock,nuab(isp)-A1+1)
         A2=A1+adim-1
cmpn
        ngaa=ngaa+1
c
cmp        write (6,*)
cmp        write (6,*) '================================='
cmp        write (6,*) ' nga ',ngaa
cmp        write (6,*) '================================='
cmp        write (6,*)
c
cmp        write (6,'(A,2(i5,2x))') 'b1,b2 = ',a1,a2
c
        call block_interf(1,1,a1,a2,
     &   ngaf,ngal,nind_ngaf,nind_ngal,
     &   ngbf,ngbl,nind_ngbf,nind_ngbl)
c
cmp        write (6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',
cmp     & ngbf,ngbl,nind_ngbf,nind_ngbl
c
c - read amplitudes T2(nv,vblock,j<i)
c
c - calculate memory requirements
c
        length1=nuab(isp)
c
          length2=0
        do i_blk=ngbf,ngbl
          length2=length2+DimGrpaR(i_blk)
        end do
c
cmp        write (6,*) 'length1, NUAB   = ',length1,NUAB(1)
cmp        write (6,*) 'length2, vblock = ',length2,vblock
c
        length=length1*length2*(((no-1)*no)/2)
cmp        write (6,*) 'length for blocked T2 amplitudes = ',length
c
c - setup memory
c
cmp        write (6,*) 'allocating t2_exp = ',length
        call GetMem('it4_exp','Allo','Real',it_exp,length)
c
c - read pertinent files and store them in the new blocked structure
c
        call GetMem('it2_tmp','Allo','Real',it2_tmp,
     & maxdim*maxdim*no*no)
        call GetMem('itmp2','Allo','Real',itmp2,
     & maxdim*maxdim*no*no)
c
        call gather_t2anti_blocked(length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & Work(it_exp),Work(it2_tmp),Work(itmp2))
c
        call GetMem('itmp2','Free','Real',itmp2,
     & maxdim*maxdim*no*no)
        call GetMem('it2_tmp','Free','Real',it2_tmp,
     & maxdim*maxdim*no*no)

cmpn
C     copies oovo
C  oovo (k>l,ai)
           K=0
           do I=2,NOAB(ISP)
           do J=1,I-1
           K=K+1
         !!do K=1,NNOAB(ISP)
            IJS=IX+(K-1)*adim*n
            !!KADT=(a1-1)*noab(isp)+IG+(K-1)*noab(isp)*nuab(isp)
            do a=A1,A2
            KADT=(J-1)*nno*nuab(isp)+(A-1)*nno + IG -1
            do r=1,noab(isp)
C    <ra|ij>
cmp            G(IJS+r-1)=G(KADT+indab(r,i))
            Work(IJS+r-1)=Work(KADT+indab(r,i))
            enddo
            KADT=(I-1)*nno*nuab(isp)+(A-1)*nno + IG -1
            do r=1,noab(isp)
C    -<ra|ji>
cmp            G(IJS+r-1)=G(IJS+r-1)-G(KADT+indab(r,j))
            Work(IJS+r-1)=Work(IJS+r-1)-Work(KADT+indab(r,j))
            enddo
!!               call dcopy_(noab(isp),G(KADT),1,G(IJS),1)
               IJS=IJS+N
!!               KADT=KADT+noab(isp)
            enddo
C now the T2
cmpn
            !!KADT=ISCR+(K-1)*NUAB(ISP)*NUAB(ISP) +(A1-1)*NUAB(ISP)
cmp            KADT=IT+(K-1)*NUAB(ISP)*NUAB(ISP) +(A1-1)*NUAB(ISP)
            IJS=NOAB(ISP)+IX+(K-1)*adim*n
c
        a_tmp=a1-nind_ngbf
cmp        write (6,*) 'K, a1, a_tmp',k,a1,a_tmp
c
        KADT=it_exp+(K-1)*NUAB(ISP)*length2
     & +(a_tmp-1)*NUAB(ISP)
c
            do a=A1,A2
cmp               call dcopy_(NUAB(isp),G(KADT),1,G(IJS),1)
               call dcopy_(NUAB(isp),Work(KADT),1,Work(IJS),1)
cmp        write (6,*) (Work(KADT+a_tmp),a_tmp=0,NUAB(isp)-1)
      !!write(6,'(A,2I3,11D10.4)')'OT',K,a,(G(r),r=IJS-noab(isp)
     !!$,IJS+nuab(isp)-1)
c
               KADT=KADT+NUAB(ISP)
               IJS=IJS+N
            enddo               !A
            enddo               !J
            enddo               !I
         !!enddo                  ! K
       !!write(6,'(A,2I5,4x,5D15.10)')
     !!$'block-m:a1,IAS,ddot',a1,ias,
     !!$ddot_(N*adim*nnoab(ISP),G(IX),1,g(ix),1),(G(I),I=IX,IX+3)
cmp         call multi_wridir(G(IX),n*adim*nnoab(isp),LU,IAS, last)
cmp
cmpn        do jjj=0,n*adim*nnoab(isp)-1
cmpn        if (abs(Work(ix+jjj)).gt.10000) then
cmpn          write (*,*) 'fucko 2'
cmpn          write (*,*) jjj,Work(ix+jjj)
cmpn          stop
cmpn        end if
cmpn        end do
cmp
         call multi_wridir(Work(IX),n*adim*nnoab(isp),LU,IAS, last)
         IAS=IAS+iasblock
c
cmp
        call GetMem('it4_exp','Free','Real',it_exp,length)
cmp
      enddo                     ! A1
cmp        write (6,*) 'cast 1 ok'
cmpn

      IS2=3-ISP
C  lmat
cmp! Mozes odjbt ix, iscr, ig
cmp        write (6,*) 'test 2 na iscr ',vblock*noab(1)*noab(1)
cmp        write (6,*) 'test 2 na ig ',noab(1)*nuab(1)*nno
cmp        write (6,*) 'test 2 na ix ',noab(1)*noab(1)*vblock*n
cmp        write (6,*) 'nno (2) = ',nno
        call GetMem('c2_iscr','Free','Real',iscr,
     & vblock*noab(1)*noab(1))
        call GetMem('c2_ig','Free','Real',ig,
     & noab(1)*nuab(1)*nno)
       call GetMem('c2_ix','Free','Real',ix,
     & noab(1)*noab(1)*vblock*n)
cmp        write (*,*) 'n (2) = ',n
cmp      call w_memchk('IX klvaa ')
cmp      call w_free(g(ix),0,'klvaaix')
cmp      call w_alloc(ix,nnoab(isp)*vblock*vblock,'Ix klvaa-v')
cmp
       call GetMem('klv_oo_ix','Allo','Real',ix,
     & nnoab(isp)*vblock*vblock)
cmp
C starts <aa||oo> integrals
c      FN(1:5)='VVOOI'
      !!FN(6:6)=ich(isp)
c      FN(6:6)=ich(3)
cmp!      CALL GET3DM(FN,G(it),NNUAB(3),NNOAB(3),0)
cmp
cmp!        call w_alloc(il1,nc*no*maxdim,'IL1 klvaa-v')
cmpn       call GetMem('klv_oo_il1','Allo','Real',il1,
cmpn     & nc*no*maxdim)
cmp!        call w_alloc(itmp,
cmp!     & max(nc*no*maxdim,maxdim*maxdim*no*no),'IL2 klvaa-v')
cmpn       call GetMem('klv_oo_itmp','Allo','Real',itmp,
cmpn     & max(nc*no*maxdim,maxdim*maxdim*no*no))
cmp!        call w_alloc(il2,
cmp!     & max(nc*no*maxdim,maxdim*maxdim*no*no),'IL2 klvaa-v')
cmpn       call GetMem('klv_oo_il2','Allo','Real',il2,
cmpn     & max(nc*no*maxdim,maxdim*maxdim*no*no))
c
cmp!        call gen_vvoo(G(it),G(il1),G(itmp),G(il2))
cmpn        call gen_vvoo(Work(it),Work(il1),Work(itmp),Work(il2))

cmp!        open (unit=36,file='vvoo_moje')
cmp!        do i=0,NNUAB(3)*NNOAB(3)-1
cmp!        if (abs(G(it+i)).lt.1.0d-7) G(it+i)=0.0d0
cmp!        write (36,*) i,G(it+i)
cmp!        end do
cmp!        close (36)
c
cmp!        call w_free(G(il1),0,'IL1 klvaa-v')
cmpn       call GetMem('klv_oo_il2','Free','Real',il2,
cmpn     & max(nc*no*maxdim,maxdim*maxdim*no*no))
cmpn       call GetMem('klv_oo_itmp','Free','Real',itmp,
cmpn     & max(nc*no*maxdim,maxdim*maxdim*no*no))
cmpn       call GetMem('klv_oo_il1','Free','Real',il1,
cmpn     & nc*no*maxdim)

cmp        write(6,*)ddot_(nnoab(3)*nnuab(3),G(it),1,G(it),1)

cmp
      !!call dscal_(NNUAB(3)*NNOAB(3),-1.d0,G(it),1)
      !!call dscal_(NNUAB(isp)*NNOAB(ISP),-1.d0,G(it),1)
C
C number of blocks written in a single multiwrite
C
      iasblock=vblock*vblock*nnoab(isp)/nblock
      if((iasblock*nblock).lt.(vblock*vblock*nnoab(isp)))
     $     iasblock=iasblock+1
!!      write(6,*)'create_aa vvoo   iasblock',iasblock
C
!!      FN(1:4)='VMAT'
!!      FN(6:6)=ich(isp)
!!      FN(5:5)=ich(isp)
!!      call multi_opendir(FN,LU)
C currently using 3-dim (big field) - will be replaced after changing
C stepiv and the rest
      do nga=1,nug
         A1=(nga-1)*vblock+1
         adim=min(vblock,nuab(isp)-A1+1)
         A2=A1+adim-1
         do ngb=1,nga
            if(nga.eq.ngb)then
               maxdim2=adim*(adim-1)/2
            else
               maxdim2=adim*vblock
            endif
            B1=(ngb-1)*vblock+1
cmpn
cmp        write (6,*)
cmp        write (6,*) '================================='
cmp        write (6,*) ' nga, ngb',nga,ngb
cmp        write (6,*) '================================='
cmp        write (6,*)
c
cmpn
c - check the largest b2
c
cmp        b2_chk=b1-1+min(vblock,a2-b1)
cmpn        b2_chk=b1-1+vblock
        b2_chk=b1+min(vblock,nuab(isp)-b1+1)-1
c
c - find out which T2 blocked files will be needed
c   for particular nga, ngb
c
cmp        write (6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk
c
        call block_interf(a1,a2,b1,b2_chk,
     &   ngaf,ngal,nind_ngaf,nind_ngal,
     &   ngbf,ngbl,nind_ngbf,nind_ngbl)
c
cmp        write (6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',
cmp     & ngaf,ngal,nind_ngaf,nind_ngal
cmp        write (6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',
cmp     & ngbf,ngbl,nind_ngbf,nind_ngbl
c
c - read amplitudes from T2_ngaf_ngbf ...  T2_ngaf_ngbl, nga>=ngb
c                           ....               ....
c                        T2ngal_ngbf  ...  T2_ngal_ngbl, nga>=ngb
c
c - calculate memory requirements (consider squared T2(a',a'))
c
          length1=0
        do i_blk=ngaf,ngal
          length1=length1+DimGrpaR(i_blk)
        end do
c
          length2=0
        do j_blk=ngbf,ngbl
          length2=length2+DimGrpaR(j_blk)
        end do
c
cmp        write (6,*) 'length1, vblock = ',length1,vblock
cmp        write (6,*) 'length2, vblock = ',length2,vblock
c
        length=length1*length2*no*no
cmp        write (6,*) 'length for blocked VVOO integrals = ',length
c
c - setup memory
c
cmp        write (6,*) 'allocating t2_exp = ',length
        call GetMem('it5_exp','Allo','Real',it_exp,length)
c
c - read pertinent files and generate block of vvoo integrals
c
       call GetMem('vvooil1','Allo','Real',il1,
     & nc*no*maxdim)
       call GetMem('vvooitmp','Allo','Real',itmp,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
       call GetMem('vvooil2','Allo','Real',il2,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
cmp
        call gen_vvoo_blocked(Work(it_exp),Work(il1),
     & Work(itmp),Work(il2),
     & length1,length2,ngaf,ngal,ngbf,ngbl)
cmp
       call GetMem('vvooil2','Free','Real',il2,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
       call GetMem('vvooitmp','Free','Real',itmp,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
       call GetMem('vvooil1','Free','Real',il1,
     & nc*no*maxdim)
cmpn
            do a=a1,a2
               B2=B1-1+min(vblock,A-B1)
               NSTEP=B2-B1+1
               if(nstep.ne.0)then
                  if(nga.eq.ngb)then
                     IJS=(a-a1-1)*(a-a1)/2+IX
                  else
                     IJS=(a-a1)*vblock+IX
                  endif
                  R=0
                  do I=2,NOAB(ISP)
                  do J=1,I-1
                   !!R=R+1
                  !!do R=1,NNOAB(ISP)
                     !!KADT=(R-1)*NNUAB(ISP)
cmpn                  R=(J-1)*noab(isp)+I
cmpn                     KADT=(R-1)*NNUAB(3)
cmpn                     KADT=KADT+(a-1)*NUAB(ISP) +B1 +IT -1
c
        a_tmp=a-nind_ngaf
        b1_tmp=b1-nind_ngbf
c
cmp        write (6,'(A,5(i4,x))') 'a,b1,I,J,nstep ',a,b1,I,J,nstep
cmp        write (6,'(A,5(i4,x))') 'a_tmp,b1_tmp   ',a_tmp,b1_tmp
        do j_tmp=0,nstep-1
c
          KADT=(i-1)*noab(isp)*length1*length2+
     & (j-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+
     & a_tmp+it_exp-1
c
        Work(IJS+j_tmp)=Work(KADT)
c
        end do
c
C address and block for the A1-A2x B1-B2
                     !!KADT=KADT+(a-1)*(a-2)/2 +B1 +IT -1
C    VO >>> G(IX)
cmpn                     call dcopy_(NSTEP,Work(KADT),1,Work(IJS),1)
cmpn                  R=(I-1)*noab(isp)+J
cmpn                     KADT=(R-1)*NNUAB(3)
cmpn                     KADT=KADT+(a-1)*NUAB(ISP) +B1 +IT -1
        do j_tmp=0,nstep-1
c
        KADT=(j-1)*noab(isp)*length1*length2+
     & (i-1)*length1*length2+(b1_tmp+j_tmp-1)*length1+
     & a_tmp+it_exp-1
c
        Work(IJS+j_tmp)=Work(IJS+j_tmp)-1.0d0*Work(KADT)
c
        end do
cmpn
cmpn                     call daxpy_(NSTEP,-1.d0,Work(KADT),1,Work(IJS),1)
       !!write(6,'(A,2I3,8D15.8)')'T',(I-1)*(I-2)/2+j
     !!$,a,(G(r),r=IJS,IJS+NSTEP-1)
                     IJS=IJS+maxdim2
                  enddo         ! RI
                  enddo         ! RJ
               endif
            enddo               ! a1
      !!write(6,'(A,4I5,4x,5D15.10)')
     !!$'block-m: a1,b1,IAS,ddot',a1,b1,ias,maxdim2,
     !!$ ddot_(nnoab(isp)*maxdim2,G(IX),1,G(IX),1),(G(I),I=IX,IX+3)
            IF(maxdim2.eq.0)then
               maxdim2=1
            endif
cmp            call multi_wridir(G(IX),nnoab(isp)*maxdim2,LU,IAS, last)
cmp
cmpn        do jjj=0,nnoab(isp)*maxdim2-1
cmpn        if (abs(Work(ix+jjj)).gt.10000) then
cmpn          write (*,*) 'fucko 3'
cmpn          write (*,*) jjj,Work(ix+jjj)
cmpn          stop
cmpn        end if
cmpn        end do
cmp
            call multi_wridir(Work(IX),nnoab(isp)*maxdim2,LU,IAS, last)
            ias=ias+iasblock
c
cmp        write (6,*) 'deallocating t2_exp = ',length
        call GetMem('it5_exp','Free','Real',it_exp,length)
c
         enddo                  ! ngb
      enddo                     ! nga
cmp        write (6,*) 'cast 2 ok'
cmp      call w_memchk('all klvaa ')
cmp      call w_free(g(ix),0,'klvaa ')
cmp
        call GetMem('klv_oo_ix','Free','Real',ix,
     & nnoab(isp)*vblock*vblock)
cmp
      !!call w_free(g(it),0,'klvaa ')
cmp      call w_alloc(ix,nnoab(3)*vblock*vblock,'ix-vvoo')
cmp
        call GetMem('klv_oo_ix','Allo','Real',ix,
     & nnoab(3)*vblock*vblock)
cmp
      !!call w_alloc(ig,nnoab(3)*nnuab(3),'ig-vvoo')
cmpn       ig=it
C   from now on as for uhf
      iasblock=nnoab(3)*vblock*vblock/nblock
      if((iasblock*nblock).lt.(nnoab(3)*(vblock**2)))iasblock=iasblock+1
!!      FN(1:5)='VVOOI'
!!      FN(6:6)=ICH(3)
!!      CALL GET3DM(FN,G(IG),NNUAB(3),NNOAB(3),0)
C   (c>d|AK)

cmpn
        nga=0
cmpn
      DO A1=1,NUAB(ISP),vblock
         A2=A1+min(vblock,nuab(isp)-A1+1)-1
         adim=A2-A1+1
cmpn
        nga=nga+1
        ngb=0
cmpn
         do B1=1,NUAB(IS2),vblock
cmpn
        ngb=ngb+1
cmpn
            NSTEP=min(vblock,nuab(is2)-B1+1)
cmpn
cmpn
cmp        write (6,*)
cmp        write (6,*) '================================='
cmp        write (6,*) ' nga, ngb',nga,ngb
cmp        write (6,*) '================================='
cmp        write (6,*)
c
c - check the largest b2
c
cmp        b2_chk=b1-1+min(vblock,nuab(is2)-B1+1)
cmpn        b2_chk=b1-1+vblock
        b2_chk=b1+min(vblock,nuab(is2)-B1+1)-1
c
c - find out which T2 blocked files will be needed
c   for particular nga, ngb
c
cmp        write (6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk
c
        if (nga.lt.ngb) then
cmp        write (6,*) 'switching nga, ngb',ngb,nga
c
        call block_interf(b1,b2_chk,a1,a2,
     &   ngaf,ngal,nind_ngaf,nind_ngal,
     &   ngbf,ngbl,nind_ngbf,nind_ngbl)
c
        else
c
        call block_interf(a1,a2,b1,b2_chk,
     &   ngaf,ngal,nind_ngaf,nind_ngal,
     &   ngbf,ngbl,nind_ngbf,nind_ngbl)
c
        end if
c
cmp        write (6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',
cmp     & ngaf,ngal,nind_ngaf,nind_ngal
cmp        write (6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',
cmp     & ngbf,ngbl,nind_ngbf,nind_ngbl
c
c - read amplitudes from T2_ngaf_ngbf ...  T2_ngaf_ngbl, nga>=ngb
c                           ....               ....
c                        T2ngal_ngbf  ...  T2_ngal_ngbl, nga>=ngb
c
c - calculate memory requirements (consider squared T2(a',a'))
c
          length1=0
        do i_blk=ngaf,ngal
          length1=length1+DimGrpaR(i_blk)
        end do
c
          length2=0
        do j_blk=ngbf,ngbl
          length2=length2+DimGrpaR(j_blk)
        end do
c
cmp        write (6,*) 'length1, vblock = ',length1,vblock
cmp        write (6,*) 'length2, vblock = ',length2,vblock
c
        length=length1*length2*no*no
cmp        write (6,*) 'length for blocked VVOO integrals = ',length
c
c - setup memory
c
cmp        write (6,*) 'allocating t2_exp = ',length
        call GetMem('it6_exp','Allo','Real',it_exp,length)
c
c - read pertinent files and generate block of vvoo integrals
c
       call GetMem('vvooil1','Allo','Real',il1,
     & nc*no*maxdim)
       call GetMem('vvooitmp','Allo','Real',itmp,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
       call GetMem('vvooil2','Allo','Real',il2,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
cmp
        call gen_vvoo_blocked(Work(it_exp),Work(il1),
     & Work(itmp),Work(il2),
     & length1,length2,ngaf,ngal,ngbf,ngbl)
cmp
       call GetMem('vvooil2','Free','Real',il2,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
       call GetMem('vvooitmp','Free','Real',itmp,
     & max(nc*no*maxdim,maxdim*maxdim*no*no))
       call GetMem('vvooil1','Free','Real',il1,
     & nc*no*maxdim)
cmpn
!!    bdim=NSTEP
C  mv T2(B,A,I,K) >> G(ix)
            KI=0
            do K=1,noab(isp)
               KADT=(K-1)*NNUAB(3)*(NOAB(2)*(2-ISP)+ISP-1)
               DO I=1,NOAB(IS2)
                  KI=KI+1
cmpn                  IADT=(I-1)*NNUAB(3)*(NOAB(2)*(2-IS2)+IS2-1)
cmpn                  BADT=(2-ISP)*B1+(ISP-1)*(B1-1)*NUAB(2)
                  do A=A1,A2
cmpn
        if (nga.lt.ngb) then
          a_tmp=a-nind_ngbf
          b1_tmp=b1-nind_ngaf
        else
          a_tmp=a-nind_ngaf
          b1_tmp=b1-nind_ngbf
        end if
cmpn
cmpn                     AADT=IADT+KADT+BADT+A*(ISP-1)+(2-ISP)
cmpn     $                    *(A-1)*NUAB(2)+IG-1
C  T2 <A,B,J,K> for isp=1 T2 <B,A,K,J> for isp=2
!!
                     RAD=(KI-1)*adim*nstep+(A-A1)*nstep +IX
cmpn                     ISTEP=(ISP-1)*NUAB(2)+2-ISP
cmpn
        if (nga.ge.ngb) then ! nga> ngb

cmp        write (6,'(A,4(i5,2x),3x,i5)') '(I)  a_tmp,b1_tmp,k,i  nstep = ',
cmp     & a_tmp,b1_tmp,k,i,nstep
c T2(B1,A,I,K) =? T2(A,B1,K,I)
          do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim
c
          AADT=(I-1)*length1*length2*NOAB(2)+(K-1)*length1*length2+
     &   (b1_tmp-1+j_tmp)*length1+a_tmp+it_exp-1
c
          RAD_tmp=RAD+j_tmp
c
          Work(RAD_tmp)=Work(AADT)
c
          end do
c
        else ! nga < ngb
c
cmp        write (6,'(A,4(i5,2x),3x,i5)') '(II) b1_tmp,a_tmp,k,i   nstep = ',
cmp     & b1_tmp,a_tmp,i,k,nstep
c T2(B1,A,I,K)
          do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim
c
          AADT=(k-1)*length1*length2*NOAB(2)+(i-1)*length1*length2+
     &   (a_tmp-1)*length1+b1_tmp+j_tmp+it_exp-1
c
          RAD_tmp=RAD+j_tmp
c
          Work(RAD_tmp)=Work(AADT)
c
          end do
        end if
cmpn
cmpn                     call dcopy_(NSTEP,Work(AADT),ISTEP,Work(RAD),1)

                  enddo         ! A
               enddo            ! I
            enddo               ! K
      !!write(6,'(A,3I5,4x,5D15.10)')
     !!$'block-v: a1,b1,IAS,ddot',a1/vblock+1,b1/vblock+1,ias,
     !!$ ddot_(NNOAB(3)*adim*nstep,G(IX),1,G(IX),1),(G(I),I=IX,IX+3)
cmp            call multi_wridir(G(IX),NNOAB(3)*adim*nstep,LU,IAS, last)
cmp
cmpn        do jjj=0,NNOAB(3)*adim*nstep-1
cmpn        if (abs(Work(ix+jjj)).gt.10000) then
cmpn          write (*,*) 'fucko 1'
cmpn          write (*,*) jjj,Work(ix+jjj)
cmpn          stop
cmpn        end if
cmpn        end do
cmp
            call multi_wridir(Work(IX),NNOAB(3)*adim*nstep,LU,IAS, last)
            ias=ias+iasblock
cmpn
        call GetMem('it6_exp','Free','Real',it_exp,length)
cmpn
         enddo                  ! B1
      enddo                     ! A1
cmp        write (6,*) 'cast 3 ok'
      !!close (LU)   in calling routine
!     write(6,*) FN, isp, IAS
      !!dupblk(ndup)=last   in calling routine
cmp      call w_memchk('all klvaa ')
      !!call w_free(g(ix),0,'IT klvaa ')   in calling routine
cmp
        call GetMem('klv_oo_ix','Free','Real',ix,
     & nnoab(3)*vblock*vblock)
cmp
cmpn
        if (printkey.gt.1) then
           write (6,*) 'VVOO integrals regenerated from MOLCAS'
        end if

cmpn
      call xflush(6)
      !!stop
c
      return
      end
