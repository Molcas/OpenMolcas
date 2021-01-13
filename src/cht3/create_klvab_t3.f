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
cmp      SUBROUTINE create_klvab_t3(G,vblock)
      SUBROUTINE create_klvab_t3(vblock)
C
C  creates K(alpha-beta,alpha-beta),K(beta-alpha,alpha-beta)
C  DA files KMATBA and KMATAB, LMATBA and LMATAB
C  creates L(alpha-beta,alpha-beta),K(beta-alpha,alpha-beta)
C  max G at this place
C  both matrices have to be available at a time (unfortunately)
C
C  parallelization irrelevant at the moment
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
        integer it_exp,RAD_tmp
        integer a_tmp,b1_tmp,j_tmp
        integer nga,ngb
cmp
cmpn        integer AADT_tmp,it
c       integer jjj
cmp
cmpn
cmp      real*8 G(*),ddot_
c     real*8 ddot_
cmp      integer it,ix,ig,iscr, KADT, IJS, RAD, AADT, IADR
      integer ix,ig,iscr, IJS, RAD, AADT, IADR
      integer isp,is2,ias,vblock,n,i,j,k,lu,iasblock,ias_aa
      INTEGER A,A1,A2,B1,NSTEP,ISTEP
      CHARACTER FN*6,ich*1
      INTEGER IOPT,NOAB,NNOAB,NUAB,NNUAB,NNU,IUHF,NNO,ISPA
      integer adim, last,last_aa,nug
c     integer bdim
      COMMON/UHF/NOAB(2),NNOAB(3),NUAB(2),NNUAB(3),ICH(3)
      COMMON/IOIND/IOPT(96)
cmp
        integer itmp,il1_1,il2_1,il0,il1,it2_tmp,itmp2
        logical switch
cmp
        if (printkey.ge.10) then
        write (6,*)
        write (6,*) '------ DimGrpaR ------'
        write (6,'(8(i5,2x))') (DimGrpaR(a_tmp),a_tmp=1,NvGrp)
        write (6,*)
        end if
c
c - calculate overall memory requirements
c
        call check_create_klvab_t3_mem (vblock)
c
cmp      call w_rescope(G,'G create KL')
cmp      call w_free(g,0,'G klvab ')
      N=noab(1)+nuab(1)
c      NNRED=NOAB(1)*(NOAB(1)+1)/2
      IUHF=0
      IF(IOPT(76).NE.0)THEN
         IUHF=1
c         NNRED=NNOAB(3)
      ENDIF
      LU=98

cmp      call w_alloc(it,NNOAB(3)*NNUAB(3),'IT klvab ')
cmp      call w_alloc(ix,vblock*vblock*n,'IX klvab ')
cmpn        call GetMem('create_it','Allo','Real',it,
cmpn     & NNOAB(3)*NNUAB(3))
        call GetMem('c1_ix','Allo','Real',ix,
     & vblock*vblock*n)
cmp!
c!      FN(1:6)='T2OLDC'
c!      CALL GET3DM(FN,G(it),NNUAB(3),NNRED,0)
cmp!

cmp        call w_alloc(it2_tmp,maxdim*maxdim*no*no,'ITMP klvab ')
cmp        call w_alloc(itmp,maxdim*maxdim*no*no,'ITMP klvab ')
cmpn        call GetMem('create_it2_tmp','Allo','Real',it2_tmp,
cmpn     & maxdim*maxdim*no*no)
cmpn        call GetMem('create_itmp','Allo','Real',itmp,
cmpn     & maxdim*maxdim*no*no)
c
cmp        call gather_t2(G(it),G(it2_tmp),G(itmp))
cmpn        call gather_t2(Work(it),Work(it2_tmp),Work(itmp))
c
c
cmp        call w_free(G(it2_tmp),0,'ITMP klvab ')
cmpn        call GetMem('create_itmp','Free','Real',itmp,
cmpn     & maxdim*maxdim*no*no)
cmpn        call GetMem('create_it2_tmp','Free','Real',it2_tmp,
cmpn     & maxdim*maxdim*no*no)

cmpn        write (6,*) 'T2 regenerated from MOLCAS'
cmpn        write (6,*)
c
cmp      call dscal_(NNUAB(3)*NNRED,-1.d0,G(it),1)
cmp         call dscal_(NNUAB(3)*NNOAB(3),-1.d0,G(it),1)
cmpn         call dscal_(NNUAB(3)*NNOAB(3),-1.d0,Work(it),1)
cmp      IF(IUHF.EQ.0) call decomp2ind(G(it),NNUAB(3),noab(1),NUAB(1))
       !!write(6,*)ddot_(nnoab(3)*nnuab(3),G(it),1,G(it),1)
c
C number of blocks written in a single multiwrite
C
      iasblock=vblock*vblock*N/nblock
      if((iasblock*nblock).lt.(vblock*vblock*N))iasblock=iasblock+1
C
      do isp=1,IUHF+1
         is2=3-isp
         FN(1:4)='KMAT'
         FN(6:6)=ich(isp)
         FN(5:5)=ich(3-isp)
         Write (6,*) 'FN,LU=',FN,LU
         call multi_opendir(FN,LU)
         ndup=ndup+1
         if (ndup.gt.ndupmx)
     $        call barf('create_klvab_t3 -- ndupmx exceeded')
         !!write(6,*) FN, isp,ndup
         dupfil(ndup)=FN
         if(IUHF.eq.0)then
         FN(6:6)=ich(isp)
         FN(5:5)=ich(isp)
         call multi_opendir(FN,LU+1)
         !!ndup=ndup+1
         if (ndup.gt.ndupmx)
     $        call barf('create_klvab_t3 -- ndupmx exceeded')
         !!write(6,*) FN, isp,ndup+1
         dupfil(ndup+1)=FN
         endif
C currently using 3-dim (big field) - will be replaced after changing
C stepiv and the rest
         nnu=(nuab(is2)*(nuab(is2)+1))/2
cmp         call w_alloc(ig,(nuab(isp)*nnu),'IG klvab')
cmp         call w_alloc(iscr,nuab(is2)*nuab(is2),'IG iscr')
cmp
        call GetMem('c1_ig','Allo','Real',ig,
     & nuab(isp)*nnu)
        call GetMem('c1_iscr','Allo','Real',iscr,
     & nuab(is2)*nuab(is2))
cmp
         IAS=1
         IAS_AA=1
c
         DO K=1,noab(isp)
!!            IF(IUHF.EQ.1)THEN
!!               KADT=(K-1)*NNUAB(3)*(NOAB(2)*(2-ISP)+ISP-1)
!!            ELSE
!!               KADT=0
!!            ENDIF

cmp!            FN(1:5)='VVVAI'
cmp!            FN(6:6)=ICH(ISP)
cmp!            CALL GET3DM(FN,G(IG),NNU,NUAB(ISP),K)
        if (printkey.gt.1) then
        write (6,*) 'Regenerating VVVo integrals for o = ',K
        end if
cmp
cmp         call w_alloc(il1_1,nc*maxdim,'IL1_1 iscr')
cmp         call w_alloc(il2_1,nc*maxdim*maxdim,'IL2_1 iscr')
cmp         call w_alloc(itmp,
cmp     & max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim),
cmp     & 'ITMP iscr')
cmp
        call GetMem('cc_il1_1','Allo','Real',il1_1,
     & nc*maxdim)
        call GetMem('cc_il2_1','Allo','Real',il2_1,
     & nc*maxdim*maxdim)
        call GetMem('cc_itmp','Allo','Real',itmp,
     & max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim))
cmp
cmp        call gen_vvvo(K,G(IG),G(il1_1),G(il2_1),G(itmp))
        call gen_vvvo(K,Work(IG),
     & Work(il1_1),Work(il2_1),Work(itmp))
cmp       write(6,*) ddot_(nnu*nuab(isp),G(ig),1,G(ig),1)
cmp        call zeroma (G(ig),1,NNU*NUAB(ISP))
c
cmp        call w_free(G(il1_1),0,'IL1_1 iscr')
        call GetMem('cc_itmp','Free','Real',itmp,
     & max(nc*maxdim*maxdim,nc*no*maxdim,maxdim*maxdim*maxdim))
        call GetMem('cc_il2_1','Free','Real',il2_1,
     & nc*maxdim*maxdim)
        call GetMem('cc_il1_1','Free','Real',il1_1,
     & nc*maxdim)
cmp
            call delf(FN,K,K)
       if(iuhf.eq.0)
     $call klvaa_vvv(ix,ig,iscr,vblock,N,nug,
     $LU+1,last_aa,iasblock,K,ias_aa)
cmpn     $call klvaa_vvv(ix,it,ig,iscr,vblock,N,nug,
cmpn     $LU+1,last_aa,iasblock,K,ias_aa)
cmp     $call klvaa_vvv(G,ix,it,ig,iscr,vblock,N,nug,
cmp     $LU+1,last_aa,iasblock,K,ias_aa)
c

!!       call xflush(6)
C   (c>d|AK)

cmpn
        nga=0
cmpn
            DO A1=1,NUAB(ISP),vblock
               A2=A1+min(vblock,nuab(isp)-A1+1)-1
               adim=A2-A1+1
        nga=nga+1
        ngb=0
               do B1=1,NUAB(IS2),vblock
        ngb=ngb+1
                  NSTEP=min(vblock,nuab(is2)-B1+1)
!!    bdim=NSTEP
                  IJS=(A1-1)*NNU+IG
cmpn
cmpn        write (6,*)
cmpn        write (6,*) '================================='
cmpn        write (6,*) ' nga, ngb',nga,ngb
cmpn        write (6,*) '================================='
cmpn        write (6,*)
c
c - check the largest b2
c
cmpn        b2_chk=b1-1+vblock
        b2_chk=b1-1+min(vblock,nuab(is2)-B1+1)
c
c - find out which T2 blocked files will be needed
c   for particular nga, ngb
c
cmpn        write (6,'(A,4(i5,2x))') 'a1,a2,b1,b2_chk = ',a1,a2,b1,b2_chk
c
        switch=.false.
        if (nga.lt.ngb) then
cmpn        write (6,*) 'switching nga, ngb',ngb,nga
        switch=.true.
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
cmpn        write (6,'(A,4(i5,2x))') 'ngaf, ngal, nind_ngaf, nind_ngal',
cmpn     & ngaf,ngal,nind_ngaf,nind_ngal
cmpn        write (6,'(A,4(i5,2x))') 'ngbf, ngbl, nind_ngbf, nind_ngbl',
cmpn     & ngbf,ngbl,nind_ngbf,nind_ngbl
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
cmpn        write (6,*) 'length1, vblock = ',length1,vblock
cmpn        write (6,*) 'length2, vblock = ',length2,vblock
c
        length=length1*length2*no*no
cmpn        write (6,*) 'length for blocked T2 amplitudes = ',length
c
c - setup memory
c
cmpn        write (6,*) 'allocating t2_exp = ',length
        call GetMem('it2_exp','Allo','Real',it_exp,length)
c
c - read pertinent files and store them in the new blocked structure
c
        call GetMem('cd_it2tmp','Allo','Real',it2_tmp,
     & maxdim*maxdim*no*no)
        call GetMem('cd_itmp','Allo','Real',itmp,
     & maxdim*maxdim*no*no)
c
        call gather_t2_blocked(length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & Work(it_exp),Work(it2_tmp),Work(itmp),
     & switch)
c
        call GetMem('cd_itmp','Free','Real',itmp,
     & maxdim*maxdim*no*no)
        call GetMem('cd_it2tmp','Free','Real',it2_tmp,
     & maxdim*maxdim*no*no)

        call dscal_(length,-1.d0,Work(it_exp),1)
c
                  do A=A1,A2
cmp                     CALL EXPA1_UHF(G(IJS),1,NUAB(IS2),1,G(ISCR))
                     CALL EXPA1_UHF(Work(IJS),1,NUAB(IS2),1,Work(ISCR))
                     IJS=IJS+NNU
C  Gix       <ng*ng|R,k>
!! not needed      B2=B1+min(ng,nuab(is2)-B1+1))-1
cmpn
        if (nga.lt.ngb) then
          a_tmp=a-nind_ngbf
          b1_tmp=b1-nind_ngaf
        else
          a_tmp=a-nind_ngaf
          b1_tmp=b1-nind_ngbf
        end if
cmp        write (6,'(A,2(i5,2x),A,i5,2x)') 'b1,   a     ,i,k = ',b1,a,
cmp     & '    i',k
C  mv T2(B,A,I,K) >> G(ix)
cmpn
                     DO I=1,NOAB(IS2)
                        ISPA=ISP
                        !!IF(IUHF.EQ.1)THEN
                        !!   IADT=(I-1)*NNUAB(3)*(NOAB(2)*(2-IS2)+IS2-1)
                        !!ELSE
                        !!   IADT=(MAX(K,I)-1)*(MAX(K,I))/2+MIN(K,I)
                        !!   IADT=(IADT-1)*NNUAB(3)
                        !!   IF(K.LT.I)ISPA=IS2
                        !!ENDIF
cmpn                        BADT=(2-ISPA)*B1+(ISPA-1)*(B1-1)*NUAB(2)
cmpn!!!
cmpn                        AADT_tmp=IADT+KADT+BADT+A*(ISPA-1)
cmpn     $                       +(2-ISPA)*(A-1)*NUAB(2)+IT-1
cmpn!!!
C  T2 <A,B,J,K> for isp=1 T2 <B,A,K,J> for isp=2
!!      RAD=(I-1)*vblock*vblock+(A-A1)*vblock+IX
                        RAD=(I-1)*adim*nstep+(A-A1)*nstep +IX
                        ISTEP=(ISPA-1)*NUAB(2)+2-ISPA
cmpn
        if (nga.ge.ngb) then ! nga> ngb

cmp        write (6,'(A,4(i5,2x),3x,i5)') '(I)  a_tmp,b1_tmp,k,i  nstep = ',
cmp     & a_tmp,b1_tmp,k,i,nstep
c T2(B,A,I,K) =? T2(A,B1,K,I)
          do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim
c
          AADT=(I-1)*length1*length2*NOAB(2)+(K-1)*length1*length2+
     &   (b1_tmp-1+j_tmp)*length1+a_tmp+it_exp-1
c
          RAD_tmp=RAD+j_tmp
c
          Work(RAD_tmp)=Work(AADT)
cmp
cmpn        if (abs(Work(AADT)-Work(AADT_tmp+j_tmp)).
cmpn     & gt.0.00001d0) then
cmpn           write (*,*) 'halohaha 1',AADT,AADT_tmp+j_tmp,
cmpn     & Work(AADT),Work(AADT_tmp+j_tmp)
cmpn           stop
cmpn        end if
cmp
          end do
c
        else ! nga < ngb
c
cmp        write (6,'(A,4(i5,2x),3x,i5)') '(II) b1_tmp,a_tmp,k,i   nstep = ',
cmp     & b1_tmp,a_tmp,i,k,nstep
c T2(B,A,I,K)
          do j_tmp=0,NSTEP-1  ! istep je 1 ak dobre tusim
c
cmp          AADT=(k-1)*length1*length2*NOAB(2)+(i-1)*length1*length2+
cmp!          AADT=(i-1)*length1*length2*NOAB(2)+(k-1)*length1*length2+
          AADT=(K-1)*length1*length2*NOAB(2)+(I-1)*length1*length2+
     &   (a_tmp-1)*length1+b1_tmp+j_tmp+it_exp-1
c
          RAD_tmp=RAD+j_tmp
c
          Work(RAD_tmp)=Work(AADT)
cmp
cmpn        if (abs(Work(AADT)-Work(AADT_tmp+j_tmp)).
cmpn     & gt.0.00001d0) then
cmpn           write (*,*) 'halohaha 2',AADT,AADT_tmp+j_tmp,
cmpn     & Work(AADT),Work(AADT_tmp+j_tmp)
cmpn           stop
cmpn        end if
cmp
c
          end do
        end if
c
cmpn                        call dcopy_(NSTEP,Work(AADT),ISTEP,Work(RAD),1)
cmpn
                     enddo      ! I
                     RAD=noab(is2)*adim*nstep+(A-A1)*nstep+IX
                     DO IADR=ISCR+B1-1,ISCR+NUAB(IS2)*NUAB(IS2)-1,
     $                    NUAB(IS2)
cmp                        call dcopy_(NSTEP,G(IADR),1,G(RAD),1)
                        call dcopy_(NSTEP,Work(IADR),1,Work(RAD),1)
                        RAD = RAD+adim*nstep
                     enddo      ! IADR
                  enddo         ! A
      !!write(6,'(A,4I5,4x,D15.10)')
     !!$'block-w: K,a1,b1,IAS,ddot',K,a1,b1,ias,
     !!$ ddot_(N*vblock*vblock,G(IX),1,G(IX),1)
cmp                  call multi_wridir(G(IX),N*vblock*vblock,LU,IAS,last)

cmp
cmp!!        do jjj=0,N*vblock*vblock-1
cmp!!        if (abs(Work(ix+jjj)).gt.10000) then
cmp!!          write (6,*) 'prasa 1 ',jjj,Work(ix+jjj)
cmp!!          stop
cmp!!        end if
cmp!!        end do
cmp
                  call multi_wridir(Work(IX),N*vblock*vblock,
     & LU,IAS,last)
        !!write (6,*) 'N*vblock*vblock,LU,last ',
     !!& N*vblock*vblock,LU,last
                  ias=ias+iasblock
cmp
        call GetMem('it2_exp','Free','Real',it_exp,length)
cmp
               enddo            ! B1
            enddo               ! A1
         enddo                  ! K
        if (printkey.gt.1) then
        write (6,*) 'VVVo integrals regenerated from MOLCAS'
        write (6,*)
        end if
cmp
         close (LU)
         dupblk(ndup)=last
         if(IUHF.EQ.0)then
         close(LU+1)
         ndup=ndup+1
         dupblk(ndup)=last_aa
         endif
!         write(6,*) FN, isp, IAS
cmp         call w_memchk('IG klvab ')
cmp         call w_free(g(ig),0,'IG klvab ')
cmp
        call GetMem('c1_iscr','Free','Real',iscr,
     & nuab(is2)*nuab(is2))
        call GetMem('c1_ig','Free','Real',ig,
     & nuab(isp)*nnu)
cmp
      enddo                     ! ISP

cmp      call dscal_(NNUAB(3)*NNOAB(3),-1.d0,G(it),1)
cmpn      call dscal_(NNUAB(3)*NNOAB(3),-1.d0,Work(it),1)
cmp??
        call GetMem('c1_ix','Free','Real',ix,
     & vblock*vblock*n)
cmp??
      do isp=1,IUHF+1
cmp         call w_memchk('IX klvab ')
cmp         call w_free(g(ix),0,'IX klvab ')
cmp
         is2=3-isp
         iasblock=nnoab(3)*vblock*N/nblock
         if((iasblock*nblock).lt.(nnoab(3)*vblock*N))iasblock=iasblock+1
         FN(1:4)='LMAT'
         FN(5:5)=ich(3-isp)
         FN(6:6)=ich(isp)
         call multi_opendir(FN,LU)
         ndup=ndup+1
         if (ndup.gt.ndupmx)
     $        call barf('create_klvab_t3 -- ndupmx exceeded')
!         write(6,*) FN, isp,ndup
         dupfil(ndup)=FN
C
         FN(1:5)='OOVAI'
         FN(6:6)=ICH(ISP)
cmp         call w_alloc(ix,noab(isp)*noab(IS2)*vblock*n,'IX klvabo ')
        call GetMem('c2_ix','Allo','Real',ix,
     & noab(isp)*noab(IS2)*vblock*n)
         nno=noab(is2)*(noab(is2)+1)/2
cmp         call w_alloc(ig,noab(isp)*nuab(isp)*nno,'IG klvabo')
cmp         call w_alloc(iscr,vblock*noab(IS2)*noab(IS2),'ISCRo klvabo ')
        call GetMem('c2_ig','Allo','Real',ig,
     & noab(isp)*nuab(isp)*nno)
        call GetMem('c2_iscr','Allo','Real',iscr,
     & vblock*noab(IS2)*noab(IS2))
cmp
cmp!         CALL GET3DM(FN,G(ig),NNO,NUAB(ISP)*NOAB(ISP),0)
cmp
cmp          call w_alloc(il0,nc*nno,'IL0 klvabo')
cmp          call w_alloc(il1,nc*no*nv,'IL1 klvabo')
cmp          call w_alloc(itmp,
cmp     & max(nc*nno,nc*no*maxdim,nc*no*nv),'ITMP klvabo')
cmp
        call GetMem('cr_il0','Allo','Real',il0,
     & nc*nno)
        call GetMem('cr_il1','Allo','Real',il1,
     & nc*no*nv)
        call GetMem('cr_itmp','Allo','Real',itmp,
     & max(nc*nno,nc*no*maxdim,nc*no*nv))
cmp
        call gen_oovo (Work(ig),Work(il0),Work(il1),Work(itmp))
cmp        call gen_oovo (G(ig),G(il0),G(il1),G(itmp))
cmp        write(6,*)ddot_(NNO*NUAB(ISP)*NOAB(ISP),G(ig),1,G(ig),1)
cmp        call zeroma (G(ig),1,NNO*NUAB(ISP)*NOAB(ISP))
c
cmp        call w_free(G(il0),0,'IL0 klvab')
        call GetMem('cr_itmp','Free','Real',itmp,
     & max(nc*nno,nc*no*maxdim,nc*no*nv))
        call GetMem('cr_il1','Free','Real',il1,
     & nc*no*nv)
        call GetMem('cr_il0','Free','Real',il0,
     & nc*nno)

        if (printkey.gt.1) then
         write (6,*) 'OOVO integrals regenerated from MOLCAS'
        end if
cmp
         IAS=1
cmpn
        nga=0
cmpn
         do A1=1,NUAB(ISP),vblock
            A2=A1+min(vblock,nuab(isp)-(A1-1))-1
            NSTEP=min(vblock,nuab(isp)-(A1-1))
cmpn
        nga=nga+1
c
cmp        write (6,*)
cmp        write (6,*) '================================='
cmp        write (6,*) ' nga ',nga
cmp        write (6,*) '================================='
cmp        write (6,*)
c
cmp        write (6,'(A,2(i5,2x))') 'b1,b2 = ',a1,a2
c
cmp        call block_interf(1,1,a1,a2,
        call block_interf(1,nuab(1),a1,a2,
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
        length=length1*length2*no*no
cmp        write (6,*) 'length for blocked T2 amplitudes = ',length
c
c - setup memory
c
cmp        write (6,*) 'allocating t2_exp = ',length
        call GetMem('it3_exp','Allo','Real',it_exp,length)
c
c - read pertinent files and store them in the new blocked structure
c
        call GetMem('it2_tmp','Allo','Real',it2_tmp,
     & maxdim*maxdim*no*no)
        call GetMem('itmp2','Allo','Real',itmp2,
     & maxdim*maxdim*no*no)
c
        call gather_t2_fblocked(length1,length2,
     & ngaf,ngal,ngbf,ngbl,
     & Work(it_exp),Work(it2_tmp),Work(itmp2))
c
        call GetMem('itmp2','Free','Real',itmp2,
     & maxdim*maxdim*no*no)
        call GetMem('it2_tmp','Free','Real',it2_tmp,
     & maxdim*maxdim*no*no)
c
cmpn
            do k=1,noab(isp)
               !!IF(IUHF.EQ.1) THEN
               !!   KADT=(K-1)*NNUAB(3)*(NOAB(2)*(2-ISP)+ISP-1)
               !!ELSE
               !!   KADT=0
               !!ENDIF
               IJS=ig+(A1-1)*nno+(k-1)*nuab(isp)*nno
cmp               CALL EXPA1_UHF(G(IJS),nstep,NOAB(IS2),1,G(iscr))
               CALL EXPA1_UHF(Work(IJS),nstep,NOAB(IS2),1,Work(iscr))
               do I=1,noab(is2)
                  ISPA=ISP
                  !!IF(IUHF.EQ.1) THEN
                   !!   IADT=(I-1)*NNUAB(3)*(NOAB(2)*(2-IS2)+IS2-1)
                  !!ELSE
                   !!  IADT=(MAX(K,I)-1)*(MAX(K,I))/2+MIN(K,I)
                   !!  IADT=(IADT-1)*NNUAB(3)
                   !!  IF(K.LT.I)ISPA=IS2
                  !!ENDIF
c
                  DO A=A1,A2

cmpn
        a_tmp=a-nind_ngbf
c
        AADT=(K-1)*noab(2)*NUAB(2)*length2+
     &       (i-1)*NUAB(2)*length2+
     &       (a_tmp-1)*NUAB(2)+it_exp
cmpn
C copies T
cmpn                     BADT=2-ISPA+(ISPA-1)*NUAB(2)
C  T2 <A,B,J,K> for isp=1 T2 <B,A,K,J> for isp=2
                     ISTEP=(ISPA-1)*NUAB(2)+2-ISPA
                     RAD=(I-1)*nstep*n+(A-A1)*n+IX+noab(is2)
     $                    +(K-1)*noab(is2)*nstep*N
cmp                     call dcopy_(nuab(is2),G(AADT),ISTEP,G(RAD),1)
cmp
                     call dcopy_(nuab(is2),Work(AADT),ISTEP,Work(RAD),1)
                  enddo         ! A
               enddo            ! I
C copies OOVO
               DO J=1,noab(is2)
                  DO A=1,NSTEP
                     RAD=IX+(A-1)*n+(J-1)*nstep*N+(K-1)*noab(is2)
     $                    *nstep*N
                     IADR=ISCR+(A-1)*noab(is2)*noab(is2)+(j-1)*noab(is2)
cmp                     call dcopy_(noab(is2),G(IADR),1,G(RAD),1)
                     call dcopy_(noab(is2),Work(IADR),1,Work(RAD),1)
                  enddo         ! A
               enddo            ! J
            enddo               ! K
!!      write(6,'(A,2I5,4x,D15.10)')
!!     $'block-w:a1,IAS,ddot',a1,ias,ddot_(N*vblock*nnoab(3),g(ix),1,g(ix),1)
!!      write(6,'(a,a,2I4,D16.8)')'block-w',ich(isp),(A1/vblock)+1,IAS,
!!     $ddot_(N*nnoab(3)*NSTEP,g(ix),1,g(ix),1)
cmp            call multi_wridir(G(IX),N*nstep*nnoab(3),LU,IAS,last)
cmp
cmp!!        do jjj=0,N*nstep*nnoab(3)-1
cmp!!        if (abs(Work(ix+jjj)).gt.10000) then
cmp!!          write (*,*) 'prasa 2 ',jjj,Work(ix+jjj)
cmp!!          stop
cmp!!        end if
cmp!!        end do
cmp
            call multi_wridir(Work(IX),N*nstep*nnoab(3),LU,IAS,last)
cmp
            IAS=IAS+iasblock
cmpn
        call GetMem('it3_exp','Free','Real',it_exp,length)
cmpn
         enddo                  ! A1
         close (LU)
!         write(6,*) FN, isp, IAS
         dupblk(ndup)=last
         if(IUHF.EQ.0)then
         FN(1:4)='LMAT'
         FN(6:6)=ich(isp)
         FN(5:5)=ich(isp)
         call multi_opendir(FN,LU)
         ndup=ndup+1
         if (ndup.gt.ndupmx)
     $        call barf('create_klvab_t3 -- ndupmx exceeded')
!         write(6,*) FN, isp,ndup
         dupfil(ndup)=FN
C   this is to ensure correct copy to slaves
         IAS_AA=1
cmp        write (6,*) 'test 1 na iscr ',vblock*noab(IS2)*noab(IS2)
cmp        write (6,*) 'test 1 na ig ',noab(isp)*nuab(isp)*nno
cmp        write (6,*) 'test 1 na ix ',noab(isp)*noab(IS2)*vblock*n
cmp         call klvaa_oovo(G,ix,it,ig,iscr,vblock,N,nug,
cmpn         call klvaa_oovo(ix,it,ig,iscr,vblock,N,nug,
         call klvaa_oovo(ix,ig,iscr,vblock,N,nug,
     $LU,last_aa,ias_aa)
cmp
cmp        write (6,*) 'klvaa_oovo finished'
cmp
         close(LU)
         ndup=ndup
         dupblk(ndup)=last_aa
         endif
      enddo                     ! ISP
cmp      call w_memchk('all klvab ')
cmp      call w_free(g(it),0,'IT klvab ')
c iscr a ig sa odalokuju v klva_oovo
cmp @@@        call GetMem('create_iscr','Free','Real',iscr,
cmp @@@     & vblock*noab(IS2)*noab(IS2))
cmp @@@        call GetMem('create_ig','Free','Real',ig,
cmp @@@     & noab(isp)*nuab(isp)*nno)
cmp ???
cmp??       call GetMem('c2_ix','Free','Real',ix,
cmp??     & noab(isp)*noab(IS2)*vblock*n)
cmpn        call GetMem('create_it','Free','Real',it,
cmpn     & NNOAB(3)*NNUAB(3))
      call xflush(6)
      return
      end
