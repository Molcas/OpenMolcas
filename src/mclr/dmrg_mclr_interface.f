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

! *********************************************************************
      subroutine read_dmrg_parameter_for_mclr()

#include "Input.fh"
#include "dmrginfo_mclr.fh"

        integer ierr,i
        open(unit=100,file="dmrg_for_mclr.parameters",
     &       status='OLD',action='READ',iostat=ierr)
          if(ierr.ne.0)then
            doDMRG=.false.
            doMCLR=.false.
          else
            read(100,"(11X,L,4X)")   doDMRG
            read(100,"(4X,I8,4X)")nele_RGLR
            read(100,"(4X,I8,4X)") ms2_RGLR
!>          write(6,*)doDMRG,dmrg_state%nactel,dmrg_state%ms2
            do i=1,8
              read(100,"(4X,I3)",advance='no')RGras2(i)
            end do
            read(100,*)
            do i=1,8
              read(100,"(4X,I3)",advance='no')LRras2(i)
            end do
!>          write(6,*)RGras2
!>          write(6,*)LRras2
            read(100,*)
            read(100,"(4X,I8,4X)")nstates_RGLR
!            allocate(checkpoint(nstates_RGLR))
!            checkpoint=""
            do i=1,nstates_RGLR
              read(100,*)
              read(100,'(G20.12)')ERASSCF(i)
              write(6,*)"RASSCF energy", ERASSCF(i)
            end do
!> It is redundant
            doMCLR=.true.
          end if
        close(100)

        write(6,*)"doDMRG, nele_dmrg, ms2_dmrg"
        write(6,*) doDMRG, nele_rglr, ms2_rglr
        call xflush(6)

      end subroutine read_dmrg_parameter_for_mclr

! *********************************************************************
      Subroutine dmrg_spc_change_mclr(orbspc,lrspc)

       integer::orbspc(8)
       integer::lrspc(8)

       lrspc=0

       lrspc=orbspc

      End Subroutine dmrg_spc_change_mclr

! *********************************************************************

      Subroutine dmrg_dim_change_mclr(orbspc,ndim,iflag)

        integer::orbspc(8)
        integer::iflag
        integer::ndim
        integer i,n1,n2

!        write(6,*)"==================================================="
!        write(6,*)" Currently, only valid for no symmetry calculation"
!        write(6,*)"==================================================="

!        I remember it was already valid for all symmetry.
!                         -- yma, need check it again 2015.5.14

!        write(6,*)"orbspc",orbspc(1:8)
        ndim=0

        n1 = 0

        if(iflag.eq.0)then
          do i=1,8
            ndim=ndim+orbspc(i)
          end do
        else if(iflag.eq.1)then
          do i=1,1
            n1=n1+orbspc(i)
            ndim=n1*n1
          end do
        else if(iflag.eq.2)then
          do i=1,1
            n1=n1+orbspc(i)
            ndim=n1**4
          end do
        else if(iflag.eq.3)then
          do i=1,1
            n1=n1+orbspc(i)
            ndim=(n1+1)*n1/2
          end do
        else if(iflag.eq.4)then
          do i=1,1
            n1=n1+orbspc(i)
            n2=n1*n1
            ndim=(n2+1)*n2/2
          end do
        else
          write(6,*)"unknow iflag"
          call Quit_OnUserError()
        end if

      End Subroutine dmrg_dim_change_mclr
! *********************************************************************

      Subroutine ci_reconstruct(istate,nSDET,vector,indexSD)

#include "dmrginfo_mclr.fh"
        character*100,allocatable :: checkpoint(:)  ! for many states

        integer :: istate,nsdet

        integer i,j,idet,idx_det
        integer norb,norbLR

        integer neletol
        integer nele_mod
        integer nele_alpha,nele_beta

        integer irrep_pre,iorbLR0,iorbLR
        integer irrep_diff(8)

! need to be rewrittren using the Molcas getmem,
!                                    or mma_allocate and mma_deallocate
! soon ..
        integer,allocatable::ele_orb_alpha(:)
        integer,allocatable::ele_orb_beta(:)

        integer,allocatable::pre_ele(:)

        integer ndets_mclr

        character*200 tmp_run

        real*8 dtmp

        logical IFFILE
        integer rc

        integer :: indexSD(nsdet)            ! index
        real*8  :: vector(nsdet)             ! determinants

        type Slater_determinant
          integer :: itype                        = 1 ! excitation type
          integer :: inum                         = 1 ! determinant number
          integer :: isign                        = 1 ! determinant phase
          integer,allocatable :: electron(:)          ! determinant
          integer,allocatable :: ele_conf(:)          ! electron configuration
          real*8 ,allocatable :: dV(:)                ! for many states
        end type Slater_determinant

        type (Slater_determinant), allocatable::SD_DMRG(:)

          !> The total electrons
          neletol=0
          neletol= nele_RGLR

          !> Check if there is single electron
          nele_mod=mod(neletol,2)

          !> electrons in alpha or beta orbitals
          if(nele_mod.eq.0)then
            !> If no single electron
            nele_alpha= neletol/2 + ms2_RGLR/2
            nele_beta = neletol/2 - ms2_RGLR/2
          else
            nele_alpha= neletol/2 + ms2_RGLR/2+1
            nele_beta = neletol/2 - ms2_RGLR/2
          end if

! Read in all the name of checkpoint file
          lcheckpoint=20
          lcheckpoint=isFreeUnit(lcheckpoint)
          Call Molcas_Open(lcheckpoint,'dmrg_for_mclr.parameters')
          allocate(checkpoint(nstates_RGLR))
          checkpoint=""
          do i=1,6
            read(lcheckpoint,*)
          end do
          do i=1,nstates_RGLR
            read(lcheckpoint,*)checkpoint(i)
            read(lcheckpoint,*)
          end do
          close(lcheckpoint)

          ! reconstructing dets for current state
          do i=istate,istate
            write(6,*)trim(checkpoint(i))
          end do

          ! preparing for point group symmetry
          norb=0
          norbLR=0
          irrep_diff=0
          do j=1,8
             norb          =  norb      +   RGras2(j)
             norbLR        =  norbLR    +   LRras2(j)
             irrep_diff(j) =  RGras2(j) -   LRras2(j)
          end do

          allocate(pre_ele(norbLR)); pre_ele=0
          iorbLR0=1
          iorbLR =0
          irrep_pre=0
          do i=1,8
            iorbLR=iorbLR+LRras2(i)
            if(i.eq.1)then
            else
              irrep_pre=irrep_pre+irrep_diff(i-1)
            end if
            do j=iorbLR0,iorbLR,1
              pre_ele(j)=irrep_pre
              write(6,*)"j,pre_ele(j)",j,pre_ele(j)
            end do
            iorbLR0=iorbLR+1 ! At least work for C1
          end do

          write(6,*)"pre_ele, ndets_RGLR", pre_ele, ndets_RGLR
          write(6,*)"nalpha,  nbeta     ", nele_alpha, nele_beta

          !> DETs read from mclr_dets.initial
          open(UNIT=117,file="mclr_dets.initial",status="OLD")
          allocate(SD_DMRG(ndets_RGLR))
          do idet=1,ndets_RGLR
            allocate(SD_DMRG(idet)%electron(neletol))
            allocate(SD_DMRG(idet)%ele_conf(norb))
            allocate(SD_DMRG(idet)%dv(nstates_RGLR))
            SD_DMRG(idet)%electron=0
            SD_DMRG(idet)%ele_conf=0
            SD_DMRG(idet)%dv=0.0d0
            read(117,"(1X,I8,6X)",advance='no')SD_DMRG(idet)%ITYPE
            do i=1,neletol
              read(117,"(1X,I5)",advance='no')SD_DMRG(idet)%electron(i)
            end do
            read(117,"(5X,I3)",advance='no')SD_DMRG(idet)%isign
            read(117,"(11X,I20)")SD_DMRG(idet)%inum
          end do
          close(117)

          write(6,*)"before get the executable file"

          !> get the executable file
         call systemf("cp /home/eth/yma/Maquis_MPSLR/build/applications/
     &srcas/srcas $PWD",rc)

          write(6,*)"before get the executable file"

          allocate(ele_orb_alpha(norb))
          allocate(ele_orb_beta(norb))
          !> All of the DETs into Maquis format
          open(unit=118,file="dets.mclr")
          do idet=1,ndets_RGLR
!            write(6,*)"idet",idet
            ele_orb_alpha=0
            ele_orb_beta=0
            do i=1,neletol
              if(SD_DMRG(idet)%electron(i)>0)then ! With preconditioner
                j=abs(SD_DMRG(idet)%electron(i))
                ele_orb_alpha(j+pre_ele(j))=1
              else
                j=abs(SD_DMRG(idet)%electron(i))
                ele_orb_beta(j+pre_ele(j))=1
              end if
            end do
            !> The same style also used in Maquis input
            do i=1,norb
              if(ele_orb_alpha(i).eq.1.and.ele_orb_beta(i).eq.1)
     &           SD_DMRG(idet)%ele_conf(i)=4
              if(ele_orb_alpha(i).eq.1.and.ele_orb_beta(i).eq.0)
     &           SD_DMRG(idet)%ele_conf(i)=3
              if(ele_orb_alpha(i).eq.0.and.ele_orb_beta(i).eq.1)
     &           SD_DMRG(idet)%ele_conf(i)=2
              if(ele_orb_alpha(i).eq.0.and.ele_orb_beta(i).eq.0)
     &           SD_DMRG(idet)%ele_conf(i)=1
              !write(6,*)"SD_DMRG(idet)%ele_conf(i)",SD_DMRG(idet)%ele_conf(i)
            end do
            do i=1,norb
              write(118,"(I1)",advance='no')SD_DMRG(idet)%ele_conf(i)
!              write(6,*)"SD_DMRG(idet)%ele_conf(",i,")",
!     &                   SD_DMRG(idet)%ele_conf(i)
            end do
            write(118,*)
          end do
          close(118)

          write(6,*)"After write dets.mclr file"

          !> Test for part of DETS
          !> The way of "open file" need to be written
          !                                   Yingjin 2015.8.13
          call systemf("wc -l dets.mclr > dets.mclr.info",rc)
          open(unit=118,file="dets.mclr.info")
            read(118,*)ndets_total
          close(118)
          !> If too many determinants,
          !> try to use the single, double, triple gradually untill 9999 (as the maximum)
          if(ndets_total.gt.9999)then
            call systemf("head -9999 dets.mclr > ELE_CISR_FOR_MCLR",rc)
          end if

          !> Recover the determinants, off-diagional multiply 2
          !> ========= should be improved by Hash etc. ==========
          do i=istate,istate
             !write(6,*)"SD_DMRG(idet)%dv(i)",i
            open(unit=118,file="GET_COEFF_IN_LIST")
              inquire(file='ELE_CISR_FOR_MCLR',exist=IFFILE)
              if(IFFILE)then
                tmp_run="./srcas "//trim(checkpoint(i))//
     &       " dets.mclr 1.0 1.0 0 ELE_CISR_FOR_MCLR > "//"CIRE.scratch"
              else
                tmp_run="./srcas "//trim(checkpoint(i))//
     &               " dets.mclr 1.0 1.0 0 dets.mclr > "//"CIRE.scratch"
              end if
              write(118,*)trim(tmp_run)
            close(118)
            call systemf("chmod +x GET_COEFF_IN_LIST",rc)
            call systemf("./GET_COEFF_IN_LIST",rc)
            !> read in the dets-coefficients
            open(unit=118,file="det_coeff.tmp")
            read(118,*)ndets_mclr
            do idet=1,ndets_mclr
              read(118,*)idx_det,SD_DMRG(idx_det)%dv(i)
!              write(6,*)idx_det,SD_DMRG(idx_det)%dv(i)
            end do
            close(118)
            !>  off-diagional multiply 2
            do idet=1,ndets_RGLR
              if(SD_DMRG(idet)%ITYPE.eq.1)then
                SD_DMRG(idet)%dv(i)=-1.0d0*SD_DMRG(idet)%dv(i)
              else
                dtmp=sqrt((SD_DMRG(idet)%dv(i)**2)*2.0d0)
                SD_DMRG(idet)%dv(i)=
     &         -1.0d0*dsign(dtmp,SD_DMRG(idet)%dv(i))
              end if
!              write(6,*)"i,idet,dv",i,idet,SD_DMRG(idet)%dv(i)
            end do
          end do

          dtmp=0.0d0
          indexSD=0
          vector=0.0d0
          do i=1,ndets_RGLR
            indexSD(i)=i !SD_DMRG(i)%inum*SD_DMRG(i)%isign
            vector(i)=SD_DMRG(i)%dv(istate)
            dtmp=dtmp+vector(i)**2
          end do

          write(6,*)"nele_alpha","nele_beta",nele_alpha,nele_beta
          write(6,*)"Total CI weight is ", dtmp
          write(6,*)" ================================================ "
          write(6,*)"  IF for frequency, the CI weight must be  "
          write(6,*)"      very close to 1 (i.e. 0.9999)"
          write(6,*)" ------------------------------------------------ "
          write(6,*)"  IF gradients in state-averaged case,  "
          write(6,*)"      even very few is stil OK (e.g. 0.1)"
          write(6,*)"      however, better around 0.9 "
          write(6,*)" ================================================ "
          call xflush(6)

!          stop

          deallocate(pre_ele)
          deallocate(SD_DMRG)

      end Subroutine ci_reconstruct

