************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016,2017, Giovanni Li Manni                           *
*               2019, Oskar Weser                                      *
************************************************************************

      module fciqmc_read_RDM
      use stdalloc, only : mma_allocate, mma_deallocate
      use general_data, only : iSpin, nActEl
      use rasscf_data, only : nAc, nAcPar
! Note that two_el_idx_flatten has also out parameters.
      use index_symmetry, only : two_el_idx_flatten
      private
      public :: read_neci_RDM, cleanup, read_neci_GUGA_RDM
      contains

!>  @brief
!>    Start and control FCIQMC.
!>
!>  @author Werner Dobrautz
!>
!>  @details
!>  Read the spin-free TwoRDM file written by the GUGA-NECI
!>  implementation and transfer them to Molcas.
!>  for now it also outputs 'fake' additional spin-dependent files
!>  because I don't know yet how I have to handle them in the further
!>  Molcas control flow

!>  @paramin[out] DMAT Average spin-free 1 body density matrix
!>  @paramin[out] DSPN 'fake' spin-dependent 1-RDM
!>  @paramin[out] PSMAT Average spin-free 2 body density matrix
!>  @paramin[out] PAMAT 'fake' Average antisymm. 2-dens matrix
      subroutine read_neci_GUGA_RDM(DMAT, DSPN, PSMAT, PAMAT)
      implicit none
#include "para_info.fh"
#include "output_ras.fh"
      real*8, intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
      integer :: iUnit, isfreeunit, i, iread, iprlev
      logical :: tExist
      real*8 :: RDMval
      parameter(routine = 'read_neci_GUGA_RDM')

      Call qEnter(routine)

      iprlev = iprloc(1)
      if(iprlev == debug) then
          write(6,*) 'Rank of process: ', MyRank
      end if

!     NOTE(Giovanni, Oskar): The suffix ".1" corresponds to the root.
!     For state averaged the .1 shall be replaced by iroot. (irdm in NECI)
      if(myRank /= 0) then
          call bcast_2RDM("PSMAT")
          call bcast_2RDM("PAMAT")
          call bcast_2RDM("DMAT")
      end if
**********************************************************************************
******************************** existency check *********************************
**********************************************************************************
      call f_Inquire('PSMAT',tExist)
      if(.not.tExist) goto 123
      call f_Inquire('PAMAT',tExist)
      if(.not.tExist) goto 123
      call f_Inquire('DMAT',tExist)
      if(.not.tExist) goto 123

      PSMAT(:) = 0.0d0
      PAMAT(:) = 0.0d0
      DMAT(:) = 0.0d0

**************************************************************
******************* read PSMAT *******************************
**************************************************************
      iUnit = IsFreeUnit(11)
      call Molcas_Open(iUnit, 'PSMAT')
      Rewind(iUnit)

      IF(IPRLEV >= DEBUG) THEN
          write(6,*) '    i     PSMAT   '
        write(6,*) ' ********************************************** '
      end if
      do
        read(iUnit, "(I6,G25.17)", iostat=iread) i, RDMval
        if(iread /= 0) exit
        psmat(i) = RDMval
        IF(IPRLEV >= DEBUG) THEN
           write(6,'(I6,G25.17)') i, RDMval
        end if
      end do
      close(iunit)

**************************************************************
******************* read PAMAT *******************************
**************************************************************
      iUnit = IsFreeUnit(11)
      call Molcas_Open(iUnit, 'PAMAT')
      Rewind(iUnit)

      IF(IPRLEV >= DEBUG) THEN
          write(6,*) '    i     PAMAT   '
        write(6,*) ' ********************************************** '
      end if
      do
        read(iUnit, "(I6,G25.17)", iostat=iread) i, RDMval
        if(iread /= 0) exit
        pamat(i) = RDMval
        IF(IPRLEV >= DEBUG) THEN
           write(6,'(I6,G25.17)') i, RDMval
        end if
      end do
      close(iunit)

**************************************************************
******************* read DMAT ********************************
**************************************************************
      iUnit = IsFreeUnit(11)
      call Molcas_Open(iUnit, 'DMAT')
      Rewind(iUnit)

      IF(IPRLEV >= DEBUG) THEN
          write(6,*) '    i     DMAT   '
        write(6,*) ' ********************************************** '
      end if
      do
        read(iUnit, "(I6,G25.17)", iostat=iread) i, RDMval
        if(iread /= 0) exit
        dmat(i) = RDMval
        IF(IPRLEV >= DEBUG) THEN
           write(6,'(I6,G25.17)') i, RDMval
        end if
      end do
      close(iunit)

      IF(IPRLEV >= DEBUG) THEN
        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',DMAT,NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',psmat,NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',pamat , NACPAR)
      end if

******* Clean evil non-positive semi-definite matrices. DMAT is input and output.
      ! hopefully this just works for GUGA-RDMs
      call cleanMat(DMAT)

      Call qExit(routine)
      Return

123   continue
          write(6,*) 'RDM file not found!'
          write(6,*) 'Probably file not generated by NECI?'
          Call QTrace()
          Call Abend()
      end subroutine read_neci_GUGA_RDM

!>  @brief
!>    Start and control FCIQMC.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @details
!>  Read TwoRDM files written by NECI and transfer them to Molcas.
!>  Neci can have some intermediate spin-resolved/spin-free RDMs where basically aaaa contains
!>  average of aaaa and bbbb, abab contains average of abab and baba...
!>  This is ok for CASSCF but not ok for spin-resolved properties, in which case the completely
!>  spin-resolved RDMs need to be read-in.
!>  In principle, NECI could also evaluate and store completely spin-free matrices.
!>  In that case only a reordering following Molcas convention is necessary.
!>
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] DSPN Average spin 1-dens matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
      subroutine read_neci_RDM(DMAT, DSPN, PSMAT, PAMAT)
      implicit none
#include "para_info.fh"
#include "output_ras.fh"
      real*8, intent(out) :: DMAT(:), DSPN(:), PSMAT(:), PAMAT(:)
      integer :: iUnit, isfreeunit, p, q, r, s, pq, rs, ps, rq, psrq,
     &  pqrs, iread, Nalpha, norb, iprlev
      logical :: tExist, switch
      real*8 :: fac, RDMval, fcalpha, fcbeta, fcnacte
      real*8 :: D_alpha(size(DMAT)), D_beta(size(DMAT))
      parameter(routine = 'read_neci_RDM')

      Call qEnter(routine)

      iprlev = iprloc(1)
      if(iprlev == debug) then
        write(6,*) 'Rank of process: ', MyRank
      end if
      switch = .true.
*     ^  This variable must become a keyword for discriminating spin-resolved from spin-free input RDMs.
*     ^  For now when .false. it is assumed that 3 files (only   aaaa, abab and abba) are fed.
*     ^  ....... when .true.  it is assumed that 6 files (adding bbbb, baba and baab) are fed.
*********************************************************************************
* Broadcasting TwoRDM generated by QMC code in master node into all processors. *
*********************************************************************************
! NOTE(Giovanni, Oskar): The suffix ".1" corresponds to the root.
!   For state averaged the .1 shall be replaced by iroot. (irdm in NECI)
      if(myRank /= 0) then
        call bcast_2RDM("TwoRDM_aaaa.1")
        call bcast_2RDM("TwoRDM_aaaa.1")
        call bcast_2RDM("TwoRDM_abab.1")
        call bcast_2RDM("TwoRDM_abba.1")
        call bcast_2RDM("TwoRDM_bbbb.1")
        call bcast_2RDM("TwoRDM_baba.1")
        call bcast_2RDM("TwoRDM_baab.1")
      end if
**********************************************************************************
******************************** existency check *********************************
**********************************************************************************
      call f_Inquire('TwoRDM_aaaa.1',tExist)
      if(.not.tExist) goto 123
      call f_Inquire('TwoRDM_abab.1',tExist)
      if(.not.tExist) goto 123
      call f_Inquire('TwoRDM_abba.1',tExist)
      if(.not.tExist) goto 123
      if(switch) then
        call f_Inquire('TwoRDM_bbbb.1',tExist)
        if(.not.tExist) goto 123
        call f_Inquire('TwoRDM_baba.1',tExist)
        if(.not.tExist) goto 123
        call f_Inquire('TwoRDM_baab.1',tExist)
        if(.not.tExist) goto 123
      end if

      D_alpha(:) = 0.0d0
      D_beta(:) = 0.0d0
      PSMAT(:) = 0.0d0
      PAMAT(:) = 0.0d0

      Nalpha = (nactel + iSpin - 1) / 2
      fac = merge(0.5d0, 1.0d0, switch)
      fcalpha = 1.0d0 / dble(nalpha - 1)
      fcbeta = 1.0d0 / dble(nactel - nalpha - 1)
      fcnacte = 1.0d0 / dble(nactel - 1)

*******************************************************************************************
*************************** Processing TwoRDM-AAAA ****************************************
*******************************************************************************************
      iUnit = IsFreeUnit(11)
      call Molcas_Open(iUnit, 'TwoRDM_aaaa.1')
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
        write(6,*) '    p     q     r     s    pq    rs   pqrs        ',
     &  'RDMval                PSMAT                   PAMAT'
        write(6,*) ' ********************** AAAA ****************** '
      end if
      do
******************* processing as PQRS ***********************
**************************************************************
        read(iUnit, "(4I6,G25.17)", iostat=iread) s, q, r, p, RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
        if (r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
        if (p > q.and.r > s) PAMAT(pqrs) = PAMAT(pqrs) + fac * RDMval
        if (p > q.and.r < s) PAMAT(pqrs) = PAMAT(pqrs) - fac * RDMval
        IF(IPRLEV >= DEBUG) THEN
           write(6,'(7I6,3G25.17)')
     &     p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        END IF
******* Contribution to D_alpha (not final):
        if (p == q) D_alpha(rs) = D_alpha(rs) + RDMval
        if (r == s) D_alpha(pq) = D_alpha(pq) + RDMval
******************* processing as PSRQ ***********************
**************************************************************
        psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
        ! WD: why twice?? on purpose or bug?
        psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
******* Contribution to PSMAT and PAMAT:
        if (r <= q) then
          PSMAT(psrq) = PSMAT(psrq) - fac * RDMval
          if (r /= q) PAMAT(psrq) = PAMAT(psrq) + fac * RDMval
        end if
        if (r > q) then
          PSMAT(psrq) = PSMAT(psrq) - fac * RDMval
          PAMAT(psrq) = PAMAT(psrq) - fac * RDMval
        end if
        IF(IPRLEV >= DEBUG) THEN
          write(6,'(7I6,3G25.17)')
     &      p,s,r,q,ps,rq,psrq, RDMval,PSMAT(psrq),PAMAT(psrq)
        END IF
******* Contribution to D_alpha (not final):
* The minus sign comes from the fact that in NECI these elements have opposite sign
* compared to the element in normal order, that is d_pqrs = -d_psrq.
        if (p == s) D_alpha(rq) = D_alpha(rq) - RDMval
        if (r == q) D_alpha(ps) = D_alpha(ps) - RDMval
      end do
      close(iunit)

*******************************************************************************************
*************************** Processing TwoRDM-BBBB ****************************************
*******************************************************************************************
      if (switch) then
        iUnit=IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_bbbb.1')
        Rewind(iUnit)
        if (IPRLEV >= DEBUG) then
         write(6,*) '    p     q     r     s    pq    rs   pqrs       ',
     &   'RDMval                PSMAT                   PAMAT'
         write(6,*) ' ********************** BBBB ****************** '
        end if
        do
******************* processing as PQRS ***********************
**************************************************************
          read(iUnit,"(4I6,G25.17)",iostat=iread) s, q, r, p, RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
          if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
          if(p > q.and.r > s) PAMAT(pqrs) = PAMAT(pqrs) + fac * RDMval
          if(p > q.and.r < s) PAMAT(pqrs) = PAMAT(pqrs) - fac * RDMval
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &         p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
******* Contribution to D_beta (not final):
          if (p == q) D_beta(rs) = D_beta(rs) + RDMval
          if (r == s) D_beta(pq) = D_beta(pq) + RDMval
******************* processing as PSRQ ***********************
**************************************************************
          psrq = two_el_idx_flatten(p, s, r, q, ps, rq)
******* Contribution to PSMAT and PAMAT:
          if(r <= q) then
            PSMAT(psrq) = PSMAT(psrq) - fac*RDMval
            if(r /= q) PAMAT(psrq) = PAMAT(psrq) + fac*RDMval
          end if
          if(r > q) then
            PSMAT(psrq) = PSMAT(psrq) - fac*RDMval
            PAMAT(psrq) = PAMAT(psrq) - fac*RDMval
          end if
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &        p,s,r,q,ps,rq,psrq, RDMval,PSMAT(psrq),PAMAT(psrq)
          END IF
******* Contribution to D_beta (not final):
* The minus sign comes from the fact that in NECI these elements have opposite sign
* compared to the element in normal order, that is d_pqrs = -d_psrq.
          if(p == s) D_beta(rq)=D_beta(rq)-RDMval
          if(r == q) D_beta(ps)=D_beta(ps)-RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
*******************************************************************************************
*************************** Processing TwoRDM-ABAB ****************************************
*******************************************************************************************
      iUnit=IsFreeUnit(11)
      Call Molcas_Open(iUnit,'TwoRDM_abab.1')
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
       write(6,*) ' ********************** ABAB ****************** '
      END IF
      do
        read(iUnit,"(4I6,G25.17)",iostat=iread) s,q,r,p,RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
        if(r > s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
        if(r < s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
        if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
        if (IPRLEV >= DEBUG) then
          write(6,'(7I6,3G25.17)')
     &        p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        end if
******* Contribution to D_alpha and D_beta (not final):
        if (p == q) D_alpha(rs) = D_alpha(rs) + RDMval
        if (r == s .and. p /= r) D_beta(pq) = D_beta(pq) + RDMval
      end do

******* Copy D_beta to D_alpha and clean D_beta again for further use:
      if (.not. switch) then
        D_alpha(:) = D_beta(:) + D_alpha(:)
        D_beta(:) = 0.0d0
      end if
      close(iunit)
*******************************************************************************************
*************************** Processing TwoRDM-BABA ****************************************
*******************************************************************************************
      if (switch) then
        iUnit = IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_baba.1')
        rewind(iUnit)
        if (IPRLEV >= DEBUG) then
           write(6,*) ' ********************** BABA ****************** '
        end if
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) s,q,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) + fac * RDMval
          if(r > s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
          if(r < s.and.p /= q) PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
          if(r /= s.and.p == q) PSMAT(pqrs) = PSMAT(pqrs) + fac*RDMval
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &          p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
******* Contribution to D_alpha (not final):
          if(p == q) D_beta(rs) = D_beta(rs) + RDMval
          if(r == s.and.p /= r) D_alpha(pq) = D_alpha(pq)+RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
*******************************************************************************************
*************************** Processing TwoRDM-ABBA ****************************************
*******************************************************************************************
      iUnit=IsFreeUnit(11)
      Call Molcas_Open(iUnit,'TwoRDM_abba.1')
      Rewind(iUnit)
      IF(IPRLEV >= DEBUG) THEN
        write(6,*) ' ********************** ABBA ****************** '
      END IF
      do
        read(iUnit,"(4I6,G25.17)",iostat=iread) q,s,r,p,RDMval
        if(iread /= 0) exit
        pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
        PSMAT(pqrs) = PSMAT(pqrs) - fac*RDMval
        if(r < s) then
          PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
        end if
        if(r > s) then
          PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
        end if
        IF(IPRLEV >= DEBUG) THEN
          write(6,'(7I6,3G25.17)')
     &        p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
        END IF
******* Contribution to D_alpha (not final):
          if(r == s) D_alpha(pq)=D_alpha(pq)-RDMval
          if(p == q) D_beta(rs)=D_beta(rs)-RDMval
      end do
      if (.not. switch) then
        D_alpha(:) = D_beta(:) + D_alpha(:)
        D_beta(:) = 0.0d0
      end if
      close(iunit)
*******************************************************************************************
*************************** Processing TwoRDM-BAAB ****************************************
*******************************************************************************************
      if(switch) then
        iUnit=IsFreeUnit(11)
        Call Molcas_Open(iUnit,'TwoRDM_baab.1')
        Rewind(iUnit)
        IF(IPRLEV >= DEBUG) THEN
          write(6,*) ' ********************** BAAB ****************** '
        END IF
        do
          read(iUnit,"(4I6,G25.17)",iostat=iread) q,s,r,p,RDMval
          if(iread /= 0) exit
          pqrs = two_el_idx_flatten(p, q, r, s, pq, rs)
******* Contribution to PSMAT and PAMAT:
          PSMAT(pqrs) = PSMAT(pqrs) - fac*RDMval
          if(r < s) then
            PAMAT(pqrs) = PAMAT(pqrs) + fac*RDMval
          end if
          if(r > s) then
            PAMAT(pqrs) = PAMAT(pqrs) - fac*RDMval
          end if
          IF(IPRLEV >= DEBUG) THEN
            write(6,'(7I6,3G25.17)')
     &          p,q,r,s,pq,rs,pqrs, RDMval,PSMAT(pqrs),PAMAT(pqrs)
          END IF
******* Contribution to D_alpha (not final):
          if(p == q) D_alpha(rs)=D_alpha(rs)-RDMval
          if(r == s) D_beta(pq)=D_beta(pq)-RDMval
        end do
        close(iunit)
      end if ! End statement for spin-resolved RDMs.
*******************************************************************************************
***************************   Final Updates to RDMs  **************************************
*******************************************************************************************
      if (.not.switch) D_beta(:) = D_alpha(:)
      D_alpha(:) = fcnacte * D_alpha(:)
      D_beta(:) = fcnacte * D_beta(:)
      DSPN(:) = D_Beta(:) - D_alpha(:)
      DMAT(:) = D_Beta(:) + D_alpha(:)

******* Clean evil non-positive semi-definite matrices. DMAT is input and output.
      call cleanMat(DMAT)

      IF(IPRLEV >= DEBUG) THEN
       norb  = (int(sqrt(dble(1 + 8 * size(DMAT)))) - 1) / 2
       call triprt('D_alpha in neci2molcas',' ',D_alpha,norb)
       call triprt('D_beta  in neci2molcas',' ',D_beta ,norb)
       call triprt('DMAT in neci2molcas',' ',DMAT,norb)
       call triprt('DSPN in neci2molcas',' ',DSPN,norb)
      END IF
      Call qExit(routine)
      Return

123   continue
          write(6,*) 'RDM files not found!'
          write(6,*) 'Probably file not generated by NECI?'
          Call QTrace()
          Call Abend()
      end subroutine read_neci_RDM

      subroutine bcast_2RDM(InFile)
        use filesystem, only : symlink_, strerror_, get_errno_
        implicit none
        character(*), intent(in) :: InFile
        character(1024) :: master
        integer :: lmaster1, err

        call prgmtranslate_master(InFile, master, lmaster1)
        call symlink_(trim(master), trim(InFile), err)
        if (err == 0) write(6, *) strerror_(get_errno_())
      end subroutine bcast_2RDM

      Subroutine CleanMat(MAT)
************* by G. Li Manni Stuttgart April 2016 *************
*
* MAT: One-body density matrix in MO basis as passed by QMC calculation.

* It could well be an average matrix in SA calculation.
*
* It has following shape:
*        11
*        12 22
*        ** ** 33
*        ** ** ** 44
*        ** ** ** 45 55
*        ** ** ** 46 56 66
*        ** ** ** 47 57 67 77
*        ** ** ** ** ** ** ** 88
*        ** ** ** ** ** ** ** 89  99
*        ** ** ** ** ** ** ** 810 910 1010
*        """""""""""""""""""""""""""""""""""
* mimicking a system with (2 0 0 1 4 3 0 0)  actice orbitals (blocked by Irreps)

*           DMAT will be destroyed and replaced with a positive semi-definite one.
*           N-representability will be preserved.

      implicit none
#include "WrkSpc.fh"
* NACPAR = NAC*(NAC+1)/2 with NAC total number of active orbitals
      real*8, intent(inout) :: MAT(NacPar)
      real*8, allocatable :: EVC(:), Tmp(:), Tmp2(:), MAT_copy(:)
      integer :: rc, i, j
      real*8 :: trace
      character(12), parameter :: routine = 'CleanMat'
      logical :: cleanup_required

      Call qEnter(routine)

      rc = 0
      If (nacpar .lt. 1) then
        rc= -1
        write(6,*) 'matrix size < 1.'
        Go To 10
      end if

      call mma_allocate(MAT_copy, NacPar)
      MAT_copy(:) = MAT(:)

* Allocate memory for eigenvectors and new DMAT
      call mma_allocate(EVC, NAC**2)
* Initialize eigenvectors
      Call dCopy_(NAC**2, [0.0d0], 0, EVC, 1)
* set eigenvector array to identity for this version of JACOB
      Call dCopy_(NAC, [1.0d0], 0, EVC, NAC + 1)

* Step 1: Diagonalize MAT. Eigenvalues are stored in diagonal of MAT
      trace = 0.0d0
      do i = 1, nac
         trace = trace + mat(i * (i + 1) / 2)
      end do
      CALL JACOB(MAT_copy, EVC, NAC, NAC)

#ifdef _DEBUG_
      write(6,*) 'eigenvalues: '
      do i=1,nac
         write(6,*) MAT_copy(I*(I+1)/2)
      end do
      write(6,*) 'eigenvectors: '
      do i=1, nac
        write(6,*) (EVC(i * NAC + j), j = 0, NAC)
      end do
#endif
* Set to zero negative eigenvalue and to TWO values larger than 2.0d0.
      cleanup_required = .false.
      do j = 1, nac
        if (MAT_copy(j * (j + 1) / 2) > 2.0d0) then
          MAT_copy(j * (j + 1) / 2) = 2.0d0
          cleanup_required = .true.
        end if
        if (MAT_copy(j * (j + 1) / 2) < 1.0d-12) then
          MAT_copy(j * (j + 1) / 2) = 0.0d0
          cleanup_required = .true.
        end if
      end do

      if (cleanup_required) then
        trace = 0.0d0
        do i = 1, nac
          trace = trace + MAT_copy(I * (I + 1) / 2)
        end do
        write(6,*) 'trace after removing negative eigenvalues =', trace
* Combine pieced to form the output MAT
* blas routine for square*triangular operation
        call mma_allocate(Tmp, nac**2)
        call mma_allocate(Tmp2, nac**2)
        Call dCopy_(nac**2, [0.0d0], 0, Tmp, 1)
        Call dCopy_(nac**2, [0.0d0], 0, Tmp2, 1)
        do i = 1, nac
          do j = 1, nac
            Tmp(j + (i - 1) * nac) =
     &          EVC(j + (i - 1) * NAC) * MAT_copy(I * (I + 1) / 2)
          end do
        end do
        Call DGEMM_('N','T',nac,nac,nac,
     &              1.0d0, Tmp, nac, EVC, nac,
     &              0.0d0, Tmp2, nac)
* Copy back to MAT
        do i = 1, nac
          do j = 1, i
            MAT(j + (i - 1) * i / 2) = Tmp2(j + (i - 1) * nac)
          end do
        end do
#ifdef _DEBUG_
        write(6,*) 'trace after recombination:'
        trace = 0.0d0
        do i = 1, nac
           trace = trace + MAT(i * (i + 1) / 2)
        end do
#endif
        call mma_deallocate(tmp)
        call mma_deallocate(tmp2)
      end if
      call mma_deallocate(MAT_copy)
      call mma_deallocate(EVC)
****************** Exit ****************
10    Continue
      Call qExit(routine)
      return
      end subroutine cleanMat

      subroutine cleanup()
        implicit none
        ! Add your deallocations here.
        ! This routine will be called when exiting rasscf.
        continue
      end subroutine
      end module fciqmc_read_RDM
