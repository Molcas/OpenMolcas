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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

subroutine RSBB1E_MCLR(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,ISEL1,ISEL3,ICEL1,ICEL3,SB,CB,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,SSCR,CSCR,I1, &
                       XI1S,H,NSM,SGN)
! One electron excitations on column strings
!
! =====
! Input
! =====
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP        : String group of columns
! NROW        : Number of rows in S and C block
! ISEL1(3)    : Number of electrons in RAS1(3) for S block
! ICEL1(3)    : Number of electrons in RAS1(3) for C block
! CB          : Input C block
! NTSOB       : Number of orbitals per type and symmetry
! IBTSOB      : base for orbitals of given type and symmetry
! IBORB       : Orbitals of given type and symmetry
! NSM         : Number of symmetries of orbitals
! MAXI        : Largest Number of "spectator strings" treated simultaneously
! MAXK        : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! SB : updated sigma block
!
! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the largest
!              number of orbital pairs of given symmetries and types.
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! H          : Space for one electron integrals
!
! Jeppe Olsen, Winter of 1991

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISCSM, ISCTP, ICCSM, ICCTP, IGRP, NROW, ISEL1, ISEL3, ICEL1, ICEL3, NTSOB(3,*), IBTSOB(3,*), &
                                 ITSOB(*), MAXI, MAXK, NSM
real(kind=wp), intent(inout) :: SB(*)
real(kind=wp), intent(in) :: CB(*), SGN
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), XI1S(*), H(*)
integer(kind=iwp), intent(_OUT_) :: I1(*)
integer(kind=iwp) :: IBORB, IBOT, ICGOFF, IIORB, IJSM, IJTP, IORB, IPART, ISBOFF, ISM, ITOP, ITP(3), ITYP, JBORB, JJORB, JORB, &
                     JSM, JTP(3), JTYP, KBOT, KEND, KTOP, NIBTC, NIK, NIORB, NIPART, NJORB, NKBTC, NSXTP

! Type of single excitations that connects the two column strings

call SXTYP(NSXTP,ITP,JTP,ISEL1,ISEL3,ICEL1,ICEL3)

! Symmetry of single excitation that connects IBSM and JBSM

IJSM = Mul(ISCSM,ICCSM)
if (IJSM == 0) return

do IJTP=1,NSXTP
  ITYP = ITP(IJTP)
  JTYP = JTP(IJTP)
  do ISM=1,NSM
    JSM = Mul(ISM,IJSM)
    if (JSM == 0) cycle
    NIORB = NTSOB(ITYP,ISM)
    NJORB = NTSOB(JTYP,JSM)
    if ((NIORB == 0) .or. (NJORB == 0)) cycle
    ! Loop over partitionings of the row strings
    NIPART = NROW/MAXI
    if (NIPART*MAXI /= NROW) NIPART = NIPART+1
    do IPART=1,NIPART
      IBOT = (IPART-1)*MAXI+1
      ITOP = min(IBOT+MAXI-1,NROW)
      NIBTC = ITOP-IBOT+1
      ! Loop over partitionings of N-1 strings
      KBOT = 1-MAXK
      KTOP = 0
      do
        KBOT = KBOT+MAXK
        KTOP = KTOP+MAXK
        JBORB = IBTSOB(JTYP,JSM)
        NJORB = NTSOB(JTYP,JSM)
        ! Obtain CSCR(JORB,I,K) = SUM(J)<K!A JORB!J>C(I,J)
        ! Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
        do JJORB=1,NJORB
          !if (JJORB == 1) then
          !  JOFF = 1
          !else
          !  JOFF = IOFF+NK
          !  ! Note : NK is not not known until first call to ADST
          !end if
          JORB = ITSOB(JBORB-1+JJORB)
          ! set up I1(K) =  XI1S(K) a JORB !J STRING >
          call ADST(JORB,1,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND)
          ! Gather  C Block
          ! First index : JORB, second index : JaKb
          ICGOFF = 1+(JJORB-1)*NKBTC*NIBTC
          call MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,NKBTC,I1,XI1S)
        end do
        ! Obtain one electron integrals (JSM,JTP,ISM,ITP)
        call NGETH1(SSCR,ISM,ITYP,JSM,JTYP)
        call TRNSPS(NIORB,NJORB,SSCR,H)

        !Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
        NIK = NIBTC*NKBTC

        call DGEMM_('N','N',NIK,NIORB,NJORB,One,CSCR,max(NIK,1),H,max(1,NJORB),Zero,SSCR,max(1,NIK))
        !S(I,a+K) = S(I,a+K)+sgn*Sscr(I,K,i)
        IBORB = IBTSOB(ITYP,ISM)
        do IIORB=1,NIORB
          IORB = ITSOB(IBORB-1+IIORB)
          ! set up I1(IORB,K) = a IORB !I STRING >
          call ADST(IORB,1,ISCTP,ISCSM,IGRP,KBOT,KTOP,I1,XI1S,1,NKBTC,KEND)
          ! Well, someplace the minus must come in
          if (SGN == -One) XI1S(1:NKBTC) = -XI1S(1:NKBTC)
          ISBOFF = 1+(IIORB-1)*NKBTC*NIBTC
          call MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,NKBTC,I1,XI1S)
        end do
        ! end of this K partitioning
        if (KEND /= 0) exit
      end do
    end do
    ! End of loop over I partitioninigs
  end do
  ! (end of loop over symmetries)
end do

end subroutine RSBB1E_MCLR
