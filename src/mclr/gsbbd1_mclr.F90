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
! Copyright (C) 1991,1995, Jeppe Olsen                                 *
!***********************************************************************

subroutine GSBBD1_MCLR(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,ISEL,ICEL,SB,CB,NOBPTS,IOBPTS,MAXI,MAXK, &
                       SSCR,CSCR,I1,XI1S,I2,XI2S,NSMOB,RHO1S)
! Contributions to one electron density matrix from column excitations
!
! GAS version, August 95, Jeppe Olsen
!
! =====
! Input
! =====
! RHO1        : One body density matrix to be updated
! NACOB       : Number of active orbitals
! ISCSM,ISCTP : Symmetry and type of sigma columns
! ICCSM,ICCTP : Symmetry and type of C     columns
! IGRP        : String group of columns
! NROW        : Number of rows in S and C block
! NGAS        : Number of active spaces
! ISEL        : Number of electrons per AS for S block
! ICEL        : Number of electrons per AS for C block
! CB          : Input C block
! NOBPTS      : Number of orbitals per type and symmetry
! IOBPTS      : base for orbitals of given type and symmetry
! IBORB       : Orbitals of given type and symmetry
! NSMOB,NSMDX : Number of symmetries of orbitals, double excitations
! MAXI        : Largest Number of "spectator strings" treated simultaneously
! MAXK        : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! RHO1 : Updated density block
!
! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the largest
!              number of orbital pairs of given symmetries and types.
! I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
! I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
!              type and symmetry
! RHO1S      : Space for one electron density
!
! Jeppe Olsen, Winter of 1991
! Updated for GAS, August '95

use Symmetry_Info, only: Mul
use Constants, only: Zero, One

implicit real*8(A-H,O-Z)
! General input
integer NOBPTS(3,*), IOBPTS(3,*)
!integer NTSOB(3,*),IBTSOB(3,*)
! Input
integer ISEL(NGAS), ICEL(NGAS)
dimension CB(*), SB(*)
! Output
dimension RHO1(*)
! Scatch
dimension SSCR(*), CSCR(*), RHO1S(*)
dimension I1(*), XI1S(*)
dimension I2(*), XI2S(*)
! Local arrays
dimension ITP(3*3), JTP(3*3)

! Type of single excitations that connects the two column strings
call SXTYP_GAS(NSXTP,ITP,JTP,3,ISEL,ICEL)
! Symmetry of single excitation that connects IBSM and JBSM
IJSM = Mul(ISCSM,ICCSM)
if (IJSM == 0) return
do IJTP=1,NSXTP
  ITYP = ITP(IJTP)
  JTYP = JTP(IJTP)
  do ISM=1,NSMOB
    ! new i and j so new intermediate strings
    JSM = Mul(ISM,IJSM)
    if (JSM == 0) cycle
    NIORB = NOBPTS(ITYP,ISM)
    NJORB = NOBPTS(JTYP,JSM)
    IBIORB = IOBPTS(ITYP,ISM)
    IBJORB = IOBPTS(JTYP,JSM)
    if ((NIORB == 0) .or. (NJORB == 0)) cycle

    !OLD Loop over partitionings of the row strings
    NIPART = NROW/MAXI
    if (NIPART*MAXI /= NROW) NIPART = NIPART+1
    ! Loop over partitionings of N-1 strings
    KBOT = 1-MAXK
    KTOP = 0
    do
      KBOT = KBOT+MAXK
      KTOP = KTOP+MAXK
      !EAWBEGIN970207
      !     -1 -> MAXK
      !     KTOP -> -1
      !     KTOP=-1
      ! Single excitation information independent of I strings

      ! set up I1(K) =  XI1S(K) a JORB !J STRING >
      call ADST(IBJORB,NJORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND)
      ! set up I2(K) =  XI1S(K) a JORB !J STRING >
      call ADST(IBIORB,NIORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,I2,XI2S,MAXK,NKBTC,KEND)
      !EAWEND
      ! Appropriate place to start partitioning over I strings
      ! Loop over partitionings of the row strings
      do IPART=1,NIPART
        IBOT = (IPART-1)*MAXI+1
        ITOP = min(IBOT+MAXI-1,NROW)
        NIBTC = ITOP-IBOT+1

        ! Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
        ! Gather  C Block
        do JJORB=1,NJORB
          ICGOFF = 1+(JJORB-1)*NKBTC*NIBTC
          call MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,NKBTC,I1(1+(JJORB-1)*MAXK),XI1S(1+(JJORB-1)*MAXK))
        end do

        ! Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
        do IIORB=1,NIORB
          ! Gather S Block
          ISGOFF = 1+(IIORB-1)*NKBTC*NIBTC
          call MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,NKBTC,I2(1+(IIORB-1)*MAXK),XI2S(1+(IIORB-1)*MAXK))
        end do
        NKI = NKBTC*NIBTC
        if (NKI*NIORB*NJORB /= 0) then
          call DGEMM_('T','N',NIORB,NJORB,NKI,One,SSCR,NKI,CSCR,NKI,Zero,RHO1S,NIORB)
        else
          RHO1S(1:NIORB*NJORB) = Zero
        end if
        ! Scatter out to complete matrix
        do JJORB=1,NJORB
          JORB = IBJORB-1+JJORB
          do IIORB=1,NIORB
            IORB = IBIORB-1+IIORB
            RHO1((JORB-1)*NACOB+IORB) = RHO1((JORB-1)*NACOB+IORB)+RHO1S((JJORB-1)*NIORB+IIORB)
          end do
        end do

      end do
      ! /\ end of this I partitioning
      ! end of this K partitioning
      if (KEND /= 0) exit
    end do
    ! End of loop over I partitioninigs
  end do
  ! (end of loop over symmetries)
end do

return

end subroutine GSBBD1_MCLR
