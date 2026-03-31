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

subroutine NEWSSTAB(ORBTAB)

use rassi_global_arrays, only: SSTAB
use Cntrl, only: MORSBITS
use Symmetry_Info, only: MUL, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ORBTAB(*)
integer(kind=iwp) :: I, IERR, IKMS2, IKSCR, IKSYM, IMS2, IPOP, ISBS, ISBS1, ISBS2, ISBS3, ISBS4, ISBS5, ISBS6, ISBS7, ISBS8, ISCR, &
                     ISGN, ISORB, ISPART, ISSTP, ISYM, ITYPE, J, KMRSSBS, KMS2, KOINFO, KSBSANN, KSBSCRE, KSBSMRS, KSSTANN, &
                     KSSTCRE, KSSTP, KSYM, LPOS, MRS, MS2, MXMS2, MXPOP, N, NASPO, NASPRT, NEWMRS, NEWSBS, NEWSST, NO, NSBS, &
                     NSBSTOT, NSCR, NSSTP, NTAB
integer, external :: MORSANN, MORSCRE, MORSPOP, MORSSPIN, MORSSYMM
integer, allocatable :: NOSUB(:), OSPN(:), OSYM(:), SCR(:), SCR2(:)

! Table type ID:
ITYPE = 19
! Pick up some data from the orbital table:
NASPO = ORBTAB(4)
nIrrep = ORBTAB(5)
NASPRT = ORBTAB(9)
KOINFO = 19
! Make temporary arrays for spin label and symmetry of each orbital
call mma_allocate(OSPN,NASPO,Label='OSPN')
call mma_allocate(OSYM,NASPO,Label='OSYM')
! Make a temporary array, which gives the number of spin orbitals in
! each subpartition:
call mma_allocate(NOSUB,NASPRT,Label='NOSUB')
NOSUB(:) = 0
do ISORB=1,NASPO
  ISYM = OrbTab(KOINFO+1+(ISORB-1)*8)
  OSYM(ISORB) = ISYM
  MS2 = OrbTab(KOINFO+3+(ISORB-1)*8)
  OSPN(ISORB) = MS2
  ISPART = OrbTab(KOINFO+6+(ISORB-1)*8)
  N = 1+NOSUB(ISPART)
  NOSUB(ISPART) = N
end do
! We need a temporary array, NSBSSCR(nIrrep,0:NPOP,-MXMS2:MXMS2,NASPRT),
! to keep the number of substrings of different kind:
MXPOP = MORSBITS
MXMS2 = MORSBITS
NSCR = nIrrep*(1+MXPOP)*(2*MXMS2+1)*NASPRT
call mma_allocate(SCR,NSCR,Label='SCR')
! Addressing will be through the cumbersome formula
! ISCR=ISYM+nIrrep*(IPOP+(1+MXPOP)*(MXMS2+MS2+(2*MXMS2+1)*(ISPART-1)))
! Initialize counter of substrings:
do ISPART=1,NASPRT
  do IMS2=-MXMS2,MXMS2
    do IPOP=0,MXPOP
      do ISYM=1,nIrrep
        ! NSBSSCR(ISYM,IPOP,IMS2,ISPART)=0:
        ISCR = ISYM+nIrrep*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*(ISPART-1)))
        SCR(ISCR) = 0
      end do
    end do
  end do
  ! NSBSSCR(ISYM=1,IPOP=0,IMS2=0,ISPART)=1:
  ISCR = 1+nIrrep*(0+(1+MXPOP)*(MXMS2+0+(2*MXMS2+1)*(ISPART-1)))
  SCR(ISCR) = 1
end do
! Compute number of substrings:
ISORB = 0
do ISPART=1,NASPRT
  NO = NOSUB(ISPART)
  do I=1,NO
    ISORB = ISORB+1
    KSYM = OSYM(ISORB)
    KMS2 = OSPN(ISORB)
    do IPOP=I,1,-1
      do IMS2=-IPOP,IPOP
        IKMS2 = IMS2-KMS2
        if (abs(IKMS2) <= IPOP-1) then
          do ISYM=1,nIrrep
            IKSYM = MUL(ISYM,KSYM)
            ISCR = ISYM+nIrrep*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*(ISPART-1)))
            N = SCR(ISCR)
            IKSCR = IKSYM+nIrrep*(IPOP-1+(1+MXPOP)*(MXMS2+IKMS2+(2*MXMS2+1)*(ISPART-1)))
            N = N+SCR(IKSCR)
            SCR(ISCR) = N
          end do
        end if
      end do
    end do
  end do
end do
NSSTP = 0
NSBSTOT = 0
do ISPART=1,NASPRT
  NO = NOSUB(ISPART)
  do IPOP=0,NO
    do ISYM=1,nIrrep
      do IMS2=-NO,NO
        ISCR = ISYM+nIrrep*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*(ISPART-1)))
        N = SCR(ISCR)
        if (N > 0) then
          NSSTP = NSSTP+1
          NSBSTOT = NSBSTOT+N
        end if
      end do
    end do
  end do
end do
! We finally know the number of substring types and the number of
! substrings. Transfer non-zero entries to the final table.

! Size of table, and offsets:
KSSTP = 15
KSSTANN = KSSTP+5*NSSTP
KSSTCRE = KSSTANN+MORSBITS*NSSTP
KSBSMRS = KSSTCRE+MORSBITS*NSSTP
KMRSSBS = KSBSMRS+2*NSBSTOT
KSBSANN = KMRSSBS+2*(2**MORSBITS)*NASPRT
KSBSCRE = KSBSANN+NSBSTOT*MORSBITS
NTAB = KSBSCRE+NSBSTOT*MORSBITS-1
call mma_allocate(SSTAB,NTAB,Label='SSTAB')
SSTAB(:) = 0
! The header data
SSTAB(1) = NTAB
SSTAB(2) = ITYPE
SSTAB(3) = -1 ! Not used
SSTAB(4) = nIrrep
SSTAB(5) = NASPRT
SSTAB(6) = MORSBITS
SSTAB(7) = NSSTP
SSTAB(8) = NSBSTOT
SSTAB(9) = KSSTANN
SSTAB(10) = KSSTCRE
SSTAB(11) = KSBSMRS
SSTAB(12) = KMRSSBS
SSTAB(13) = KSBSANN
SSTAB(14) = KSBSCRE
! Fill in the Substring Type table
! Change the counter array into an array of offsets. Also allocate
! a translation table (POP,MS2,ISYM) to Substring Type.
call mma_allocate(SCR2,NSCR,Label='SCR2')
ISSTP = 0
ISBS = 0
do ISPART=1,NASPRT
  NO = NOSUB(ISPART)
  do IPOP=0,NO
    do ISYM=1,nIrrep
      do IMS2=-NO,NO
        ISCR = ISYM+nIrrep*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*(ISPART-1)))
        NSBS = SCR(ISCR)
        SCR(ISCR) = -1
        SCR2(ISCR) = -1
        if (NSBS > 0) then
          ISSTP = ISSTP+1
          SCR(ISCR) = ISBS
          SCR2(ISCR) = ISSTP
          ISBS = ISBS+NSBS
          SSTAB(KSSTP+0+5*(ISSTP-1)) = NSBS
          SSTAB(KSSTP+1+5*(ISSTP-1)) = IPOP
          SSTAB(KSSTP+2+5*(ISSTP-1)) = ISYM
          SSTAB(KSSTP+3+5*(ISSTP-1)) = IMS2
          SSTAB(KSSTP+4+5*(ISSTP-1)) = ISPART
        end if
      end do
    end do
  end do
end do

! Now produce all possible substrings:
!TEST write(u6,*) ' Test in NEWSSTAB, producing substrings.'
!TEST write(u6,'(1x,a,8i5)') 'OSYM(ISORB):',(OSYM(ISORB),ISORB=1,NASPO)
!TEST write(u6,'(1x,a,8i5)') 'OSPN(ISORB):',(OSPN(ISORB),ISORB=1,NASPO)
ISORB = 1
do ISPART=1,NASPRT
  NO = NOSUB(ISPART)
  do MRS=0,2**NO-1
    IPOP = MorsPop(MRS)
    ISYM = MorsSymm(MRS,OSYM(ISORB))
    IMS2 = MorsSpin(MRS,OSPN(ISORB))
    ! Which substring is this?
    ISCR = ISYM+nIrrep*(IPOP+(1+MXPOP)*(MXMS2+IMS2+(2*MXMS2+1)*(ISPART-1)))
    ISBS = 1+SCR(ISCR)
    !TEST write(u6,'(1x,a,8i5)') 'MRS,IPOP,IMS2,ISYM,ISBS:',MRS,IPOP,IMS2,ISYM,ISBS
    SCR(ISCR) = ISBS
    ISSTP = SCR2(ISCR)
    ! Fill in the Substring/Morsel translation arrays:
    LPOS = KSBSMRS+2*(ISBS-1)
    SSTAB(LPOS) = MRS
    SSTAB(LPOS+1) = ISSTP
    LPOS = KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
    SSTAB(LPOS) = ISBS
    SSTAB(LPOS+1) = ISSTP
  end do
  ISORB = ISORB+NO
end do
call mma_deallocate(SCR)
call mma_deallocate(SCR2)
call mma_deallocate(OSPN)
call mma_deallocate(OSYM)
do ISPART=1,NASPRT
  NO = NOSUB(ISPART)
  do MRS=0,2**NO-1
    LPOS = KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
    ISBS = SSTAB(LPOS)
    ISSTP = SSTAB(LPOS+1)
  end do
end do

! Create the Substring Annihilator and Creator arrays, and
! also Substring Type Ann/Cre arrays.
!TEST write(u6,*) ' Making annih and creat arrays:'
do ISPART=1,NASPRT
  NO = NOSUB(ISPART)
  do MRS=0,2**NO-1
    LPOS = KMRSSBS+2*(MRS+(2**MORSBITS)*(ISPART-1))
    ISBS = SSTAB(LPOS)
    ISSTP = SSTAB(LPOS+1)
    do I=1,NO
      NEWMRS = MORSANN(MRS,I)
      if (NEWMRS /= 999999) then
        ISGN = 1
        if (NEWMRS < 0) ISGN = -1
        NEWMRS = ISGN*NEWMRS
        ! Position in Morsel-to-Substring Table
        LPOS = KMRSSBS+2*(NEWMRS+(2**MORSBITS)*(ISPART-1))
        ! New substring times sign factor
        NEWSBS = ISGN*SSTAB(LPOS)
        NEWSST = SSTAB(LPOS+1)
        ! Put new substring in the table
        LPOS = KSBSANN-1+I+MORSBITS*(ISBS-1)
        SSTAB(LPOS) = NEWSBS
        ! Put new substring type in the Substring Type Annihil Table
        LPOS = KSSTANN-1+I+MORSBITS*(ISSTP-1)
        SSTAB(LPOS) = NEWSST
      end if
    end do
    ! Create: Very similar to the Annihilate code above.
    do I=1,NO
      NEWMRS = MORSCRE(MRS,I)
      if (NEWMRS /= 999999) then
        ISGN = 1
        if (NEWMRS < 0) ISGN = -1
        NEWMRS = ISGN*NEWMRS
        ! Position in Morsel-to-Substring Table
        LPOS = KMRSSBS+2*(NEWMRS+(2**MORSBITS)*(ISPART-1))
        ! New substring times sign factor
        NEWSBS = ISGN*SSTAB(LPOS)
        NEWSST = SSTAB(LPOS+1)
        ! Put new substring in the table
        LPOS = KSBSCRE-1+I+MORSBITS*(ISBS-1)
        SSTAB(LPOS) = NEWSBS
        ! Put new substring type in the Substring Type Creat Table
        LPOS = KSSTCRE-1+I+MORSBITS*(ISSTP-1)
        SSTAB(LPOS) = NEWSST
      end if
    end do
  end do
end do

! Test section, may be removed later..
if (.false.) then
  ! For all strings, check the anticommutation rules.
  IERR = 0
  do ISBS=1,NSBSTOT
    ! Its substring type
    LPOS = KSBSMRS+2*(ISBS-1)
    ISSTP = SSTAB(LPOS+1)
    ! Its subpartition
    ISPART = SSTAB(KSSTP+4+5*(ISSTP-1))
    NO = NOSUB(ISPART)
    do I=1,NO
      ! Annihilate orbital nr I in the subpartition.
      LPOS = KSBSANN-1+I+MORSBITS*(ISBS-1)
      ISBS1 = SSTAB(LPOS)
      ! Create orbital nr I in the subpartition.
      LPOS = KSBSCRE-1+I+MORSBITS*(ISBS-1)
      ISBS2 = SSTAB(LPOS)
      if ((ISBS1 /= 0) .and. (ISBS2 /= 0)) IERR = IERR+1
      if ((ISBS1 == 0) .and. (ISBS2 == 0)) IERR = IERR+1
      do J=1,I
        ! Same, for orbital J instead.
        LPOS = KSBSANN-1+J+MORSBITS*(ISBS-1)
        ISBS3 = SSTAB(LPOS)
        LPOS = KSBSCRE-1+J+MORSBITS*(ISBS-1)
        ISBS4 = SSTAB(LPOS)
        ! SBS5=A(I-)*SBS4:
        ISBS5 = 0
        if (ISBS4 /= 0) then
          ISGN = 1
          if (ISBS4 < 0) ISGN = -1
          LPOS = KSBSANN-1+I+MORSBITS*(abs(ISBS4)-1)
          ISBS5 = ISGN*SSTAB(LPOS)
        end if
        ! SBS6=A(J+)*SBS1:
        ISBS6 = 0
        if (ISBS1 /= 0) then
          ISGN = 1
          if (ISBS1 < 0) ISGN = -1
          LPOS = KSBSCRE-1+J+MORSBITS*(abs(ISBS1)-1)
          ISBS6 = ISGN*SSTAB(LPOS)
        end if
        ! SBS7=A(I+)*SBS3:
        ISBS7 = 0
        if (ISBS3 /= 0) then
          ISGN = 1
          if (ISBS3 < 0) ISGN = -1
          LPOS = KSBSCRE-1+I+MORSBITS*(abs(ISBS3)-1)
          ISBS7 = ISGN*SSTAB(LPOS)
        end if
        ! SBS8=A(J-)*SBS2:
        ISBS8 = 0
        if (ISBS2 /= 0) then
          ISGN = 1
          if (ISBS2 < 0) ISGN = -1
          LPOS = KSBSANN-1+J+MORSBITS*(abs(ISBS2)-1)
          ISBS8 = ISGN*SSTAB(LPOS)
        end if
        ! Now check:
        if (I /= J) then
          if (ISBS5+ISBS6 /= 0) IERR = IERR+1
          if (ISBS7+ISBS8 /= 0) IERR = IERR+1
        else
          if (ISBS5+ISBS6 /= ISBS) IERR = IERR+1
          if (ISBS7+ISBS8 /= ISBS) IERR = IERR+1
        end if
      end do
    end do
  end do
  write(u6,*) ' NEWSSTAB: Nr of test errors IERR=',IERR
end if

call mma_deallocate(NOSUB)

end subroutine NEWSSTAB
