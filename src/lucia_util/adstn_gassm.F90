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

!#define _DEBUGPRINT_
subroutine ADSTN_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,IOFFKSP,ICREORB,ICRESTR,IORBTSF,IORBTF,NORBTS,NSTAK,NSTAI,NSTAKTS,NELB, &
                       NACGSOB,ISTMAP,SGNMAP,SCLFAC)
! Creation mappings from K-strings of given sym in each gasspace
!
! Input
! NSTB    : Number of strings before active gasspace
! NSTA    : Number of strings after active gasspace
! IOFFK   : Offset for K group of strings in active gasspace, i.e. start of
!           this symmetry of active K group strings
! IOFFI   : Offset for I group of strings in active gasspace, i.e. start of
!           this symmetry of active I group strings
! IOFFISP : Offset for this symmetrydistribution of active I supergroup strings
! IOFFKSP : Offset for this symmetrydistribution of active K supergroup strings
! ICREORB : Orbital part of creation map for active K groupstrings
! ICRESTR : String  part of creation map for active K groupstrings
! IORBTSF : First active orbital (first orbital in in active GASspace
!           with required sym)
! IORBTF  : First orbital in active gas space, (can have any sym)
! NORBTS  : Number of orbitals of given symmetry and type
! NSTAK   : Number of K groupstrings with given correct symmetry
! NSTAKTS : Total Number of K supergroup strings with correct symmetry
! NSTAI   : Number of I groupstrings in active gasspace

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSTB, NSTA, IOFFK, IOFFI, IOFFISP, IOFFKSP, NACGSOB, NSTAK, ICREORB(NACGSOB,NSTAK+IOFFK-1), &
                                 ICRESTR(NACGSOB,NSTAK+IOFFK-1), IORBTSF, IORBTF, NORBTS, NSTAI, NSTAKTS, NELB
integer(kind=iwp), intent(inout) :: ISTMAP(NSTAKTS,IORBTSF+NORBTS)
real(kind=wp), intent(inout) :: SGNMAP(NSTAKTS,IORBTSF+NORBTS)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: IA, IADRI0, IADRK0, IB, IORB, IORBRT, IORBRTS, ISTR, KSTR, NSTAINSTA, NSTAKNSTA
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IORBR, NK
#endif
real(kind=wp) :: SIGN0, SGN

!write(u6,*) ' ADSTN_GASSM : NSTA, NSTB, NSTAK',NSTA,NSTB,NSTAK
!write(u6,*) ' IOFFISP,IOFFKSP',IOFFISP,IOFFKSP
!write(u6,*) ' IORBTSF IORBTF ',IORBTSF,IORBTF

!SIGN0 = SCLFAC*(-One)**NELB
if (mod(NELB,2) == 0) then
  SIGN0 = SCLFAC
else
  SIGN0 = -SCLFAC
end if
!write(u6,*) ' NELB sign0 = ',NELB,SIGN0
do KSTR=IOFFK,NSTAK+IOFFK-1
  do IORB=IORBTSF,IORBTSF-1+NORBTS
    ! Relative to Type-symmetry start
    IORBRTS = IORB-IORBTSF+1
    ! Relative to type start
    IORBRT = IORB-IORBTF+1
    !write(u6,*) 'IORB IORBRT KSTR ',IORB,IORBRT,KSTR
    !write(u6,*) 'ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)',ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)
    if (ICREORB(IORBRT,KSTR) > 0) then
      ! Excitation is open, corresponding active I string
      if (ICRESTR(IORBRT,KSTR) > 0) then
        SGN = SIGN0
        ISTR = ICRESTR(IORBRT,KSTR)
      else
        SGN = -SIGN0
        ISTR = -ICRESTR(IORBRT,KSTR)
      end if
      ! Relative to start of given symmetry for this group
      ISTR = ISTR-IOFFI+1
      ! This Creation is active for all choices of strings in supergroup
      ! before and after the active type. Store the corrsponding mappings
      IADRK0 = (KSTR-IOFFK)*NSTA+IOFFKSP-1
      IADRI0 = (ISTR-1)*NSTA+IOFFISP-1
      !write(u6,*) ' ISTR IADRK0 IADRI0 = ',ISTR,IADRK0,IADRI0

      NSTAINSTA = NSTAI*NSTA
      NSTAKNSTA = NSTAK*NSTA
      do IB=1,NSTB
        ISTMAP(IADRK0+1:IADRK0+NSTA,IORBRTS) = [(IADRI0+IA,IA=1,NSTA)]
        SGNMAP(IADRK0+1:IADRK0+NSTA,IORBRTS) = SGN
        IADRI0 = IADRI0+NSTAINSTA
        IADRK0 = IADRK0+NSTAKNSTA
      end do
    !else
    !  SGN = Zero
    !  ISTR = 0
    !  ! This Creation is inactive for all choices of strings in supergroup
    !  ! before and after the active type.
    !  IADRK0 = (KSTR-IOFFK)*NSTA+IOFFKSP-1
    !  !write(u6,*) ' ISTR IADRK0 = ',ISTR,IADRK0
    !
    !  do IB=1,NSTB
    !    do IA=1,NSTA
    !      KBKA = IADRK0+(IB-1)*NSTAK*NSTA+IA
    !      !write(u6,*) ' IBKA, KBKA ',IBKA,KBKA
    !      if (ISTMAP(KBKA,IORBRTS) /= 999) then
    !        write(u6,*) ' overwriting ???'
    !        write(u6,*) ' Element ',(IORBRTS-1)*NSTAKTS+KBKA
    !        stop
    !      end if
    !      ISTMAP(KBKA,IORBRTS) = ISTR
    !      SGNMAP(KBKA,IORBRTS) = SGN
    !
    !    end do
    !  end do
    end if
    ! This Creation is active for all choices of strings in supergroup
    ! before and after the active type. Store the corrsponding mappings
    !OLD IADRK0 = (KSTR-IOFFK)*NSTA+IOFFKSP-1
    !OLD IADRI0 = (ISTR-1)*NSTA+IOFFISP-1
    !write(u6,*) ' ISTR IADRK0 IADRI0 = ',ISTR,IADRK0,IADRI0

    !OLD do IB=1,NSTB
    !OLD  do IA=1,NSTA
    !OLD    IBKA = IADRI0+(IB-1)*NSTAI*NSTA+IA
    !OLD    KBKA = IADRK0+(IB-1)*NSTAK*NSTA+IA
    !?      write(u6,*) ' IBKA, KBKA ',IBKA,KBKA
    !OLD    ISTMAP(KBKA,IORBRTS) = IBKA
    !OLD    SGNMAP(KBKA,IORBRTS) = SGN
    !OLD  end do
    !OLD end do
    !end if

  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from ADSTN_GASSM'
write(u6,*) ' ======================='
NK = NSTB*NSTAK*NSTA
write(u6,*) ' Number of K strings accessed ',NK
if (NK /= 0) then
  do IORB=IORBTSF,IORBTSF+NORBTS-1
    IORBR = IORB-IORBTSF+1
    write(u6,*) ' Update Info for orbital ',IORB
    write(u6,*) ' Excited strings and sign'
    call IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
    call WRTMAT(SGNMAP(1,IORBR),1,NK,1,NK)
  end do
end if
#endif

end subroutine ADSTN_GASSM
