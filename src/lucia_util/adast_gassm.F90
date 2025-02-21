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
subroutine ADAST_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,IOFFKSP,ICREORB,ICRESTR,IORBTSF,IORBTF,NORBTS,NSTAK,NSTAI,NSTAKTS,NELB, &
                       ISTMAP,SGNMAP,SCLFAC,IAC,LROW_IN,IEC)
! Annihilation/Creation mappings from K-strings of given sym in each gasspace
!
! Input
! NSTB    : Number of strings before active gasspace
! NSTA    : Number of strings after accive gasspace
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

use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSTB, NSTA, IOFFK, IOFFI, IOFFISP, IOFFKSP, LROW_IN, NSTAK, ICREORB(LROW_IN,NSTAK+IOFFK-1), &
                                 ICRESTR(LROW_IN,NSTAK+IOFFK-1), IORBTSF, IORBTF, NORBTS, NSTAI, NSTAKTS, NELB, IAC, IEC
integer(kind=iwp), intent(inout) :: ISTMAP(NSTAKTS,NORBTS)
real(kind=wp), intent(inout) :: SGNMAP(NSTAKTS,NORBTS)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: I_AM_ACTIVE, IA, IADRI0, IADRK0, IB, IORB, IORBRT, IORBRTS, IROW, ISTR, KSTR, NSTAINSTA, NSTAKNSTA
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IORBR, NK
#endif
real(kind=wp) :: SIGN0, SGN

! Some dummy initializations
SGN = Zero ! jwk-cleanup
ISTR = 0 ! jwk-cleanup

!write(u6,*) ' ICRESTR'
!call IWRTMA(ICRESTR,LROW_IN,NSTAK,LROW_IN,NSTAK)
!write(u6,*) ' IOFFI = ',IOFFI
!PAM2009 SIGN0 = SCLFAC*(-One)**NELB
if (mod(NELB,2) == 0) then
  SIGN0 = SCLFAC
else
  SIGN0 = -SCLFAC
end if
do KSTR=IOFFK,NSTAK+IOFFK-1
  do IORB=IORBTSF,IORBTSF-1+NORBTS
    ! Relative to Type-symmetry start
    IORBRTS = IORB-IORBTSF+1
    ! Relative to type start
    IORBRT = IORB-IORBTF+1
    !write(u6,*) ' IORBRTS IORBRT',IORBRTS,IORBRT
    ! Change of active group
    I_AM_ACTIVE = 0
    if (IAC == 2) then
      !write(u6,*) ' ICREORB = ',ICREORB(IORBRT,KSTR)
      !write(u6,*) ' ICRESTR = ',ICRESTR(IORBRT,KSTR)
      if (ICREORB(IORBRT,KSTR) > 0) then
        ! Creation is nonvanishing
        I_AM_ACTIVE = 1
        if (ICRESTR(IORBRT,KSTR) > 0) then
          SGN = SIGN0
          ISTR = ICRESTR(IORBRT,KSTR)
        else
          SGN = -SIGN0
          ISTR = -ICRESTR(IORBRT,KSTR)
        end if
      end if
    else if (IAC == 1) then
      if (IEC == 1) then
        ! Expanded map
        if (ICREORB(IORBRT,KSTR) < 0) then
          ! Annihilation is non-vanishing
          I_AM_ACTIVE = 1
          if (ICRESTR(IORBRT,KSTR) > 0) then
            SGN = SIGN0
            ISTR = ICRESTR(IORBRT,KSTR)
          else
            SGN = -SIGN0
            ISTR = -ICRESTR(IORBRT,KSTR)
          end if
        end if
      else
        ! Compressed map
        do IROW=1,LROW_IN
          !write(u6,*) ' IROW, ICREORB(IROW,KSTR)',IROW,ICREORB(IROW,KSTR)
          !OLD if (ICREORB(IROW,KSTR)  == -IORBRT) then
          if (ICREORB(IROW,KSTR) == -IORB) then
            ! Annihilation is non-vanishing
            I_AM_ACTIVE = 1
            if (ICRESTR(IROW,KSTR) > 0) then
              SGN = SIGN0
              ISTR = ICRESTR(IROW,KSTR)
            else
              SGN = -SIGN0
              ISTR = -ICRESTR(IROW,KSTR)
            end if
          end if
        end do
      end if
      ! End of expanded/compact switch
    end if
    ! End of Creation/annihilation switch

    if (I_AM_ACTIVE == 1) then
      ! Excitation is open, corresponding active I string
      ! Relative to start of given symmetry for this group
      ISTR = ISTR-IOFFI+1
      !write(u6,*) ' ISTR, relative = ',ISTR
      ! This Creation is active for all choices of strings in supergroup
      ! before and after the active type. Store the corrsponding mappings
      IADRK0 = (KSTR-IOFFK)*NSTA+IOFFKSP-1
      IADRI0 = (ISTR-1)*NSTA+IOFFISP-1
      !write(u6,*) ' IADRK0 IOFFK IOFFKSP ',IADRK0,IOFFK,IOFFKSP
      !write(u6,*) ' IADRI0, IOFFISP ',IADRI0,IOFFISP

      NSTAINSTA = NSTAI*NSTA
      NSTAKNSTA = NSTAK*NSTA

      !write(u6,*) ' ISTR NSTA NSTB ',ISTR,NSTA,NSTB
      !write(u6,*) ' NSTAI,NSTAK',NSTAI,NSTAK
      do IB=1,NSTB
        ISTMAP(IADRK0+1:IADRK0+NSTA,IORBRTS) = [(IADRI0+IA,IA=1,NSTA)]
        SGNMAP(IADRK0+1:IADRK0+NSTA,IORBRTS) = SGN
        IADRI0 = IADRI0+NSTAINSTA
        IADRK0 = IADRK0+NSTAKNSTA
      end do
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from ADAST_GASSM'
write(u6,*) ' ======================='
NK = NSTB*NSTAK*NSTA
write(u6,*) ' Number of K strings accessed ',NK
if (NK /= 0) then
  do IORB=IORBTSF,IORBTSF+NORBTS-1
    IORBR = IORB-IORBTSF+1
    write(u6,*) ' Update Info for orbital ',IORB
    write(u6,*) ' Mapped strings and sign'
    call IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
    call WRTMAT(SGNMAP(1,IORBR),1,NK,1,NK)
  end do
end if
#endif

end subroutine ADAST_GASSM
