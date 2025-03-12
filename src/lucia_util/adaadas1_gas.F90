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
subroutine ADAADAS1_GAS(NK,I1,XI1S,LI1,IORB,NIORB,IAC,JORB,NJORB,JAC,KSTR,NKEL,NKSTR,IREO,IZ,NOCOB,KMAX,KMIN,IEND,SCLFAC,NSTRI)
! Obtain I1(KSTR) = +/- a+/a  IORB a+/a JORB !KSTR>
! Only orbital pairs IOB > JOB are included
!
! KSTR is restricted to strings with relative numbers in the
! range KMAX to KMIN
!
! =====
! Input
! =====
! IORB : First I orbital to be added
! NIORB : Number of orbitals to be added : IORB to IORB-1+NIORB
!        are used. They must all be in the same TS group
! JORB : First J orbital to be added
! LORB : Number of orbitals to be added : JORB to JORB-1+NJORB
!        are used. They must all be in the same TS group
! KMAX : Largest allowed relative number for K strings
! KMIN : Smallest allowed relative number for K strings
!
! ======
! Output
! ======
!
! NK      : Number of K strings
! I1(KSTR,JORB) : /= 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
! XI1S(KSTR,JORB) : above +/-
!               : == 0    a + JORB !KSTR> = 0
! Offset is KMIN
!
! ISCR : Local scratch, at most 1000 orbitals in a given TS block)
!
! L.R. Jan 20, 1998

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LI1, IORB, NIORB, IAC, JORB, NJORB, JAC, NKEL, NKSTR, KSTR(NKEL,NKSTR), IREO(*), NOCOB, &
                                 IZ(NOCOB,NKEL+2), KMAX, KMIN, NSTRI
integer(kind=iwp), intent(out) :: NK, I1(LI1,NIORB*NJORB), IEND
real(kind=wp), intent(out) :: XI1S(LI1,NIORB*NJORB)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: IACT, IDIAG, IEL, IIEL, IIELMAX, IIELMIN, IIORB, IJ, IJOFF, ILEX, ILEX0, ILEX1, ILEX2, IORB1, IORB2, IORBMAX, &
                     IORBMIN, ISCR(1000), JEL, JJ, JJEL, JJELMAX, JJELMIN, JJORB, JORB1, JORB2, JORBMAX, JORBMIN, KEL, KEND, &
                     KKORB, KKSTR
real(kind=wp) :: SIGNIJ, SIGNJ, SGN, SGNARR(0:63)

!PAM2009 Array of values, replacing expressions such as (-One)**INT:
SGNARR(0::2) = One
SGNARR(1::2) = -One

! Some dummy initializations
SGN = Zero

#ifdef _DEBUGPRINT_
write(u6,*) ' ===================='
write(u6,*) ' ADADS1_GAS speaking'
write(u6,*) ' ===================='
write(u6,*) ' IORB,NIORB,IAC ',IORB,NIORB,IAC
write(u6,*) ' JORB,NJORB,JAC ',JORB,NJORB,JAC

!write(u6,*) ' Kstrings in action (el,string)'
!write(u6,*) ' ==============================='
write(u6,*) ' 24 elements of reorder array'
call IWRTMA(IREO,1,24,1,24)
!    IWRTMA(KSTR,NKEL,NKSTR,NKEL,NKSTR)
#endif

IORBMIN = IORB
IORBMAX = IORB+NIORB-1

JORBMIN = JORB
JORBMAX = JORB+NJORB-1

KEND = min(NKSTR,KMAX)
if (KEND < NKSTR) then
  IEND = 0
else
  IEND = 1
end if
NK = KEND-KMIN+1

if ((IAC == 2) .and. (JAC == 2)) then

  ! ==========================
  ! Creation- creation mapping
  ! ==========================

  do KKSTR=KMIN,KEND
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Occupation of string ',KKSTR
    call IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
#   endif
    ! Loop over electrons after which JORB can be added
    do JEL=0,NKEL

      if (JEL == 0) then
        JORB1 = JORBMIN-1
      else
        JORB1 = max(JORBMIN-1,KSTR(JEL,KKSTR))
      end if
      if (JEL == NKEL) then
        JORB2 = JORBMAX+1
      else
        JORB2 = min(JORBMAX+1,KSTR(JEL+1,KKSTR))
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
#     endif

      if ((JEL > 0) .and. (JORB1 >= JORBMIN) .and. (JORB1 <= JORBMAX)) then
        ! vanishing for any IORB
        IJOFF = (JORB1-JORBMIN)*NIORB
        do IIORB=1,NIORB
          IJ = IJOFF+IIORB
          I1(KKSTR-KMIN+1,IJ) = 0
          XI1S(KKSTR-KMIN+1,IJ) = Zero
        end do
      end if

      if ((JORB1 < JORBMAX) .and. (JORB2 > JORBMIN)) then
        ! Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
        !PAM2009 SIGNJ = SCLFAC*(-One)**JEL)
        SIGNJ = SGNARR(JEL)*SCLFAC
        ! reverse lexical number of the first JEL ELECTRONS
        ILEX0 = 1
        do JJEL=1,JEL
          ILEX0 = ILEX0+IZ(KSTR(JJEL,KKSTR),JJEL)
        end do
        do JJORB=JORB1+1,JORB2-1
          ! And electron JEL + 1
          ILEX1 = ILEX0+IZ(JJORB,JEL+1)
          ! Add electron IORB
          do IEL=JEL,NKEL
            if (IEL == JEL) then
              IORB1 = max(JJORB,IORBMIN-1)
            else
              IORB1 = max(IORBMIN-1,KSTR(IEL,KKSTR))
            end if
            if (IEL == NKEL) then
              IORB2 = IORBMAX+1
            else
              IORB2 = min(IORBMAX+1,KSTR(IEL+1,KKSTR))
            end if
#           ifdef _DEBUGPRINT_
            write(u6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
#           endif
            if ((IEL > JEL) .and. (IORB1 >= IORBMIN) .and. (IORB1 <= IORBMAX)) then
              IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
              I1(KKSTR-KMIN+1,IJ) = 0
              XI1S(KKSTR-KMIN+1,IJ) = Zero
            end if
            if ((IORB1 < IORBMAX) .and. (IORB2 > IORBMIN)) then
              ! Orbitals IORB1+1 - IORB2 -1 can be added after ELECTRON IEL in KSTR
              ! Reverse lexical number of the first IEL+1 electrons
              ILEX2 = ILEX1
              do IIEL=JEL+1,IEL
                ILEX2 = ILEX2+IZ(KSTR(IIEL,KKSTR),IIEL+1)
              end do
              ! add terms for the last electrons
              do IIEL=IEL+1,NKEL
                ILEX2 = ILEX2+IZ(KSTR(IIEL,KKSTR),IIEL+2)
              end do
              IJOFF = (JJORB-JORBMIN)*NIORB
              !PAM2009 SIGNIJ = SIGNJ*(-One)**(IEL+1)
              SIGNIJ = SIGNJ*SGNARR(IEL+1)
              do IIORB=IORB1+1,IORB2-1
                IJ = IJOFF+IIORB-IORBMIN+1
                ILEX = ILEX2+IZ(IIORB,IEL+2)
                IACT = IREO(ILEX)
#               ifdef _DEBUGPRINT_
                write(u6,*) 'IIORB JJORB',IIORB,JJORB
                write(u6,*) ' ILEX IACT ',ILEX,IACT
#               endif
                I1(KKSTR-KMIN+1,IJ) = IACT
                XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
              end do
            end if
          end do
        end do
      end if
    end do
  end do
else if ((IAC == 1) .and. (JAC == 1)) then

  ! =============================================
  ! annihilation - annihilation mapping (a i a j)
  ! =============================================

  do KKSTR=KMIN,KEND
    ! Active range for electrons
    IIELMIN = 0
    IIELMAX = 0
    JJELMIN = 0
    JJELMAX = 0

    do KEL=1,NKEL
      KKORB = KSTR(KEL,KKSTR)
      if ((IIELMIN == 0) .and. (KKORB >= IORBMIN)) IIELMIN = KEL
      if ((JJELMIN == 0) .and. (KKORB >= JORBMIN)) JJELMIN = KEL
      if (KKORB <= IORBMAX) IIELMAX = KEL
      if (KKORB <= JORBMAX) JJELMAX = KEL
    end do
    if (IIELMIN == 0) IIELMIN = NKEL+1
    if (JJELMIN == 0) JJELMIN = NKEL+1

#   ifdef _DEBUGPRINT_
    write(u6,*) ' Occupation of string ',KKSTR
    call IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
#   endif
    ! Loop over first electron to be removed
    !do JEL=1,NKEL
    do JEL=JJELMIN,JJELMAX
      JJORB = KSTR(JEL,KKSTR)
      ! Loop over second electron to be removed
      !do IEL=JEL+1,NKEL
      do IEL=max(JEL+1,IIELMIN),IIELMAX
        IIORB = KSTR(IEL,KKSTR)
#       ifdef _DEBUGPRINT_
        write(u6,*) ' IEL JEL IORB JORB ',IEL,JEL,IORB,JORB
#       endif
        if ((IIORB >= IORBMIN) .and. (IIORB <= IORBMAX) .and. (JJORB >= JORBMIN) .and. (JJORB <= JORBMAX)) then
          !PAM09 SGN = SCLFAC*(-One)**(IEL+JEL-1)
          SGN = SGNARR(IEL+JEL+1)*SCLFAC
          ! reverse lexical number of the double annihilated string
          ILEX = 1
          do JJEL=1,JEL-1
            ILEX = ILEX+IZ(KSTR(JJEL,KKSTR),JJEL)
          end do
          do JJEL=JEL+1,IEL-1
            ILEX = ILEX+IZ(KSTR(JJEL,KKSTR),JJEL-1)
          end do
          do JJEL=IEL+1,NKEL
            ILEX = ILEX+IZ(KSTR(JJEL,KKSTR),JJEL-2)
          end do
          IACT = IREO(ILEX)
          if ((IACT <= 0) .or. (IACT > NSTRI)) then
            write(u6,*) ' IACT out of bounds, IACT =  ',IACT
            !stop ' IACT out of bounds'
            call SYSABENDMSG('lucia_util/adaadas1_gas','Internal error','')
          end if
#         ifdef _DEBUGPRINT_
          write(u6,*) ' ILEX and IACT ',ILEX,IACT
#         endif

          IJ = (JJORB-JORB)*NIORB+IIORB-IORB+1
          I1(KKSTR-KMIN+1,IJ) = IACT
          XI1S(KKSTR-KMIN+1,IJ) = SGN
        end if
        ! End if orbitals are in correct range
      end do
      ! End of loop over IEL
    end do
    ! End of loop over JEL
  end do
  ! End of loop over Kstrings

else if ((IAC == 2) .and. (JAC == 1)) then

  ! ==================================
  ! Creation-annihilation map a+ i a j
  ! ==================================

  !do KKSTR=1,NKSTR
  do KKSTR=KMIN,KEND
    ! Indicate where a given orbital i should be added in KKSTR
    ISCR(IORBMIN:IORBMIN+NIORB-1) = 0
    IIEL = 1
    do IIORB=IORBMIN,IORBMAX
      do
        if (IIEL <= NKEL) then
          if (IIORB < KSTR(IIEL,KKSTR)) then
            ISCR(IIORB) = IIEL
            exit
          else if (IIORB == KSTR(IIEL,KKSTR)) then
            ISCR(IIORB) = -IIEL
            IIEL = IIEL+1
            exit
          else if (IIORB > KSTR(IIEL,KKSTR)) then
            IIEL = IIEL+1
          end if
        else if (IIEL == NKEL+1) then
          ISCR(IIORB) = IIEL
          exit
        end if
      end do
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ISCR from IORBMIN array for KKSTR = ',KKSTR
    write(u6,*) ' IORBMIN, NIORB',IORBMIN,NIORB
    call IWRTMA(ISCR(IORBMIN),1,NIORB,1,NIORB)
#   endif
    do JEL=1,NKEL
      JJORB = KSTR(JEL,KKSTR)
      if ((JJORB >= JORBMIN) .and. (JJORB <= JORBMAX)) then
        do IIORB=IORBMIN,IORBMAX
          IEL = ISCR(IIORB)
          !write(u6,*) ' JEL IEL JJORB IIORB',JEL,IEL,JJORB,IIORB
          IACT = 0
          if ((IEL > 0) .and. (IIORB > JJORB)) then
            ! New string is  a+1 ... a+ jel-1 a+jel+1 ..a+iel-1 a+iiorb a+iel+1 ...
            ! Lexical number of new string
            ILEX = 1
            do KEL=1,JEL-1
              ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
            end do
            do KEL=JEL+1,IEL-1
              ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL-1)
            end do
            ILEX = ILEX+IZ(IIORB,IEL-1)
            do KEL=IEL,NKEL
              ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
            end do
            IACT = IREO(ILEX)
            if ((IACT <= 0) .or. (IACT > NSTRI)) then
              write(u6,*) ' 1: IACT out of bounds, IACT =  ',IACT
              write(u6,*) ' NSTRI = ',NSTRI
              write(u6,*) 'IIORB,JJORB ',IIORB,JJORB
              write(u6,*) ' Kstring :'
              call IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
              write(u6,*) ' ILEX = ',ILEX
              write(u6,*) 'IZ matrix'
              call IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
              !stop ' IACT out of bounds'
              call SYSABENDMSG('lucia_util/adaadas1_gas','Internal error','')
            end if
            !PAM2009 SGN = SCLFAC*(-One)**(IEL+JEL-1)
            SGN = SGNARR(IEL+JEL+1)*SCLFAC
          else if ((IEL > 0) .and. (IIORB < JJORB)) then
            ! New string is  a+1 ... a+ iel-1 a+ iiorb a+iel+1 ..a+jel-1 a+jel+1 ...
            ILEX = 1
            do KEL=1,IEL-1
              ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
            end do
            ILEX = ILEX+IZ(IIORB,IEL)
            do KEL=IEL,JEL-1
              ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL+1)
            end do
            do KEL=JEL+1,NKEL
              ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
            end do
            !write(u6,*) ' ILEX =',ILEX
            IACT = IREO(ILEX)
            if ((IACT <= 0) .or. (IACT > NSTRI)) then
              write(u6,*) '2 IACT out of bounds, IACT =  ',IACT
              !stop ' IACT out of bounds'
              call SYSABENDMSG('lucia_util/adaadas1_gas','Internal error','')
            end if
            !write(u6,*) ' IACT = ',IACT
            !PAM2009 SGN = SCLFAC*(-One)**(IEL+JEL)
            SGN = SGNARR(IEL+JEL)*SCLFAC
          else if ((IEL < 0) .and. (IIORB == JJORB)) then
            ! Diagonal excitation
            SGN = SCLFAC
            IACT = KKSTR
          end if
          if (IACT /= 0) then
            IJ = (JJORB-JORB)*NIORB+IIORB-IORB+1
            I1(KKSTR-KMIN+1,IJ) = IACT
            XI1S(KKSTR-KMIN+1,IJ) = SGN
          end if
        end do
        ! End of loop over IIORB
      end if
      ! End of  active cases
    end do
    ! End of loop over electrons to be annihilated
  end do
  !^ End of loop over Kstrings
else if ((IAC == 1) .and. (JAC == 2)) then

  ! ===================================
  ! Annihilation-creation  map a i a+ j
  ! ===================================

  ! Diagonal excitations ?
  if (IORBMIN == JORBMIN) then
    IDIAG = 1
  else
    IDIAG = 0
  end if
  !do KKSTR=1,NKSTR
  do KKSTR=KMIN,KEND
    ! Indicate where a given orbital j should be added in KKSTR
    ISCR(JORBMIN:JORBMIN+NJORB-1) = 0
    JJEL = 1
    do JJORB=JORBMIN,JORBMAX
      do
        if (JJEL <= NKEL) then
          if (JJORB < KSTR(JJEL,KKSTR)) then
            ISCR(JJORB) = JJEL
            exit
          else if (JJORB == KSTR(JJEL,KKSTR)) then
            ISCR(JJORB) = -JJEL
            JJEL = JJEL+1
            exit
          else if (JJORB > KSTR(JJEL,KKSTR)) then
            JJEL = JJEL+1
          end if
        else if (JJEL == NKEL+1) then
          ISCR(JJORB) = JJEL
          exit
        end if
      end do
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) ' ISCR from JORBMIN array for KKSTR = ',KKSTR
    call IWRTMA(ISCR(JORBMIN),1,NJORB,1,NJORB)
#   endif
    do IEL=1,NKEL+IDIAG
      ! IEL = NKEL + 1 will be used for excitations a+j aj
      if (IEL <= NKEL) IIORB = KSTR(IEL,KKSTR)
      if (((IIORB >= IORBMIN) .and. (IIORB <= IORBMAX)) .or. (IEL == NKEL+1)) then
        do JJORB=JORBMIN,JORBMAX
          JEL = ISCR(JJORB)
          !write(u6,*) ' IEL IIORB JEL JJORB ',IEL,IIORB,JEL,JJORB
          IACT = 0
          if (IEL <= NKEL) then
            if ((JEL > 0) .and. (JJORB > IIORB)) then
              ! New string is  a+1 ... a+ iel-1 a+iel+1 ..a+jel-1 a+jjorb a+jel+1 ...
              ! Lexical number of new string
              ILEX = 1
              do KEL=1,IEL-1
                ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
              end do
              do KEL=IEL+1,JEL-1
                ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL-1)
              end do
              ILEX = ILEX+IZ(JJORB,JEL-1)
              do KEL=JEL,NKEL
                ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
              end do
              IACT = IREO(ILEX)
              if ((IACT <= 0) .or. (IACT > NSTRI)) then
                write(u6,*) '3 IACT out of bounds, IACT =  ',IACT
                write(u6,*) ' ILEX,KKSTR ',ILEX,KKSTR
                write(u6,*) ' occupation of KSTR'
                call IWRTMA(KSTR,1,NKEL,1,NKEL)
                write(u6,*) ' IEL JEL ',IEL,JEL
                write(u6,*) ' IIORB,JJORB',IIORB,JJORB
                !stop ' IACT out of bounds'
                call SYSABENDMSG('lucia_util/adaadas1_gas','Internal error','')
              end if
              !PAM2009 SGN = SCLFAC*(-One)**(IEL+JEL)
              SGN = SGNARR(IEL+JEL)*SCLFAC
            else if ((JEL > 0) .and. (JJORB < IIORB)) then
              ! New string is  a+1 ... a+ jel-1 a+ jjorb a+jel+1 ..a+iel-1 a+iel+1 ...
              ILEX = 1
              do KEL=1,JEL-1
                ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
              end do
              ILEX = ILEX+IZ(JJORB,JEL)
              do KEL=JEL,IEL-1
                ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL+1)
              end do
              do KEL=IEL+1,NKEL
                ILEX = ILEX+IZ(KSTR(KEL,KKSTR),KEL)
              end do
              IACT = IREO(ILEX)
              !PAM2009 SGN = SCLFAC*(-One)**(IEL+JEL-1)
              SGN = SGNARR(IEL+JEL+1)*SCLFAC
              if ((IACT <= 0) .or. (IACT > NSTRI)) then
                write(u6,*) '4 IACT out of bounds, IACT =  ',IACT
                write(u6,*) ' NSTRI = ',NSTRI
                write(u6,*) 'IIORB,JJORB ',IIORB,JJORB
                write(u6,*) ' Kstring :'
                call IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
                write(u6,*) ' ILEX = ',ILEX
                write(u6,*) 'IZ matrix'
                call IWRTMA(IZ,NOCOB,NKEL,NOCOB,NKEL)
                !stop ' IACT out of bounds'
                call SYSABENDMSG('lucia_util/adaadas1_gas','Internal error','')
              end if
            end if
            if (IACT /= 0) then
              IJ = (JJORB-JORB)*NIORB+IIORB-IORB+1
              I1(KKSTR-KMIN+1,IJ) = IACT
              XI1S(KKSTR-KMIN+1,IJ) = SGN
            end if
          else if ((IEL == NKEL+1) .and. (JEL > 0)) then
            ! Diagonal excitations aja+j
            JJ = (JJORB-JORB)*NJORB+JJORB-JORB+1
            I1(KKSTR-KMIN+1,JJ) = KKSTR
            XI1S(KKSTR-KMIN+1,JJ) = SCLFAC
          end if
        end do
        ! End of loop over JJORB
      end if
      ! End of  active cases
    end do
    ! End of loop over electrons to be annihilated
  end do
  ! End of loop over Kstrings
end if
! End of types of creation mappings

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from ADADST1_GAS'
write(u6,*) ' ====================='
write(u6,*) ' Number of K strings accessed ',NK
if (NK /= 0) then
  IJ = 0
  do JJORB=JORB,JORB+NJORB-1
    do IIORB=IORB,IORB+NIORB-1
      IJ = IJ+1
      !write(u6,*) ' IJ = ',IJ
      !if (IIORB > JJORB) then
      write(u6,*) ' Info for orbitals (iorb,jorb) ',IIORB,JJORB
      write(u6,*) ' Excited strings and sign'
      call IWRTMA(I1(1,IJ),1,NK,1,NK)
      call WRTMAT(XI1S(1,IJ),1,NK,1,NK)
      !end if
    end do
  end do
end if
#endif

end subroutine ADAADAS1_GAS
