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
subroutine ADADS1_GAS(NK,I1,XI1S,LI1,IORB,NIORB,JORB,NJORB,KSTR,NKEL,NKSTR,IREO,IZ,NOCOB,KMAX,KMIN,IEND,SCLFAC)
! Obtain I1(KSTR) = +/- A+ IORB A+ JORB !KSTR>
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
! I1(KSTR,JORB) : ne. 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
! XI1S(KSTR,JORB) : above +/-
!          : eq. 0    a + JORB !KSTR> = 0
! Offset is KMIN

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LI1, IORB, NIORB, JORB, NJORB, NKEL, NKSTR, KSTR(NKEL,NKSTR), IREO(*), NOCOB, IZ(NOCOB,NKEL+2), &
                                 KMAX, KMIN
integer(kind=iwp), intent(out) :: NK, I1(LI1,NIORB*NJORB), IEND
real(kind=wp), intent(out) :: XI1S(LI1,NIORB*NJORB)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: IACT, IEL, IIEL, IIORB, IJ, IJOFF, ILEX, ILEX0, ILEX1, ILEX2, IORB1, IORB2, IORBMAX, IORBMIN, JEL, JJEL, &
                     JJORB, JORB1, JORB2, JORBMAX, JORBMIN, KEND, KKSTR, NIJ
real(kind=wp) :: ODDIEL, ODDJEL, SIGNIJ, SIGNJ

#ifdef _DEBUGPRINT_
write(u6,*) ' ===================='
write(u6,*) ' ADADS1_GAS speaking'
write(u6,*) ' ===================='
write(u6,*) ' IORB,NIORB ',IORB,NIORB
write(u6,*) ' JORB,NJORB ',JORB,NJORB
#endif

IORBMIN = IORB
IORBMAX = IORB+NIORB-1

JORBMIN = JORB
JORBMAX = JORB+NJORB-1

NIJ = NIORB*NJORB

KEND = min(NKSTR,KMAX)
if (KEND < NKSTR) then
  IEND = 0
else
  IEND = 1
end if
NK = KEND-KMIN+1

do KKSTR=KMIN,KEND
# ifdef _DEBUGPRINT_
  write(u6,*) ' Occupation of string ',KKSTR
  call IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
# endif
  ! Loop over electrons after which JORB can be added
  !PAM2009 Added variable ODDJEL to replace (-One)**JEL
  ODDJEL = -One
  do JEL=0,NKEL
    ODDJEL = -ODDJEL

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
#   ifdef _DEBUGPRINT_
    write(u6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
#   endif

    if ((JEL > 0) .and. (JORB1 >= JORBMIN) .and. (JORB1 <= JORBMAX)) then
      ! vanishing for any IORB
      IJOFF = (JORB1-JORBMIN)*NIORB
      do IIORB=1,NIORB
        IJ = IJOFF+IIORB
        if (ij > nij) then
          write(u6,*) ' ij > nij'
          write(u6,*) ' JORB1 IIORB',JORB1,IIORB
          write(u6,*) ' ijoff ',ijoff
          !stop
          call SYSABENDMSG('lucia_util/adads1_gas','Internal error','')
        end if
        I1(KKSTR-KMIN+1,IJ) = 0
        XI1S(KKSTR-KMIN+1,IJ) = Zero
      end do
    end if

    if ((JORB1 < JORBMAX) .and. (JORB2 > JORBMIN)) then
      ! Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
      !PAM2009 SIGNJ = SCLFAC*(-One)**JEL
      SIGNJ = ODDJEL*SCLFAC
      ! reverse lexical number of the first JEL ELECTRONS
      ILEX0 = 1
      do JJEL=1,JEL
        ILEX0 = ILEX0+IZ(KSTR(JJEL,KKSTR),JJEL)
      end do
      do JJORB=JORB1+1,JORB2-1
        ! And electron JEL + 1
        ILEX1 = ILEX0+IZ(JJORB,JEL+1)
        ! Add electron IORB
        !PAM2009 Added ODDIEL to replace (-One)**IEL
        ODDIEL = -ODDJEL
        do IEL=JEL,NKEL
          ODDIEL = -ODDIEL
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
#         ifdef _DEBUGPRINT_
          write(u6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
#         endif
          if ((IEL > JEL) .and. (IORB1 >= IORBMIN) .and. (IORB1 <= IORBMAX)) then
            IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
            if (ij > nij) then
              write(u6,*) ' ij > nij'
              write(u6,*) ' JJORB IORB1',JJORB,IORB1
              write(u6,*) ' ijoff ',ijoff
              !stop
              call SYSABENDMSG('lucia_util/adads1_gas','Internal error','')
            end if
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
            SIGNIJ = SIGNJ*(-ODDIEL)
            do IIORB=IORB1+1,IORB2-1
              IJ = IJOFF+IIORB-IORBMIN+1
              if ((IJ <= 0) .or. (IJ > NIJ)) then
                write(u6,*) ' PROBLEMO ADADS1 : IJ : ',IJ
                write(u6,*) ' IJOFF IORBMIN ',IJOFF,IORBMIN
                write(u6,*) ' IIORB JJORB ',IIORB,JJORB
                !stop
                call SYSABENDMSG('lucia_util/adads1_gas','Internal error','')
              end if
              ILEX = ILEX2+IZ(IIORB,IEL+2)
              IACT = IREO(ILEX)
#             ifdef _DEBUGPRINT_
              write(u6,*) ' IIORB JJORB ',IIORB,JJORB
              write(u6,*) ' IJ ILEX,IACT',IJ,ILEX,IACT
              write(u6,*) ' ILEX0 ILEX1 ILEX2 ILEX ',ILEX0,ILEX1,ILEX2,ILEX
#             endif
              I1(KKSTR-KMIN+1,IJ) = IACT
              XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
              if (IJ < 0) then
                write(u6,*) ' NEGATIVE IJ in ADADS1'
                !stop ' NEGATIVE IJ in ADADS1'
                call SYSABENDMSG('lucia_util/adads1_gas','Internal error','')
              end if
            end do
          end if
        end do
      end do
    end if
  end do
end do

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

end subroutine ADADS1_GAS
