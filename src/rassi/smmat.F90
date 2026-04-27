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

subroutine SMMAT(PROP,PRMAT,NSS,ISONUM,ISPINCMP)

use rassi_global_arrays, only: JBNUM
use Cntrl, only: ICOMP, ISOCMP, MLTPLT, NPROP, NSTATE, PNAME, PTYPE, SOPRNM, SOPRTP
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP)
integer(kind=iwp), intent(in) :: NSS, ISONUM, ISPINCMP
real(kind=wp), intent(inout) :: PRMAT(NSS,NSS)
integer(kind=iwp) :: IFSPIN, IPRCMP, IPRNUM, IPROP, ISS, ISTATE, JOB1, JOB2, JSS, JSTATE, MPLET1, MPLET2, MSPROJ1, MSPROJ2
real(kind=wp) :: CG0, CGM, CGP, CGX, CGY, EXPKR, FACT, s1, s2, SM1, SM2, SMINUS, SPLUS, SXMER, SYMEI, SZMER
real(kind=wp), external :: DCLEBS

IPRNUM = -1
IPRCMP = 0
! IFSPIN takes values the values 0,1,2
! 0 = spin free property
! 1 = spin operator (S)
! 2 = spin dependent property, triplet operator
IFSPIN = 0

if (ISONUM == 0) then
  IPRNUM = 0
  IFSPIN = 1
  IPRCMP = ISPINCMP
else
  do IPROP=1,NPROP
    if ((PNAME(IPROP) == SOPRNM(ISONUM)) .and. (PTYPE(IPROP) == SOPRTP(ISONUM)) .and. (ICOMP(IPROP) == ISOCMP(ISONUM))) then
      IPRNUM = IPROP
      if (PNAME(IPRNUM)(1:3) == 'PSO') then
        IFSPIN = 0
        IPRCMP = ISPINCMP
      end if
      if (PNAME(IPRNUM)(1:5) == 'TMOM0') then
        IFSPIN = 2
        IPRCMP = ISPINCMP
      end if
      if ((PNAME(IPRNUM) == 'MLTPL  0') .and. (PTYPE(IPRNUM) == 'ANTITRIP')) then
        IFSPIN = 2
        IPRCMP = ISPINCMP
      end if
      if ((PNAME(IPRNUM) == 'MLTPL  1') .and. (PTYPE(IPRNUM) == 'ANTITRIP')) then
        IFSPIN = 2
        IPRCMP = ISPINCMP
      end if
      exit
    end if
  end do
end if
if (IPRNUM == -1) then
  write(u6,*) 'SMMAT, Abend IPRNUM == -1'
  write(u6,*) 'SMMAT, PRLBL=','>',PNAME(ISONUM),'<'
  call Abend()
end if

! Mapping from spin states to spin-free state and to spin:
ISS = 0
do ISTATE=1,NSTATE
  JOB1 = JBNUM(ISTATE)
  MPLET1 = MLTPLT(JOB1)
  S1 = Half*real(MPLET1-1,kind=wp)

  do MSPROJ1=-MPLET1+1,MPLET1-1,2
    SM1 = Half*real(MSPROJ1,kind=wp)
    ISS = ISS+1

    JSS = 0

    do JSTATE=1,NSTATE
      JOB2 = JBNUM(JSTATE)
      MPLET2 = MLTPLT(JOB2)
      S2 = Half*real(MPLET2-1,kind=wp)

      do MSPROJ2=-MPLET2+1,MPLET2-1,2
        SM2 = Half*real(MSPROJ2,kind=wp)
        JSS = JSS+1

        if ((IFSPIN == 0) .and. (IPRNUM /= 0)) then
          if ((MPLET1 == MPLET2) .and. (MSPROJ1 == MSPROJ2)) PRMAT(ISS,JSS) = PROP(ISTATE,JSTATE,IPRNUM)
        else if ((IFSPIN == 1) .and. (IPRNUM == 0)) then
          SXMER = Zero
          SYMEI = Zero
          SZMER = Zero
          SMINUS = Zero
          SPLUS = Zero
          if ((ISTATE == JSTATE) .and. (MPLET1 == MPLET2)) then
            if (MSPROJ1 == MSPROJ2-2) then
              SMINUS = sqrt((S1+SM2)*(S1-SM1))
              SXMER = Half*SMINUS
              SYMEI = Half*SMINUS
            else if (MSPROJ1 == MSPROJ2) then
              SZMER = SM1
            else if (MSPROJ1 == MSPROJ2+2) then
              SPLUS = sqrt((S1+SM1)*(S1-SM2))
              SXMER = Half*SPLUS
              SYMEI = -Half*SPLUS
            end if
            if (IPRCMP == 1) then
              PRMAT(ISS,JSS) = SXMER
            else if (IPRCMP == 2) then
              PRMAT(ISS,JSS) = SYMEI
            else if (IPRCMP == 3) then
              PRMAT(ISS,JSS) = SZMER
            end if
          end if
        else if (IFSPIN == 2) then
          ! 1-electron triplet operator so only Delta S =0,+-1 and Delta MS =0,+-1
          ! Notice S1,S2 and SM1,SM2 need not to be integers
          ! Hence MPLET1,MPLET2 and MSPROJ1,MSPROJ2 are used
          ! Notice SMINUS and SPLUS is interchanged compared to above
          ! Notice that the Y part is imaginary

          ! see section 3 (Spin_orbit coupling in RASSI) in
          ! P A Malmqvist, et. al CPL, 357 (2002) 230-240
          ! for details

          ! Note that we work on the x, y, and z components at this time.

          ! What follows applies only to the exact operator for the
          ! transition moment.

          ! On page 234 we have the notation V^{AB}(x), that is the
          ! potential has Cartesian components. Here, however, this is
          ! partitioned in a slightly different way since we have that
          ! V^{AB}(x)=(k x e_l)_x V^{AB}. We will only handle the
          ! V^{AB} part.

          ! Hence, we will compute the contributions to T(i), i=x,y,z
          ! here and form the inner product
          ! (k x e_l)_i V^{AB} . T(i)
          ! outside the code.

          ! Note that this code strictly follows the code of soeig where
          ! the term L.S is added to the Hamiltonian.

          FACT = One/sqrt(real(MPLET1,kind=wp))
          if (MPLET1 == MPLET2-2) FACT = -FACT
          CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
          CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
          CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
          CGX = sqrt(Half)*(CGM-CGP)
          CGY = -sqrt(Half)*(CGM+CGP)

          EXPKR = PROP(ISTATE,JSTATE,IPRNUM)

          if (IPRCMP == 1) then
            EXPKR = EXPKR*CGX
          else if (IPRCMP == 2) then
            EXPKR = EXPKR*CGY
          else if (IPRCMP == 3) then
            EXPKR = EXPKR*CG0
          end if
          PRMAT(ISS,JSS) = EXPKR
        end if

      end do !MSPROJ2
    end do !JSTATE
  end do !MSPROJ1
end do !ISTATE

end subroutine SMMAT
