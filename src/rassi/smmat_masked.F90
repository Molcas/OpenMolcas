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

subroutine SMMAT_MASKED(PROP,PRMAT,NSS,ISONUM,ISPINCMP,ISS_INDEX,IST,INUM,JST,JNUM)

use Cntrl, only: NSTATE, NPROP, ICOMP, ISOCMP, PNAME, PTYPE, SOPRNM, SOPRTP
use Constants, only: Zero, One, Half
use Definitions, only: wp, u6

implicit none
integer NSS, ISONUM, ISPINCMP
real*8 PRMAT(NSS,NSS)
integer INUM, JNUM
integer ISS_INDEX(NSTATE+1), IST(INUM), JST(JNUM)
real*8 PROP(NSTATE,NSTATE,NPROP)
real*8, external :: DCLEBS
integer IPRNUM, IPRCMP, IPROP, I, ISTATE, ISS, MPLET1, MSPROJ1, IFSPIN, J, JSTATE, JSS, MPLET2, MSPROJ2
real*8 S1, S2, SM1, SM2, SXMER, SYMEI, SZMER, SMINUS, SPLUS, FACT, CGM, CG0, CGP, CGX, CGY, EXPKR

IPRNUM = -1
IPRCMP = 0
PRMAT(:,:) = Zero
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
      if (PNAME(IPRNUM)(1:5) == 'TMOM0') then
        IFSPIN = 2
        IPRCMP = ISPINCMP
      end if
      exit
    end if
  end do
end if
if (IPRNUM == -1) then
  write(u6,*) 'SMMAT_MASKED, Abend IPRNUM == -1'
  write(u6,*) 'SMMAT_MASKED, PRLBL=','>',PNAME(ISONUM),'<'
  call Abend()
end if

! Mapping from spin states to spin-free state and to spin:
do I=1,INUM
  ISTATE = IST(I)
  ISS = ISS_INDEX(ISTATE)
  MPLET1 = ISS_INDEX(ISTATE+1)-ISS_INDEX(ISTATE)
  S1 = Half*real(MPLET1-1,kind=wp)
  do MSPROJ1=-MPLET1+1,MPLET1-1,2
    SM1 = Half*real(MSPROJ1,kind=wp)
    ISS = ISS+1

    do J=1,JNUM
      JSTATE = JST(J)
      JSS = ISS_INDEX(JSTATE)
      MPLET2 = ISS_INDEX(JSTATE+1)-ISS_INDEX(JSTATE)
      S2 = Half*real(MPLET2-1,kind=wp)
      do MSPROJ2=-MPLET2+1,MPLET2-1,2
        SM2 = Half*real(MSPROJ2,kind=wp)
        JSS = JSS+1

        if ((IFSPIN == 0) .and. (IPRNUM /= 0)) then
          if ((MPLET1 == MPLET2) .and. (MSPROJ1 == MSPROJ2)) then
            PRMAT(ISS,JSS) = PROP(ISTATE,JSTATE,IPRNUM)
          else
            PRMAT(ISS,JSS) = Zero
          end if
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
          else
            PRMAT(ISS,JSS) = Zero
          end if
        else if (IFSPIN == 2) then

          ! The code here is a replica from smmat.
          ! Look in that source for comments.

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

return

end subroutine SMMAT_MASKED
