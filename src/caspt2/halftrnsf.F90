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
! Copyright (C) Francesco Aquilante                                    *
!               2007, Per Ake Malmqvist                                *
!***********************************************************************

subroutine HALFTRNSF(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,JREDC,CMO,NCMO,ISTART,NUSE,ipChoT,Work,nWork)
!********************************************************
!   Author: F. Aquilante as subroutine cho_vtra
!   Modified PAM 2007: Use ordinary CMO array without restructuring
!   This amounts to (1) having symmetry blocks of MO:s accessed
!   untransposed, (2) not through pointers to workspace array Work()
!   but giving the index of start orbital in each symmetry instead,
!   and (3) having CMO array in call parameter list.
!
!   Purpose:  SCR(lscr) contains JNUM cholesky vectors
!             starting from JVEC1 and stored in reduced
!             sets. The routine performs an MOs half transformation
!             of these elements in a set of target
!             arrays identified by the pointers ipChoT.
!             In the target arrays, the vectors are
!             stored in full dimension and as a
!             subset of a a given NUMV number of vectors.
!             Each pointer should thereby point to a
!             location where the corresponding Cholesky
!             vector of a given unique symmetry pair
!             of indices has to be stored
!
!   Input:
!       jVref =  index of the first vector to be transformed
!                computed wrt to the first vector in the
!                target arrays (see calling routine)
!
!       Jvec1 =  first vector to be transformed
!       JNum  =  # of vectors to be transformed
!
!       NumV  =  total # of vectors in the target arrays
!
!       JREDC :  reduced set in core at the moment of
!                the first call to the routine.
!                Can be set to -1 by the calling routine
!
!********************************************************

use Symmetry_Info, only: Mul
use Cholesky, only: iBas, iiBstR, IndRed, InfVec, iRS2F, nBas, nDimRS, nnBstR, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: irc, JREDC
integer(kind=iwp), intent(in) :: lscr, jVref, JVEC1, JNUM, NUMV, JSYM, NCMO, ISTART(8), NUSE(8), ipChoT(8), nWork
real(kind=wp), intent(inout) :: Scr(lscr)
real(kind=wp), intent(in) :: CMO(NCMO)
real(kind=wp), intent(out) :: Work(nWork)
integer(kind=iwp) :: iag, ias, ibg, ibs, IOC, IOFFC(8), iRab, ISCA, ISCB, ISYM, iSyma, iSymb, iSymp, jRab, JRED, JVEC, JVTRNS, &
                     kchot, kRab, kscr, NBA, NBB, nElem, NREAD, NUSEA, NUSEB
integer(kind=iwp), external :: cho_isao
integer(kind=iwp), parameter :: iLoc = 3 ! this means 'use scratch location in reduced index arrays'

! Offset counter into CMO array:
IOC = 0
do ISYM=1,NSYM
  IOFFC(ISYM) = IOC
  IOC = IOC+NBAS(ISYM)**2
end do

do iSymp=1,nSym
  if (nUse(iSymp) /= 0) then
    iSymb = Mul(JSYM,iSymp)
    nElem = nUse(iSymp)*nBas(iSymb)*NUMV
    call dCopy_(nElem,[Zero],0,Work(ipChoT(iSymp)),1)
  end if
end do

NREAD = 0
do JVEC=1,JNUM
  ! JVTRNS=Cholesky vector to be transformed.
  JVTRNS = JVEC1-1+JVEC
  ! JRED: Which reduced set does it belong to:
  JRED = InfVec(JVTRNS,2,JSYM)

  ! Is it the same still?
  if (JRED /= JREDC) then
    ! It is not. Tables must be regenerated with information that is
    ! common to this reduced set.
    write(u6,*) ' Rats! It was assumed that the Cholesky vectors'
    write(u6,*) ' in HALFTRNSF all belonged to a given reduced'
    write(u6,*) ' set, but they don''t!'
    write(u6,*) ' JRED, JREDC:',JRED,JREDC
    write(u6,*) ' Back to the drawing board?'
    write(u6,*) ' Let the program continue and see what happens.'
    call Cho_X_SetRed(irc,iLoc,JRED)
    JREDC = JRED
  end if

  kscr = NREAD
  NREAD = NREAD+nDimRS(JSYM,JRED)

  if (JSYM == 1) then
    ! L(a,b,J)=L(b,a,J); only a >= b stored

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,3)

      ! Global address:
      iag = iRS2F(1,iRab)
      ibg = iRS2F(2,iRab)

      iSyma = cho_isao(iag)

      kscr = kscr+1

      ISYMB = ISYMA

      NUSEA = nUse(iSyma)
      NUSEB = nUse(iSymb)

      if (NUSEA /= 0) then

        ias = iag-ibas(iSyma)
        ibs = ibg-ibas(iSyma)
        !  L(p,J,b) = sum( C(a,p)* L(a,b,J), a=1..NBA), where p=1..NUSEA
        !  L(p,J,a) = sum( C(b,p)* L(b,a,J), b=1..NBB), where p=1..NUSEB
        !  ----------------------------------------

        NBA = NBAS(ISYMA)
        ISCA = IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA)-1)
        NBB = NBAS(ISYMB)
        ISCB = IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB)-1)

        kchot = ipChoT(iSymb)+NUSEB*(jVref+JVEC-2+NUMV*(ias-1))

        call DAXPY_(NUSEA,Scr(kscr),CMO(ISCB),NBB,Work(kchot),1)

        if (IAS /= IBS) then
          kchot = ipChoT(iSyma)+NUSEA*(jVref+JVEC-2+NUMV*(ibs-1))

          call DAXPY_(NUSEA,Scr(kscr),CMO(ISCA),NBA,Work(kchot),1)
        end if

        ! End of NUSE  test
      end if

      ! End of loop over basis function pair index JRAB
    end do

  else
    ! jSym /= 1

    do jRab=1,nnBstR(jSym,iLoc)

      kRab = iiBstr(jSym,iLoc)+jRab
      iRab = IndRed(kRab,3)

      ! Global address:
      iag = iRS2F(1,iRab)
      ibg = iRS2F(2,iRab)

      ! iSyma = cho_isao(iag) = symmetry block of basis function iag
      iSyma = cho_isao(iag)
      ! iSyma > isymb since jsym /= 1 and a >= b
      iSymb = Mul(jSym,iSyma)
      NUSEA = nUse(iSyma)
      NUSEB = nUse(iSymb)

      kscr = kscr+1

      ias = iag-ibas(iSyma)
      ibs = ibg-ibas(iSymb)

      if (NUSEA /= 0) then
        !  L(p,J,b) = sum( C(a,p)*L(a,b,J), a=1..NBA), where p=1..NUSEA

        NBA = NBAS(ISYMA)
        ISCA = IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA)-1)

        kchot = ipChoT(iSyma)+nUse(iSyma)*(jVref+JVEC-2+NUMV*(ibs-1))

        call DAXPY_(NUSEA,Scr(kscr),CMO(ISCA),NBA,Work(kchot),1)

        ! End of NUSE  test
      end if

      if (NUSEB /= 0) then
        ! L(p,J,a) = sum( C(b,p)*L(b,a,J), b=1..NBA), where p=1..NUSEB

        NBB = NBAS(ISYMB)
        ISCB = IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB)-1)

        kchot = ipChoT(iSymb)+NUSEB*(jVref+JVEC-2+NUMV*(ias-1))

        call DAXPY_(NUSEB,Scr(kscr),CMO(ISCB),NBB,Work(kchot),1)

        ! End of NUSE  test
      end if

      ! End of loop over basis function pair index JRAB
    end do

    ! End of JSYM /= 1 test
  end if

  ! End of JSYM loop
end do

irc = 0

end subroutine HALFTRNSF
