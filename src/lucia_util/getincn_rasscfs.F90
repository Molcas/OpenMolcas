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

subroutine GETINCN_RASSCFS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
! Obtain integrals
!
!     ICOUL = 0 :
!                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
!
!     ICOUL = 1 :
!                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
!
!     ICOUL = 2 :  XINT(IL,JK) = (IJ!KL)         for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
!
! Storing for ICOUL = 1 not working if IKSM or JLSM /= 0
!
! Version for integrals stored in TUVX
!
! If type equals zero, all integrals of given type are fetched
! (added aug 8, 98)

use lucia_data, only: IBSO, NOBPTS, NTOOBS
use wadr, only: TUVX
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: XINT(*)
integer(kind=iwp) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM, ICOUL
integer(kind=iwp) :: I, IINT, IITP, IJ, IJKL, IL, ILKJ, IMIN, IOFF, IORB, J, JJTP, JMIN, JOFF, JORB, K, KJ, KKTP, KL, KOFF, KORB, &
                     L, LLTP, LOFF, LORB

if (ITP >= 1) then
  iOrb = NOBPTS(ITP,ISM)
else
  IORB = NTOOBS(ISM)
end if

if (JTP >= 1) then
  jOrb = NOBPTS(JTP,JSM)
else
  JORB = NTOOBS(JSM)
end if

if (KTP >= 1) then
  kOrb = NOBPTS(KTP,KSM)
else
  KORB = NTOOBS(KSM)
end if

if (LTP >= 1) then
  lOrb = NOBPTS(LTP,LSM)
else
  LORB = NTOOBS(LSM)
end if
! Offsets relative to start of all orbitals, symmetry ordered
IOFF = IBSO(ISM)
do IITP=1,ITP-1
  IOFF = IOFF+NOBPTS(IITP,ISM)
end do

JOFF = IBSO(JSM)
do JJTP=1,JTP-1
  JOFF = JOFF+NOBPTS(JJTP,JSM)
end do

KOFF = IBSO(KSM)
do KKTP=1,KTP-1
  KOFF = KOFF+NOBPTS(KKTP,KSM)
end do

LOFF = IBSO(LSM)
do LLTP=1,LTP-1
  LOFF = LOFF+NOBPTS(LLTP,LSM)
end do

! Collect Coulomb terms

iInt = 0
do l=lOff,lOff+lOrb-1
  jMin = jOff
  if (JLSM /= 0) jMin = l
  do j=jMin,jOff+jOrb-1
    do k=kOff,kOff+kOrb-1
      iMin = iOff
      if (IKSM /= 0) iMin = k
      if (ICOUL == 1) then
        ! Address before integral (1,j!k,l)
        IINT = (L-LOFF)*Jorb*Korb*Iorb+(K-KOFF)*Jorb*Iorb+(J-JOFF)*Iorb
      else if (ICOUL == 2) then
        ! Address before (1L,JK)
        IINT = (K-KOFF)*JORB*LORB*IORB+(J-JOFF)*LORB*IORB+(L-LOFF)*IORB
      end if

      do i=iMin,iOff+iOrb-1
        IJ = max(I,J)*(max(I,J)-1)/2+min(I,J)
        KL = max(K,L)*(max(K,L)-1)/2+min(K,L)
        IJKL = max(IJ,KL)*(max(IJ,KL)-1)/2+min(IJ,KL)
        ! Next line inserted by Jesper: "I don't think iInt should be the same
        ! for all i"
        iInt = iInt+1
        Xint(iInt) = TUVX(ijkl)
      end do
    end do
  end do
end do

! Collect Exchange terms

if (IXCHNG /= 0) then

  iInt = 0
  do l=lOff,lOff+lOrb-1
    jMin = jOff
    if (JLSM /= 0) jMin = l
    do j=jMin,jOff+jOrb-1
      do k=kOff,kOff+kOrb-1
        iMin = iOff
        if (IKSM /= 0) iMin = k
        if (ICOUL == 1) then
          ! Address before integral (1,j!k,l)
          IINT = (L-LOFF)*Jorb*Korb*Iorb+(K-KOFF)*Jorb*Iorb+(J-JOFF)*Iorb
        else if (ICOUL == 2) then
          ! Address before (1L,JK)
          IINT = (K-KOFF)*JORB*LORB*IORB+(J-JOFF)*LORB*IORB+(L-LOFF)*IORB
        end if

        do i=iMin,iOff+iOrb-1
          IL = max(I,L)*(max(I,L)-1)/2+min(I,L)
          KJ = max(K,J)*(max(K,J)-1)/2+min(K,J)
          ILKJ = max(IL,KJ)*(max(IL,KJ)-1)/2+min(IL,KJ)
          iInt = iInt+1
          XInt(iInt) = XInt(iInt)-TUVX(ilkj)
        end do
      end do
    end do
  end do
end if

end subroutine GETINCN_RASSCFS
