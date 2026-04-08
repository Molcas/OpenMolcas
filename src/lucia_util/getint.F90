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
subroutine GETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
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

use Index_Functions, only: iTri
use lucia_data, only: IBSO, NOBPTS, NTOOBS
use wadr, only: TUVX
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: XINT(*)
integer(kind=iwp), intent(in) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM, ICOUL
integer(kind=iwp) :: I, IINT, IJKL, ILKJ, IMIN, IOFF, IORB, J, JMIN, JOFF, JORB, K, KOFF, KORB, L, LOFF, LORB
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NI, NIK, NJ, NJL, NK, NL

!write(u6,*) ' I_USE_SIMTRH in GETINT =',I_USE_SIMTRH
write(u6,*) ' GETINT : ICOUL = ',ICOUL
write(u6,*) 'ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM :'
write(u6,'(8I4)') ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM
#endif
! Read integrals in in RASSCF format

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
IOFF = IBSO(ISM)+sum(NOBPTS(1:ITP-1,ISM))

JOFF = IBSO(JSM)+sum(NOBPTS(1:JTP-1,JSM))

KOFF = IBSO(KSM)+sum(NOBPTS(1:KTP-1,KSM))

LOFF = IBSO(LSM)+sum(NOBPTS(1:LTP-1,LSM))

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
        IJKL = iTri(iTri(I,J),iTri(K,L))
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
          ILKJ = iTri(iTri(I,L),iTri(K,J))
          iInt = iInt+1
          XInt(iInt) = XInt(iInt)-TUVX(ilkj)
        end do
      end do
    end do
  end do
end if

#ifdef _DEBUGPRINT_
if (ITP == 0) then
  NI = NTOOBS(ISM)
else
  NI = NOBPTS(ITP,ISM)
end if

if (KTP == 0) then
  NK = NTOOBS(KSM)
else
  NK = NOBPTS(KTP,KSM)
end if

if (IKSM == 0) then
  NIK = NI*NK
else
  NIK = nTri_Elem(NI)
end if

if (JTP == 0) then
  NJ = NTOOBS(JSM)
else
  NJ = NOBPTS(JTP,JSM)
end if

if (LTP == 0) then
  NL = NTOOBS(LSM)
else
  NL = NOBPTS(LTP,LSM)
end if

if (JLSM == 0) then
  NJL = NJ*NL
else
  NJL = nTri_Elem(NJ)
end if
write(u6,*) ' 2 electron integral block for TS blocks'
write(u6,*) ' Ixchng :',IXCHNG
write(u6,*) ' After GETINC'
write(u6,'(1X,4(A,I2,A,I2,A))') '(',ITP,',',ISM,')','(',JTP,',',JSM,')','(',KTP,',',KSM,')','(',LTP,',',LSM,')'
call WRTMAT(XINT,NIK,NJL,NIK,NJL)
#endif

end subroutine GETINT
