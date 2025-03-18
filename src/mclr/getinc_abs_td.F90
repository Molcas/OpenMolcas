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

subroutine GETINC_ABS_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IKSM,JLSM,INTLST,IJKLOF,NSMOB,ICTL)
! Obtain integrals
! ICOUL = 0 :      XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
! ICOUL = 1 :      XINT(IJ,KL) = (IJ!KL)
!
! Version for integrals stored in INTLST

use MCLR_Data, only: NACOB, IBTSOB, NTSOB

implicit none
real*8 XINT(*)
integer ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IKSM, JLSM
real*8 Intlst(*)
integer NSMOB
integer IJKLof(NsmOB,NsmOb,NsmOB)
integer ICTL
! Local variables
integer iOrb, jOrb, kOrb, lOrb
integer iOff, jOff, kOff, lOff
integer iBas, jBas, kBas, lBas
integer iInt
integer NTASH, JMIN, IMIN, IJ, KL, IJKL, IL, JK, ILJK
! Statement function
integer i, j, itri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

iOrb = NTSOB(ITP,ISM)
jOrb = NTSOB(JTP,JSM)
kOrb = NTSOB(KTP,KSM)
lOrb = NTSOB(LTP,LSM)
iOff = IBTSOB(ITP,ISM)
jOff = IBTSOB(JTP,JSM)
kOff = IBTSOB(KTP,KSM)
lOff = IBTSOB(LTP,LSM)
ntash = nacob

! Collect Coulomb terms

if (ICTL == 1) then
  iint = 1
  do lBas=lOff,lOff+lOrb-1
    jMin = jOff
    if (JLSM /= 0) jMin = lBas
    do jBas=jMin,jOff+jOrb-1
      do kBas=kOff,kOff+kOrb-1
        iMin = iOff
        if (IKSM /= 0) iMin = kBas
        do iBas=iMin,iOff+iOrb-1
          ij = jbas+ntash*(ibas-1)
          kl = lbas+ntash*(kbas-1)
          ijkl = itri(ij,kl)
          Xint(iInt) = Intlst(ijkl)
          iInt = iInt+1
        end do
      end do
    end do
  end do
else if (ICTL == 4) then

  ! Collect Coulomb-Exchange terms

  iint = 1
  do lBas=lOff,lOff+lOrb-1
    jMin = jOff
    if (JLSM /= 0) jMin = lBas
    do jBas=jMin,jOff+jOrb-1
      do kBas=kOff,kOff+kOrb-1
        iMin = iOff
        if (IKSM /= 0) iMin = kBas
        do iBas=iMin,iOff+iOrb-1
          il = ibas+ntash*(lbas-1)
          jk = kbas+ntash*(jbas-1)
          iljk = itri(il,jk)
          ij = ibas+ntash*(jbas-1)
          kl = kbas+ntash*(lbas-1)
          ijkl = itri(ij,kl)
          XInt(iInt) = Intlst(ijkl)-Intlst(iljk)
          iInt = iInt+1
        end do
      end do
    end do
  end do
else
  call Abend()
end if

! Avoid unused argument warnings
if (.false.) call Unused_integer_array(IJKLOF)

end subroutine GETINC_ABS_td
