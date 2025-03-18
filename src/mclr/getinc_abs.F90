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

subroutine GETINC_ABS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,ICOUL)
! Obtain integrals
! ICOUL = 0 :      XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
!                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
! ICOUL = 1 :      XINT(IJ,KL) = (IJ!KL)
!
! Version for integrals stored in INTLST

use MCLR_Data, only: NACOB, IBTSOB, NTSOB

implicit none
real*8 XINT(*)
integer ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM
real*8 Intlst(*)
integer NSMOB
integer IJKLof(NsmOB,NsmOb,NsmOB)
integer ICOUL
! Local variables
integer iOrb, jOrb, kOrb, lOrb
integer iOff, jOff, kOff, lOff
integer iBas, jBas, kBas, lBas
integer iInt, jInt

iOrb = NTSOB(ITP,ISM)
jOrb = NTSOB(JTP,JSM)
kOrb = NTSOB(KTP,KSM)
lOrb = NTSOB(LTP,LSM)
iOff = IBTSOB(ITP,ISM)
jOff = IBTSOB(JTP,JSM)
kOff = IBTSOB(KTP,KSM)
lOff = IBTSOB(LTP,LSM)

! Collect Coulomb terms

if (ICOUL == 0) then
  iint = 1
  do lBas=lOff,lOff+lOrb-1
    do jBas=jOff,jOff+jOrb-1
      do kBas=kOff,kOff+kOrb-1
        do iBas=iOff,iOff+iOrb-1
          jINT = (lBas-1)*nACOB**3+(kBas-1)*nACOB**2+(jBas-1)*nACOB+iBas
          Xint(iInt) = Intlst(jint)
          iInt = iInt+1
        end do
      end do
    end do
  end do

  ! Collect Exchange terms

  if (IXCHNG /= 0) then
    iint = 1
    do lBas=lOff,lOff+lOrb-1
      do jBas=joff,jOff+jOrb-1
        do kBas=kOff,kOff+kOrb-1
          do iBas=ioff,iOff+iOrb-1
            jINT = (jBas-1)*nACOB**3+(kBas-1)*nACOB**2+(lBas-1)*nACOB+iBas
            XInt(iInt) = XInt(iInt)-Intlst(jint)
            iInt = iInt+1
          end do
        end do
      end do
    end do
  end if
else if (ICOUL /= 0) then
  iint = 0
  do lBas=lOff,lOff+lOrb-1
    do kBas=kOff,kOff+kOrb-1
      do jBas=joff,jOff+jOrb-1
        do iBas=ioff,iOff+iOrb-1
          JINT = (LBAS-1)*nACOB**3+(KBAS-1)*nACOB**2+(JBAS-1)*nACOB+IBAS
          Xint(iInt) = Intlst(jint)
          iInt = iint+1
        end do
      end do
    end do
  end do
end if

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(IKSM)
  call Unused_integer(JLSM)
  call Unused_integer_array(IJKLOF)
end if

end subroutine GETINC_ABS
