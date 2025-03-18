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

subroutine DMXMT(A,LDA,NTA,B,LDB,NTB,C,NRC,NSUM,NCC)

implicit none
integer IROW, ICOL, LDA, LDB, NRC, NSUM, NCC, IND, ISUM
character NTA, NTB
real*8 A(LDA,*), B(LDB,*), C(*), SUM

IND = 0
if ((NTA == 'N') .and. (NTB == 'N')) then
  do ICol=1,NRC
    do IRow=icol,nrc
      SUM = 0.0d0
      do ISUM=1,NSUM
        SUM = SUM+A(IROW,ISUM)*B(ISUM,iCOL)
      end do
      IND = IND+1
      C(IND) = SUM
    end do
  end do
!else if ((NTA == 'T') .and. (NTB == 'N')) then
!  do IROW=1,NRC
!    do ICOL=0,IROW
!      SUM = 0.0D0
!      do ISUM=0,NSUM-1
!        SUM = SUM+A(ISUM,IROW)*B(ISUM,iCOL)
!      end do
!      IND = IND+1
!      C(IND) = SUM
!    end do
!  end do
!else if ((NTA == 'N') .and. (NTB == 'T')) then
!  do IROW=1,NRC
!    do ICOL=0,IROW
!      SUM = 0.0D0
!      do ISUM=0,NSUM-1
!        SUM = SUM+A(IROW,ISUM)*B(ICOL,iSUM)
!      end do
!      IND = IND+1
!      C(IND) = SUM
!    end do
!  end do
!else if ((NTA == 'T') .and. (NTB == 'T')) then
!  do IROW=1,NRC
!    do ICOL=0,IROW
!      SUM = 0.0D0
!      do ISUM=0,NSUM-1
!        SUM = SUM+A(ISUM,IROW)*B(ICOL,iSUM)
!      end do
!      IND = IND+1
!      C(IND) = SUM
!    end do
!  end do
else
  call SysHalt('dmxmt')
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(NCC)

end subroutine DMXMT
