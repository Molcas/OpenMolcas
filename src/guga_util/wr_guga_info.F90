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

subroutine WR_GUGA(Lu,iOpt,iDisk,NFREF,S,N,LN,NSYM,IR1,IR2,IFIRST,INTNUM,LSYM,NREF,LN1,NRLN1,NSH,NISH,MxSym,JRC,nJRC,JJS,nJJS, &
                   NVAL,IOCR,nIOCR)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, MxSym, nJRC, nJJS, nIOCR
integer(kind=iwp), intent(inout) :: iDisk, NFREF, N, LN, NSYM, IR1, IR2, IFIRST, INTNUM, LSYM, NREF, LN1, NRLN1, NSH(MxSym), &
                                    NISH(MxSym), JRC(nJRC), JJS(nJJS), NVAL(MxSym), IOCR(nIOCR)
real(kind=wp), intent(inout) :: S
integer(kind=iwp) :: NFREF_(1), N_(1), LN_(1), NSYM_(1), IR1_(1), IR2_(1), IFIRST_(1), INTNUM_(1), LSYM_(1), NREF_(1), LN1_(1), &
                     NRLN1_(1)
real(kind=wp) :: S_(1)

if (iOpt == 1) then
  NFREF_(1) = NFREF
  S_(1) = S
  N_(1) = N
  LN_(1) = LN
  NSYM_(1) = NSYM
  IR1_(1) = IR1
  IR2_(1) = IR2
  IFIRST_(1) = IFIRST
  INTNUM_(1) = INTNUM
  LSYM_(1) = LSYM
  NREF_(1) = NREF
  LN1_(1) = LN1
  NRLN1_(1) = NRLN1
end if
call iDaFile(Lu,iOpt,NFREF_,1,iDisk)
call dDaFile(Lu,iOpt,S_,1,iDisk)
call iDaFile(Lu,iOpt,N_,1,iDisk)
call iDaFile(Lu,iOpt,LN_,1,iDisk)
call iDaFile(Lu,iOpt,NSYM_,1,iDisk)
call iDaFile(Lu,iOpt,IR1_,1,iDisk)
call iDaFile(Lu,iOpt,IR2_,1,iDisk)
call iDaFile(Lu,iOpt,IFIRST_,1,iDisk)
call iDaFile(Lu,iOpt,INTNUM_,1,iDisk)
call iDaFile(Lu,iOpt,LSYM_,1,iDisk)
call iDaFile(Lu,iOpt,NREF_,1,iDisk)
call iDaFile(Lu,iOpt,LN1_,1,iDisk)
call iDaFile(Lu,iOpt,NRLN1_,1,iDisk)
if (iOpt == 2) then
  NFREF = NFREF_(1)
  S = S_(1)
  N = N_(1)
  LN = LN_(1)
  NSYM = NSYM_(1)
  IR1 = IR1_(1)
  IR2 = IR2_(1)
  IFIRST = IFIRST_(1)
  INTNUM = INTNUM_(1)
  LSYM = LSYM_(1)
  NREF = NREF_(1)
  LN1 = LN1_(1)
  NRLN1 = NRLN1_(1)
end if
call iDaFile(Lu,iOpt,NSH,MxSym,iDisk)
call iDaFile(Lu,iOpt,NISH,MxSym,iDisk)
call iDaFile(Lu,iOpt,JRC,nJRC,iDisk)
call iDaFile(Lu,iOpt,JJS,nJJS,iDisk)
call iDaFile(Lu,iOpt,NVAL,MxSym,iDisk)
call iDaFile(Lu,iOpt,IOCR,nIOCR,iDisk)

return

end subroutine WR_GUGA
