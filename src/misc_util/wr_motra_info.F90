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

subroutine WR_MOTRA_Info(Lu,iOpt,iDisk,TCONEMO,nTCONEMO,ECOR,NSYM,NBAS,NORB,NFRO,NDEL,MxSym,BSLBL,nBSLBL)

implicit real*8(a-h,o-z)
#include "SysDef.fh"
integer TCONEMO(nTCONEMO), nBas(MxSym), nOrb(MxSym), nFro(MxSym), nDel(MxSym)
character BSLBL(nBSLBL)*1

call iDafile(Lu,iOpt,TCONEMO,nTCONEMO,iDisk)
call s_dDafile_motra(Lu,iOpt,ECor,1,iDisk)
call s_iDafile_motra(Lu,iOpt,nSym,1,iDisk)
call iDafile(Lu,iOpt,nBas,MxSym,iDisk)
call iDafile(Lu,iOpt,nOrb,MxSym,iDisk)
call iDafile(Lu,iOpt,nFro,MxSym,iDisk)
call iDafile(Lu,iOpt,nDel,MxSym,iDisk)
call cDafile(Lu,iOpt,BSLBL,nBSLBL,iDisk)

return

! This is to allow type punning without an explicit interface
contains

subroutine s_iDaFile_motra(Lu,iOpt,Buf,lBuf_,iDisk_)

  use iso_c_binding

  integer Lu, iOpt, lBuf_, iDisk_
  integer, target :: Buf
  integer, pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf),pBuf,[1])
  call iDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine s_iDaFile_motra

subroutine s_dDaFile_motra(Lu,iOpt,Buf,lBuf_,iDisk_)

  use iso_c_binding

  integer Lu, iOpt, lBuf_, iDisk_
  real*8, target :: Buf
  real*8, pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf),pBuf,[1])
  call dDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine s_dDaFile_motra

end subroutine WR_MOTRA_Info
