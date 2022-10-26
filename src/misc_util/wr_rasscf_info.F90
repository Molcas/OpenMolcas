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
subroutine WR_RASSCF_Info(Lu,iOpt,iDisk,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,mxSym,NAME,nName,NCONF,HEADER,nHeader, &
                          TITLE,nTitle,POTNUC,LROOTS,NROOTS,IROOT,mxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)

implicit real*8(a-h,o-z)
#include "SysDef.fh"
integer nFro(MxSym), nISh(MxSym), nASh(MxSym), nDel(MxSym), nBas(MxSym), iRoot(MxRoot), nRS1(MxSym), nRS2(MxSym), nRS3(MxSym)
character Name(nName)*1, Header(nHeader)*1, Title(nTitle)*1
real*8 Weight(MxRoot)

call s_iDaFile_rasscf(Lu,iOpt,nActEl,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,iSpin,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,nSym,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,lSym,1,iDisk)
call iDaFile(Lu,iOpt,nFro,MxSym,iDisk)
call iDaFile(Lu,iOpt,nISh,MxSym,iDisk)
call iDaFile(Lu,iOpt,nASh,MxSym,iDisk)
call iDaFile(Lu,iOpt,nDel,MxSym,iDisk)
call iDaFile(Lu,iOpt,nBas,MxSym,iDisk)
call cDaFile(Lu,iOpt,Name,nName,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,nConf,1,iDisk)
call cDaFile(Lu,iOpt,Header,nHeader,iDisk)
call cDaFile(Lu,iOpt,Title,nTitle,iDisk)
call s_dDaFile_rasscf(Lu,iOpt,PotNuc,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,lRoots,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,nRoots,1,iDisk)
call iDaFile(Lu,iOpt,iRoot,MxRoot,iDisk)
call iDaFile(Lu,iOpt,nRS1,MxSym,iDisk)
call iDaFile(Lu,iOpt,nRS2,MxSym,iDisk)
call iDaFile(Lu,iOpt,nRS3,MxSym,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,nHole1,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,nElec3,1,iDisk)
call s_iDaFile_rasscf(Lu,iOpt,iPT2,1,iDisk)
call dDaFile(Lu,iOpt,Weight,MxRoot,iDisk)

return

! This is to allow type punning without an explicit interface
contains

subroutine s_iDaFile_rasscf(Lu,iOpt,Buf,lBuf_,iDisk_)

  use iso_c_binding

  integer Lu, iOpt, lBuf_, iDisk_
  integer, target :: Buf
  integer, pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf),pBuf,[1])
  call iDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine s_iDaFile_rasscf

subroutine s_dDaFile_rasscf(Lu,iOpt,Buf,lBuf_,iDisk_)

  use iso_c_binding

  integer Lu, iOpt, lBuf_, iDisk_
  real*8, target :: Buf
  real*8, pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf),pBuf,[1])
  call dDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine s_dDaFile_rasscf

end subroutine WR_RASSCF_Info
