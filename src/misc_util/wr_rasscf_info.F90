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

subroutine WR_RASSCF_Info(Lu,iOpt,iDisk,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,MxSym,BNAME,nName,NCONF,HEADER,nHeader, &
                          TITLE,nTitle,POTNUC,LROOTS,NROOTS,IROOT,MxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, MxSym, nName, nHeader, nTitle, MxRoot
integer(kind=iwp), intent(inout) :: iDisk, NACTEL, ISPIN, NSYM, LSYM, nFro(MxSym), nISh(MxSym), nASh(MxSym), nDel(MxSym), &
                                    nBas(MxSym), NCONF, LROOTS, NROOTS, iRoot(MxRoot), NRS1(MxSym), NRS2(MxSym), NRS3(MxSym), &
                                    NHOLE1, NELEC3, IPT2
character, intent(inout) :: BName(nName), Header(nHeader), Title(nTitle)
real(kind=wp), intent(inout) :: POTNUC, Weight(MxRoot)
integer(kind=iwp) :: idum(1)
real(kind=wp) :: dum(1)

idum(1) = nActEl
call iDaFile(Lu,iOpt,idum,1,iDisk)
nActEl = idum(1)
idum(1) = iSpin
call iDaFile(Lu,iOpt,idum,1,iDisk)
iSpin = idum(1)
idum(1) = nSym
call iDaFile(Lu,iOpt,idum,1,iDisk)
nSym = idum(1)
idum(1) = lSym
call iDaFile(Lu,iOpt,idum,1,iDisk)
lSym = idum(1)
call iDaFile(Lu,iOpt,nFro,MxSym,iDisk)
call iDaFile(Lu,iOpt,nISh,MxSym,iDisk)
call iDaFile(Lu,iOpt,nASh,MxSym,iDisk)
call iDaFile(Lu,iOpt,nDel,MxSym,iDisk)
call iDaFile(Lu,iOpt,nBas,MxSym,iDisk)
call cDaFile(Lu,iOpt,BName,nName,iDisk)
idum(1) = nConf
call iDaFile(Lu,iOpt,idum,1,iDisk)
nConf = idum(1)
call cDaFile(Lu,iOpt,Header,nHeader,iDisk)
call cDaFile(Lu,iOpt,Title,nTitle,iDisk)
dum(1) = PotNuc
call dDaFile(Lu,iOpt,dum,1,iDisk)
PotNuc = dum(1)
idum(1) = lRoots
call iDaFile(Lu,iOpt,idum,1,iDisk)
lRoots = idum(1)
idum(1) = nRoots
call iDaFile(Lu,iOpt,idum,1,iDisk)
nRoots = idum(1)
call iDaFile(Lu,iOpt,iRoot,MxRoot,iDisk)
call iDaFile(Lu,iOpt,nRS1,MxSym,iDisk)
call iDaFile(Lu,iOpt,nRS2,MxSym,iDisk)
call iDaFile(Lu,iOpt,nRS3,MxSym,iDisk)
idum(1) = nHole1
call iDaFile(Lu,iOpt,idum,1,iDisk)
nHole1 = idum(1)
idum(1) = nElec3
call iDaFile(Lu,iOpt,idum,1,iDisk)
nElec3 = idum(1)
idum(1) = iPT2
call iDaFile(Lu,iOpt,idum,1,iDisk)
iPT2 = idum(1)
call dDaFile(Lu,iOpt,Weight,MxRoot,iDisk)

return

end subroutine WR_RASSCF_Info
