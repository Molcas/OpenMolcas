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

subroutine READCI(ISTATE,SGS,CIS,NCI,CI)

use rassi_aux, only: ipglob
use rassi_global_arrays, only: JBNUM, LROOT
use gugx, only: SGStruct, CIStruct
#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_exists_attr, mh5_fetch_attr, mh5_fetch_dset, mh5_close_file
use Cntrl, only: NROOTS
#endif
use Cntrl, only: NSTATE, PRCI, CITHR, IRREP, JBNAME, MLTPLT
use cntrl, only: iTOC15, LuIph
use Molcas, only: MxRoot
use Definitions, only: u6

implicit none
#ifdef _HDF5_
integer :: refwfn_id
integer :: root2state(mxroot), IDXCI
#endif
integer ISTATE
type(SGStruct) SGS
type(CIStruct) CIS
integer NCI
real*8 CI(NCI)
integer I, IAD, IDISK, JOB, LROOT1, LSYM

if ((ISTATE < 1) .or. (ISTATE > NSTATE)) then
  write(u6,*) 'RASSI/READCI: Invalid ISTATE parameter.'
  write(u6,*) ' ISTATE, NSTATE:',ISTATE,NSTATE
  call ABEND()
end if
JOB = JBNUM(ISTATE)
LROOT1 = LROOT(ISTATE)

#ifdef _HDF5_
if (mh5_is_hdf5(jbname(job))) then
  !*********************************************************************
  !
  ! For HDF5 formatted job files
  !
  !*********************************************************************
  refwfn_id = mh5_open_file_r(jbname(job))
  if (mh5_exists_attr(refwfn_id,'ROOT2STATE')) then
    call mh5_fetch_attr(refwfn_id,'ROOT2STATE',root2state)
    IDXCI = root2state(lroot1)
  else
    IDXCI = lroot1
  end if
  if ((IDXCI <= 0) .or. (IDXCI > NROOTS(JOB))) then
    call WarningMessage(2,'Invalid CI array index, abort!')
    call AbEnd()
  end if
  call mh5_fetch_dset(refwfn_id,'CI_VECTORS',CI,[NCI,1],[0,IDXCI-1])
  call mh5_close_file(refwfn_id)
else
#endif
  !*********************************************************************
  !
  ! For JOBIPH/JOBMIX formatted job files
  !
  !*********************************************************************
  call DANAME(LUIPH,JBNAME(JOB))
  IAD = 0
  call IDAFILE(LUIPH,2,ITOC15,30,IAD)
  IDISK = ITOC15(4)
  do I=1,LROOT1-1
    call DDAFILE(LUIPH,0,CI,NCI,IDISK)
  end do
  call DDAFILE(LUIPH,2,CI,NCI,IDISK)
  call DACLOS(LUIPH)
#ifdef _HDF5_
end if
!***********************************************************************
#endif

if ((IPGLOB > 0) .and. PRCI) then
  write(u6,*) ' READCI called for state ',ISTATE
  write(u6,*) ' This is on JobIph nr.',JOB
  write(u6,*) ' JobIph file name:',JBNAME(JOB)
  write(u6,*) ' It is root nr.',LROOT(ISTATE)
  write(u6,*) ' Its length NCI=',NCI
  write(u6,*) ' Its symmetry  =',IRREP(JOB)
  write(u6,*) ' Spin multiplic=',MLTPLT(JOB)
  LSYM = IRREP(JOB)
  call PRWF(SGS,CIS,LSYM,CI,CITHR)
end if

end subroutine READCI
