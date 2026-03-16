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
!***********************************************************************
!                                                                      *
!     Only read the number of states                                   *
!                                                                      *
!***********************************************************************

subroutine rdjob_nstates(JOB)

#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_fetch_attr, mh5_close_file
#endif
use Cntrl, only: NSTATE, ISTAT, NSTAT, JBNAME, LSYM1, NCONF1
use cntrl, only: NACTE1, MPLET1, NSYM1, NFRO1, NISH1, NASH1, NDEL1, NBAS1, NRS11, NRS21, NRS31, LROT1, NROOT1, IROOT1, NHOL11, &
                 NELE31, NAME, HEAD1, TITLE1
use cntrl, only: iTOC15, LuIph
use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
use RASDim, only: MxTit

implicit none
#ifdef _HDF5_
integer :: refwfn_id
integer :: ref_nstates
#endif
real*8 Weight(MxRoot), ENUCDUMMY
integer job, iad, ipt2

#ifdef _HDF5_
if (mh5_is_hdf5(jbname(job))) then
  !*********************************************************************
  !
  ! For HDF5 formatted job files
  !
  !*********************************************************************
  refwfn_id = mh5_open_file_r(jbname(job))
  call mh5_fetch_attr(refwfn_id,'NSTATES',ref_nstates)
  ! update the state offset, number of states, and total number of states
  ISTAT(JOB) = NSTATE+1
  NSTAT(JOB) = ref_nstates
  NSTATE = NSTATE+ref_nstates
  call mh5_close_file(refwfn_id)
else
#endif
  call DANAME(LUIPH,JBNAME(JOB))
  ! READ TABLE OF CONTENTS ON THIS JOBIPH FILE:
  IAD = 0
  call IDAFILE(LUIPH,2,ITOC15,30,IAD)
  ! SCATTER-READ VARIOUS DATA:
  IAD = ITOC15(1)
  call WR_RASSCF_Info(LUIPH,2,IAD,NACTE1,MPLET1,NSYM1,LSYM1,NFRO1,NISH1,NASH1,NDEL1,NBAS1,mxSym,NAME,(LenIn+8)*mxOrb,NCONF1,HEAD1, &
                      2*72,TITLE1,4*mxTit*18,ENUCDUMMY,LROT1,NROOT1,IROOT1,mxRoot,NRS11,NRS21,NRS31,NHOL11,NELE31,IPT2,Weight)
  ! update the state offset, number of states, and total number of states
  ISTAT(JOB) = NSTATE+1
  NSTAT(JOB) = NROOT1
  NSTATE = NSTATE+NROOT1
  call DACLOS(LUIPH)
#ifdef _HDF5_
end if
#endif

end subroutine rdjob_nstates
