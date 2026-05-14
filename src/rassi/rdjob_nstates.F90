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

use Cntrl, only: bNAME, HEAD1, IROOT1, ISTAT, iTOC15, JBNAME, LROT1, LSYM1, LuIph, MPLET1, NACTE1, NASH1, NBAS1, NCONF1, NDEL1, &
                 NELE31, NFRO1, NHOL11, NISH1, NROOT1, NRS11, NRS21, NRS31, NSTAT, NSTATE, NSYM1, TITLE1
use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
use RASDim, only: MxTit
#ifdef _HDF5_
use mh5, only: mh5_close_file, mh5_fetch_attr, mh5_is_hdf5, mh5_open_file_r
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: job
integer(kind=iwp) :: iad, ipt2
#ifdef _HDF5_
integer(kind=iwp) :: ref_nstates, refwfn_id
#endif
real(kind=wp) :: ENUCDUMMY
real(kind=wp), allocatable :: Weight(:)

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
  call mma_allocate(Weight,MxRoot,Label='Weight')
  call WR_RASSCF_Info(LUIPH,2,IAD,NACTE1,MPLET1,NSYM1,LSYM1,NFRO1,NISH1,NASH1,NDEL1,NBAS1,mxSym,bNAME,(LenIn+8)*mxOrb,NCONF1, &
                      HEAD1,2*72,TITLE1,4*mxTit*18,ENUCDUMMY,LROT1,NROOT1,IROOT1,mxRoot,NRS11,NRS21,NRS31,NHOL11,NELE31,IPT2,Weight)
  call mma_deallocate(Weight)
  ! update the state offset, number of states, and total number of states
  ISTAT(JOB) = NSTATE+1
  NSTAT(JOB) = NROOT1
  NSTATE = NSTATE+NROOT1
  call DACLOS(LUIPH)
#ifdef _HDF5_
end if
#endif

end subroutine rdjob_nstates
