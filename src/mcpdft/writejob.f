************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2020, Chen Zhou                                        *
************************************************************************
      subroutine writejobms(iadr19,L1,L2)

#ifdef _HDF5_
      use mh5, only: mh5_open_file_rw, mh5_open_dset,
     &               mh5_put_dset, mh5_close_file
#endif

      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "rasrc.fh"
#include "general.fh"
#include "wjob.fh"
      real*8 :: U(lroots,lroots)
      integer :: i,j,idisk,LCIRot,LW4,LEnergy,L1,L2
      integer :: iadr19(15),dum(1)
#ifdef _HDF5_
      integer :: refwfn_id,wfn_energy,wfn_cicoef
#endif

      CALL GETMEM('Energy','ALLO','REAL',LEnergy,mxRoot*mxIter)
      call dcopy_(mxRoot*mxIter,[0.0d0],0,Work(LEnergy),1)
      iDisk = IADR19(6)
      Do i = 1,mxIter
       Do j = 1,lroots
        work(LEnergy+mxRoot*(i-1)+j-1)=work(L1+j-1)
       End do
      End do

      Call DDafile(JOBIPH,1,Work(LEnergy),mxRoot*mxIter,iDisk)
#ifdef _HDF5_
      refwfn_id = mh5_open_file_rw('RASWFN')
      wfn_energy = mh5_open_dset(refwfn_id,'ROOT_ENERGIES')
      call mh5_put_dset(wfn_energy,ener(1,1))
#endif
      CALL GETMEM('Energy','FREE','REAL',LEnergy,mxRoot*mxIter)
      idisk = 284
      Call iDafile(JOBIPH,2,dum,1,iDisk)
      ncon=dum(1)

      Do i = 1,lroots
       Do j = 1,lroots
        U(j,i) = work(L2+(i-1)*lroots+j-1)
       End do
      End do

      CALL GETMEM('CIRot','ALLO','REAL',LCIRot,NCON*LROOTS)
      call dcopy_(NCON*LROOTS,[0.0d0],0,Work(LCIRot),1)
      CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCON)
      iDisk = IADR19(4)
      DO I=1,LROOTS
       Call DDafile(JOBIPH,2,Work(LW4),nCon,iDisk)
       DO J=1,LROOTS
        DO K=1,NCON
         ITMP=LCIRot+(J-1)*NCON+K-1
         WORK(ITMP)=WORK(ITMP)+WORK(LW4+K-1)*U(I,J)
        END DO
       END DO
      END DO

      iDisk = IADR19(4)
      DO I=1,LROOTS
       Call DDafile(JOBIPH,1,Work(LCIRot+(I-1)*NCON),nCon,iDisk)
      END DO
#ifdef _HDF5_
      wfn_cicoef = mh5_open_dset(refwfn_id,'CI_VECTORS')
      Do I = 1,LROOTS
          call mh5_put_dset(wfn_cicoef,
     $ Work(LCIRot+NCON*(I-1):LCIRot+NCON*I-1),[NCON,1],[0,i-1])
      END DO
      call mh5_close_file(refwfn_id)
#endif

      CALL GETMEM('CIVEC','FREE','REAL',LW4,NCON)
      CALL GETMEM('CIRot','FREE','REAL',LCIRot,NCON*LROOTS)

      End

      subroutine writejob(iadr19)

#ifdef _HDF5_
      use mh5, only: mh5_open_file_rw, mh5_open_dset,
     &               mh5_put_dset, mh5_close_file
#endif

      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "rasscf.fh"
#include "rasrc.fh"
#include "general.fh"
#include "wjob.fh"
      integer :: i,j,idisk,LEnergy
      integer :: iadr19(15)
#ifdef _HDF5_
      integer :: refwfn_id,wfn_energy
#endif

      CALL GETMEM('Energy','ALLO','REAL',LEnergy,mxRoot*mxIter)
      call dcopy_(mxRoot*mxIter,[0.0d0],0,Work(LEnergy),1)
      Do i = 1,mxIter
       Do j = 1,lroots
        work(LEnergy+mxRoot*(i-1)+j-1)=ener(j,1)
       End do
      End do

      iDisk = IADR19(6)
      Call DDafile(JOBIPH,1,Work(LEnergy),mxRoot*mxIter,iDisk)
#ifdef _HDF5_
      refwfn_id = mh5_open_file_rw('RASWFN')
      wfn_energy = mh5_open_dset(refwfn_id,'ROOT_ENERGIES')
      call mh5_put_dset(wfn_energy,ener(1,1))
      call mh5_close_file(refwfn_id)
#endif
      CALL GETMEM('Energy','FREE','REAL',LEnergy,mxRoot*mxIter)
      End
