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
* Copyright (C) 2015, Marcus Johansson                                 *
************************************************************************
      Subroutine fmsym_create_context(ctx)
      Integer ret
*     INT cmsym_create_context(msym_context *pctx, int *err)
      call cmsym_create_context(ctx,ret)
      if (ret.ne.0) then
         Call WarningMessage(2,'Faile to create symmetry context')
         Call Abend()
      end if
      Return
      End

      Subroutine fmsym_set_elements(ctx)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
      Character*2 Element(MxAtom)
      Call Get_nAtoms_All(nAtoms)
      Call Allocate_Work(ipCoord,3*nAtoms)
      Call Get_Coord_All(Work(ipCoord),nAtoms)
      Call Get_Name_All(Element)
      Call fmsym_set_ele_orb(ctx,nAtoms,Element,Work(ipCoord))
      Call Free_Work(ipCoord)
      Return
      End

      Subroutine fmsym_release_context(ctx)
      Integer ret
*     INT cmsym_release_context(msym_context *pctx, int*err)
      call cmsym_release_context(ctx,ret)
      Return
      End

      Subroutine fmsym_set_ele_orb(ctx,nAtoms,Element,Coord)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "periodic_table.fh"
#include "WrkSpc.fh"
#include "constants2.fh"
      Character*2 Element(nAtoms)
      Real*8 Coord(3,nAtoms)
      Character*(LENIN) AtomLabel(nAtoms)
      Integer basis_ids(4*mxBas)
      Dimension nBas(mxSym)
      Integer ret

      Call Get_LblCnt_All(AtomLabel)
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)

       if (nSym.ne.1) then
         Call WarningMessage(2,'MSYM can only be used with group c1')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      nCMO=0
      nMO=0
      do iSym=1, nSym
         nMO = nMO + nBas(iSym)
         nCMO = nCMO + nBas(iSym)**2
      end do

      Call Get_iArray('Basis IDs',basis_ids,(4*nMO))
*     INT cmsym_set_elements(msym_context *pctx, INT *pel, INT *puel, char *uelement, double xyz[][3], INT *paol, INT basis_ids[][4], int *err)
      call cmsym_set_elements(ctx,nAtoms,(LENIN),AtomLabel,Coord,
     & nMO, basis_ids, ret)
      call xflush(6)
      if (ret.ne.0) then
         Call WarningMessage(2,'Failed to set elements')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      End

      Subroutine fmsym_find_symmetry(ctx)
      Character*6 PGName
      Integer ret
*     INT cmsym_find_symmetry(msym_context *pctx, char pgname[6], int *err)
      Call cmsym_find_symmetry(ctx,PGName,ret)
      if (ret.ne.0) then
         Call WarningMessage(2,'Failed to find symmetry')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      Write (6,*) 'Found Point Group:'
      Write (6,*) PGName

      End

      Subroutine fmsym_symmetrize_molecule(ctx)
      Character*256 FN
      Integer ret
      Call PrgmTranslate('MSYMOUT',FN,lFN)
      FN = FN(1:lFN)
      FN(lFN+1:lFN+1)=Char(0)
*     INT cmsym_symmetrize_molecule(msym_context *pctx, char *outfile, INT *err)
      call cmsym_symmetrize_molecule(ctx, FN, ret)
      if (ret.ne.0) then
         Call WarningMessage(2,'Failed to symmetrize molecule')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      End

      Subroutine fmsym_generate_orbital_subspaces(ctx)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#include "mh5.fh"
      integer fileid, dsetid
#endif
      Character*80 Title
      Character(8), allocatable :: irrep_strings(:)
      Integer ret
      Dimension nBas(mxSym)

      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      nCMO=0
      nMO=0
      do iSym=1, nSym
         nMO = nMO + nBas(iSym)
         nCMO = nCMO + nBas(iSym)**2
      end do

      Call Allocate_Work(ipCAO,nCMO)
      Call Allocate_Work(ipOcc,nMO)
      Call Allocate_iWork(ipIrrIds,nMO)
      Call Allocate_iWork(ipIrrInd,nMO)
      call mma_allocate(irrep_strings, nMO)

*      INT cmsym_generate_orbital_subspaces(msym_context *pctx, INT *l, double c[*l][*l], INT irrep_ids[*l], INT irrep_ind[*l], INT *err)
      call cmsym_generate_orbital_subspaces(ctx,nMO,Work(ipCAO),
     &     iWork(ipIrrIds),iWork(ipIrrInd),nIrr,irrep_strings,ret)
      Write(6,*) 'Irrep indeces= '
      Write(6,'(5i3)') (iWork(ipIrrInd+k), k=0,nMO-1)
      Write(6,*) 'Irrep ids= '
      Write(6,'(5i3)') (iWork(ipIrrInd+k), k=0,nMO-1)

      if (ret.ne.0) then
         Call WarningMessage(2,'Failed to generate SALCs')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      Call FZero(Work(ipOcc),nMO)

      Title='Orbital Subspaces'
      Call WrVec('MSYMAORB',LuOrb,'CO',nSym,nBas,nBas,Work(ipCAO),
     &     Work(ipOcc),Dymmy,iWork(ipIrrInd),Title)
#ifdef _HDF5_
      fileid = mh5_create_file('MSYMH5')
      call run2h5_molinfo(fileid)
      call one2h5_ovlmat(fileid, nsym, nbas)
      call one2h5_fckint(fileid, nsym, nbas)
      ! mocoef
      dsetid = mh5_create_dset_real(fileid,'MO_VECTORS', 1, [nCMO])
      call mh5_init_attr(dsetid, 'description',
     $        'Coefficients of the SALCs as produced by MSYM, '//
     $        'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
      call mh5_put_dset(dsetid, Work(ipCAO))
      call mh5_close_dset(dsetid)
      ! supsym
      dsetid = mh5_create_dset_int(fileid,
     $ 'SUPSYM_IRREP_IDS', 1, [nMO])
      call mh5_init_attr(dsetid, 'description',
     $        'Super-symmetry ids as produced by MSYM, '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')
      call mh5_put_dset(dsetid, iWork(ipIrrIds))
      call mh5_close_dset(dsetid)
      dsetid = mh5_create_dset_int(fileid,
     $ 'SUPSYM_IRREP_INDICES', 1, [nMO])
      call mh5_init_attr(dsetid, 'description',
     $        'Super-symmetry indices as produced by MSYM, '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')
      call mh5_put_dset(dsetid, iWork(ipIrrInd))
      call mh5_close_dset(dsetid)
      ! irrep_labels
      dsetid = mh5_create_dset_str(fileid,
     $ 'SUPSYM_IRREP_LABELS', 1, [nIrr], 8)
      call mh5_init_attr(dsetid, 'description',
     $        'Super-symmetry labels as produced by MSYM, '//
     $        'arranged as array of size i=1,#supsym_irreps')
      call mh5_put_dset(dsetid, irrep_strings)
      call mh5_close_dset(dsetid)
      call mh5_close_file(fileid)
#endif

      call mma_deallocate(irrep_strings)
      Call Free_Work(ipCAO)
      Call Free_Work(ipOcc)
      Call Free_iWork(ipIrrInd)
      Call Free_iWork(ipIrrIds)

      End

      Subroutine fmsym_symmetrize_orbitals(ctx,CIO)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Integer ret
      Dimension CIO(*)
      Dimension nBas(mxSym)

      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)

      if (nSym.ne.1) then
         Call WarningMessage(2,'MSYM can only be used with group c1')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      nCMO=0
      nMO=0
      do iSym=1, nSym
         nMO = nMO + nBas(iSym)
         nCMO = nCMO + nBas(iSym)**2
      end do

*      INT cmsym_symmetrize_orbitals(msym_context *pctx, INT *l, double c[*l][*l], INT *err)
      call cmsym_symmetrize_orbitals(ctx,nMO,CIO,ret)
      if (ret.ne.0) then
         Call WarningMessage(2,'Failed to symmetrize orbitals')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      End

      Subroutine fmsym_symmetrize_orb_file(ctx,INPORB)

      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "constants2.fh"
      Integer ret
      Dimension nBas(mxSym)
      Character*80 Title
      Character INPORB*(*)


      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)

      if (nSym.ne.1) then
         Call WarningMessage(2,'MSYM can only be used with group c1')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if

      nCMO=0
      nMO=0
      do iSym=1, nSym
         nMO = nMO + nBas(iSym)
         nCMO = nCMO + nBas(iSym)**2
      end do

      Call Allocate_Work(ipCIO,nCMO)
      Call Allocate_Work(ipOcc,nMO)
      Call Allocate_Work(ipE,nMO)
      Call Allocate_iWork(iTIND,maxbfn)

      Call RdVec(INPORB,LuOrb,'COEI',nSym,nBas,nBas,Work(ipCIO),
     &     Work(ipOcc), Work(ipE) ,iWork(iTIND),Title,iWarn,iErr)

*     INT cmsym_symmetrize_orbitals(msym_context *pctx, INT *l, double c[*l][*l], INT *err)
      call cmsym_symmetrize_orbitals(ctx,nMO,Work(ipCIO),ret)

      if (ret.ne.0) then
         Call WarningMessage(2,'Failed to symmetrize orbitals')
*        INT cmsym_release_context(msym_context *pctx, int*err)
         call cmsym_release_context(ctx, ret)
         Call Abend()
      end if
      Title='Symmetrized Orbitals'
      Call WrVec('MSYMMORB',LuOrb,'COEI',nSym,nBas,nBas,Work(ipCIO),
     &     Work(ipOcc),Work(ipE),iWork(iTIND),Title)

      Call Free_Work(ipCIO)
      Call Free_Work(ipOcc)
      Call Free_Work(ipE)
      Call Free_iWork(iTIND)

      End
