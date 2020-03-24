************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine OccDef(Occ,mmB,nD,CMO,mBB)
#include "compiler_features.h"
#ifndef POINTER_REMAP
      Use, Intrinsic :: ISO_C_BINDING
#endif
      use OccSets
      use Orb_Type
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
      Real*8 Occ(mmB,nD)
      Real*8, Target::  CMO(mBB,nD)
      Real*8, Pointer, Dimension(:,:):: pCMO
      Real*8, Dimension(:), Allocatable:: OccTmp
      Integer, Allocatable:: iFerm(:)
*
*---- Form occupation numbers
*
*     Note, this puts only in the occupation numbers of the electronic
*     orbitals. The muonic occupation numbers are put in place in
*     NewOrb_Scf on the first iteration.
*
      If (OnlyProp) Return
*
*define _DEBUG_
#ifdef _DEBUG_
      Do iD = 1, nD
         Write (6,*) 'iD=',iD
         Write(6,'(a,8i5)') 'sorb: nOcc   ',(nOcc(i,iD),i=1,nSym)
      End Do
      If (Allocated(OccSet_e)) Then
         Do iD = 1, nD
            Write (6,*) 'iD=',iD
            iOff=1
            Do iSym = 1, nSym
               Call RecPrt('OccSet_e',' ',OccSet_e(iOff,iD),1,
     &                       nOcc(iSym,iD))
               iOff = iOff + nOcc(iSym,iD)
            End Do
         End Do
      End If
      If (Allocated(OccSet_m)) Then
         Do iD = 1, nD
            Write (6,*) 'iD=',iD
            iOff=1
            Do iSym = 1, nSym
               Call RecPrt('OccSet_m',' ',OccSet_m(iOff,iD),1,
     &                       nOcc(iSym,iD))
               iOff = iOff + nOcc(iSym,iD)
            End Do
         End Do
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     First we set up the electronic occupation numbers either as
*     by default or as specied by the user. The numbers are stored
*     in the array Occ.
*
      Call FZero(Occ,mmB*nD)
*
      If (Allocated(OccSet_e)) Go To 100
*
*     Deafault is that occupied orbitals are set to 2 and 1
*     for RHF and UHF, respectively.
*
      Do iD = 1, nD
         mOcc=0
         Do iSym = 1,nSym
            Do iOrb = 1, nOcc(iSym,iD)
               Occ(iOrb+mOcc,iD)=DBLE(3-nD) ! Value 2 or 1
            End Do
            mOcc=mOcc+nOrb(iSym)
         End Do
      End Do
      Go To 200
*
 100  Continue
*
*     User define occupation numbers according to the OCCN key word.
*
      Do iD = 1, nD
         mOcc=0
         mSet=0
         Do iSym = 1,nSym
            Do iOrb = 1, nOcc(iSym,iD)
               Occ(iOrb+mOcc,iD)=OccSet_e(iOrb+mSet,iD)
            End Do
            mOcc=mOcc+nOrb(iSym)
            mSet=mSet+nOcc(iSym,iD)
         End Do
      End Do
      call mma_deallocate(OccSet_e)
*
 200  Continue
*                                                                      *
************************************************************************
*                                                                      *
*     Here we put in the occupations from the muons
*
*     Here we will have to run through the orbitals and identify
*     if the orbitals are muonic or electronic. As we proceed we
*     will fill up the vector with the occupation number from
*     the list of either muonic or electronic occupation numbers.
*
      If (Allocated(OccSet_m)) Then
*
         Call mma_Allocate(OccTmp,mmB,Label='OccTmp')
         Call mma_allocate(iFerm,nnB,Label='iFerm')
         Call Get_iArray('Fermion IDs',iFerm,nnB)
*        Write (6,*) 'iFerm=',iFerm
*
         Do iD = 1, nD
*
*           Store the electronic occupation numbers in OccTmp
*
            Call DCopy_(mmB,Occ(1,iD),1,OccTmp,1)
*           Write (6,*) 'OccTmp=',OccTmp
*           Write (6,*) 'Occset_m=',Occset_m
            Call FZero(Occ(1,iD),mmB)
*
            nOcc_e=0     ! number of occupied electronic orbitals
            nOcc_m=0     ! number of occupied muonic orbitals
            iOff = 1     ! Offset to the symmetry sections of CMO
            jEOr = 0     ! Offset to the array with occupation numbers
            Do iSym = 1, nSym
*
*              map the relevant portion of CMO onto pCMO
*
               nB = nBas(iSym)
#ifdef POINTER_REMAP
               pCMO(1:nB,1:nB) => CMO(iOff:iOff+nB**2-1,iD)
#else
               Call C_F_POINTER(C_LOC(CMO(iOff,iD)), pCMO, [nB,nB])
#endif
*
               Do iOrb = 1, nOrb(iSym)  ! Loop over the orbitals
*
*                 Compute identifier which is different from zero
*                 if the orbital is muonic.
*
                  tmp=0.0D0
                  Do k = 1, nB
                     tmp = tmp + DBLE(iFerm(jEOr+k))
     &                         * ABS(pCMO(k,iOrb))
                  End Do
                  Muon_i=0                    ! electronic
                  If (tmp.ne.0.0D0) Muon_i= 1 ! muonic
*
                  If (Muon_i.eq.0) Then
*
*                    The orbital is electronic, put in an
*                    electronic occupation number. If the index
*                    is beyond the number of occupied orbitals
*                    put in a Zero.
*
                     nOcc_e = nOcc_e + 1
                     If (nOcc_e.le.nOcc(iSym,iD)) Then
                        Occ(jEor+iOrb,iD) = OccTmp(jEOr+nOcc_e)
                     Else
                        Occ(jEor+iOrb,iD) = 0.0D0
                     End If
*                    Write (6,*) 'Electronic:',iOrb,Occ(jEOr+iOrb,iD)
*
                  Else If (Muon_i.eq.1) Then
*
*                    The orbital is muonic, put in an
*                    muonic occupation number. If the index
*                    is beyond the number of occupied orbitals
*                    put in a Zero.
*
                     nOcc_m = nOcc_m + 1
                     If (nOcc_m.le.nOcc(iSym,id)) Then
                        Occ(jEor+iOrb,iD) = OccSet_m(jEOr+nOcc_m,iD)
                     Else
                        Occ(jEor+iOrb,iD) = 0.0D0
                     End If
                     OrbType(jEor+iOrb,iD) = 1
*                    Write (6,*) 'Muonic:',iOrb,Occ(jEOr+iOrb,iD)
*
                  End If
*
               End Do
*
               Nullify(pCMO)
               iOff = iOff + nB**2
               jEOr = jEOr + nOrb(iSym)
            End Do ! iSym
*
         End Do    ! iD
*define _SPECIAL_DEBUG_
#ifdef _SPECIAL_DEBUG_
         Call DebugCMO(CMO,mBB,nD,Occ,mmB,nBas,nOrb,nSym,iFerm,
     &                 '@ the last position')
#endif
         Call mma_deallocate(iFerm)
         Call mma_deAllocate(OccTmp)
         call mma_deallocate(OccSet_m)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Finally we will have to resort the orbitals with respect to their
*     occupation numbers such that the orbitals which are formally in
*     the occupied space but who are empty are reassigned to being
*     virtuall orbitals.
*
*     This means that the the CMO and Occ arrays will be resorted, and
*     that the nOcc array will be updated with respect to the actual
*     number of orbitals with occupation different from zero, i.e.
*     they are virtual.
*
      Do iD = 1, nD
*define _DEBUG_
#ifdef _DEBUG_
         Write (6,*) 'iD=',iD
         Write (6,*) 'nOccs(original):'
         Write (6,*) (nOcc(iSym,iD),iSym=1,nSym)
         Write (6,*) 'nOV=',nOV
#endif
         iOff=1
         jOff=0
*
         Do iSym = 1, nSym
*
            If (nOrb(iSym).eq.0) Cycle
            iOcc = 0
            Do iOrb = 1, nOrb(iSym)-1
               jOrb = iOrb + 1
               If (Occ(jOrb+jOff,iD).gt.Occ(iOrb+jOff,iD)) Then
*
                  iTmp=OrbType(iOrb+jOff,iD)
                  OrbType(iOrb+jOff,iD)=OrbType(jOrb+jOff,iD)
                  OrbType(jOrb+jOff,iD)=iTmp
*
                  Tmp=Occ(iOrb+jOff,iD)
                  Occ(iOrb+jOff,iD)=Occ(jOrb+jOff,iD)
                  Occ(jOrb+jOff,iD)=Tmp
                  Call DSwap_(nBas(iSym),
     &                        CMO(iOff+(iOrb-1)*nBas(iSym),iD),1,
     &                        CMO(iOff+(jOrb-1)*nBas(iSym),iD),1)
               End If
               If (Occ(iOrb+jOff,iD).ne.0.0D0) iOcc=iOcc+1
            End Do
            If (Occ(nOrb(iSym)+jOff,iD).ne.0.0D0) iOcc=iOcc+1
            nOcc(iSym,iD)=iOcc
*
            jOff=jOff+nOrb(iSym)
            iOff=iOff+nBas(iSym)*nOrb(iSym)
         End Do
#ifdef _DEBUG_
         Write (6,*) 'iD=',iD
         Write (6,*) 'nOccs(new):'
         Write (6,*) (nOcc(iSym,iD),iSym=1,nSym)
#endif
      End Do
*
*     Sort such that the occupied orbitals are first in each irrep.
*
      Call Sorb2CMOs(CMO,mBB,nD,Occ,mmB,nBas,nOrb,nSym,OrbType)
*                                                                      *
************************************************************************
*                                                                      *
*     Recompute sizes since the nOcc array might have changed.
*
      nOV = 0
      Do iSym = 1, nSym
         If (nD.eq.1) Then
             maxnOcc=nOcc(iSym,1)
             minnOcc=nOcc(iSym,1)
         Else
             maxnOcc=max(nOcc(iSym,1),nOcc(iSym,2))
             minnOcc=min(nOcc(iSym,1),nOcc(iSym,2))
         End If
         nOV    = nOV  + (maxnOcc-nFro(iSym))*
     &                   (nOrb(iSym)-minnOcc)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
#ifdef _DEBUG_
      Do iD = 1, nD
         iOff=1
         jOff=1
         Do iSym = 1, nSym
            Call RecPrt('Occ','(10F6.2)',Occ(iOff,iD),1,nOrb(iSym))
            Call RecPrt('CMO','(10F6.2)',CMO(jOff,iD),
     &                  nBas(iSym),nOrb(iSym))
            Do i = 0, nOrb(iSym)-1
               Write (6,*) 'i,OrbType=',i+1, OrbType(iOff+i,iD)
            End Do
            iOff = iOff + nOrb(iSym)
            jOff = jOff + nOrb(iSym)*nBas(iSym)
         End Do
      End Do
#endif
      Return
      End
#ifdef _SPECIAL_DEBUG_
      Subroutine DebugCMO(CMO,nCMO,nD,Occ,nnB,nBas,nOrb,nSym,iFerm,
     &                    Label)
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(nCMO,nD), Occ(nnB,nD)
      Integer nBas(nSym),nOrb(nSym), iFerm(nnB)
      Character*(*) Label
*
      Do iD = 1, nD
         Write (6,*)
         If (iD.eq.1) Then
            If (nD.eq.1) Then
               Write (6,*) ' RHF CMOs'
            Else
               Write (6,*) ' UHF alpha CMOs'
            End If
         Else
            Write (6,*) ' UHF beta CMOs'
         End If
         Write (6,*)
         jOff=0
         iOff=1
         Do iSym = 1, nSym
            Do iOrb = 1, nOrb(iSym)
               tmp=0.0D0
               Do k = 1, nBas(iSym)
                  tmp = tmp + DBLE(iFerm(jOff+k))
     &                      * ABS(CMO(iOff-1+(iOrb-1)*nBas(iSym)+k,iD))
               End Do
               Write (6,*)
               If (tmp.ne.0.0D0) Then
                  Write (6,*) 'Muonic Orbital:',iOrb
               Else
                  Write (6,*) 'Electronic  Orbital:',iOrb
               End If
               Write (6,*) 'Occupation number:',Occ(jOff+iOrb,iD)
               Call RecPrt('CMO',' ',CMO(iOff+(iOrb-1)*nBas(iSym),iD),
     &                      1,nBas(iSym))
            End Do
            jOff = jOff + nOrb(iSym)
            iOff = iOff + nBas(iSym)*nOrb(iSym)
         End Do
      End Do
*
      Return
      End
#endif
