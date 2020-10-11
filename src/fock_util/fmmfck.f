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
* Copyright (C) 2007, Mark A. Watson                                   *
************************************************************************
      Subroutine FMMFck(Dens,TwoHam,ndim)
#ifdef _NOT_ACTIVE_
************************************************************************
*                                                                      *
*     purpose: Generate FMM interface file and call FMM driver         *
*              to update Fock matrix with multipole-derived            *
*              Coulomb matrix elements                                 *
*                                                                      *
*     called from: Drv2El_dscf                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by: Mark A. Watson (maw)                                 *
*     University of Tokyo, 2007                                        *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*

#include "mxdm.fh"
#include "real.fh"
      Parameter(LMAX = 12)
      Real*8 Dens(ndim), TwoHam(ndim), nBas(8)
*
*---- Define local variables
      Real*8 CarMoms( ndim, (LMAX+1)*(LMAX+2)/2 , LMAX+1 )
      Real*8 SphMoms( ndim, 2*LMAX+1 , LMAX+1 )
      Real*8 Moms_batch( ndim+4 )
      Real*8 CntrX(ndim+4), CntrY(ndim+4), CntrZ(ndim+4)
      Character*8 Label
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
*
*---- Compute lengths of matrices
      lDens = 0
      nBasTot = 0
      Do iSym = 1, nSym
         lDens  = lDens  + nBas(iSym)*(nBas(iSym) + 1)/2
         nBasTot = nBasTot + nBas(iSym)
      End Do
      If (.not. lDens .eq. ndim) Then
         Write (6,*) 'ERROR in FMMFck', lDens, ndim
         Call Abend()
      End If

C      call dcopy_(ndim,Zero,0,Moments,1)
C      call dcopy_(ndim,Zero,0,CntrX,1)
C      call dcopy_(ndim,Zero,0,CntrY,1)
C      call dcopy_(ndim,Zero,0,CntrZ,1)
*
*---- Read centres
*
      iRc=-1
      iOpt=2
      iComp=1
      iSyLbl=1
      Label='FMMCnX'
      Call RdOne(iRc,iOpt,Label,iComp,CntrX,iSyLbl)
      If (iRc.ne.0) Then
         Write (6,*) 'FMMFck: Error readin ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      Label='FMMCnY'
      Call RdOne(iRc,iOpt,Label,iComp,CntrY,iSyLbl)
      If (iRc.ne.0) Then
         Write (6,*) 'FMMFck: Error readin ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      Label='FMMCnZ'
      Call RdOne(iRc,iOpt,Label,iComp,CntrZ,iSyLbl)
      If (iRc.ne.0) Then
         Write (6,*) 'FMMFck: Error readin ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
*
*---- Read moments from one-electron files
*
      Do L = 0, LMAX
      Do iComp = 1, (L+1)*(L+2)/2
         iRc=-1
         iOpt=2
         iSyLbl=1
         Write (Label,'(A,I2)') 'FMMInt', L
         Call RdOne(iRc,iOpt,Label,iComp,Moms_batch,iSyLbl)
         If (iRc.ne.0) Then
            Write (6,*) 'FMMFck: Error readin ONEINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
         Do ij = 1, ndim
            CarMoms(ij,iComp,L+1) = Moms_batch(ij)
         End Do
      End Do
      End Do
*
*---- Transform cartesian to spherical components
*
#ifndef  _NO_F90_COMPILER_
C     CALL mm_call_car_to_sph(CarMoms,SphMoms,ndim,LMAX)
#endif
*
*---- Write to FMM interface file
*
*     Write array lengths in header file
      OPEN(98,FILE='MM_DATA_HEADER',FORM='UNFORMATTED',STATUS='REPLACE')
      WRITE (98) LMAX, nBasTot, ndim, 0
      CLOSE(98,STATUS='KEEP')

*     Write multipole moments and density information
      OPEN(98,FILE='MM_DATA',FORM='UNFORMATTED',STATUS='REPLACE')

      ij = 0
      Do J = 1, nBasTot
      Do I = 1, J
         ij = ij+1
         Do L = 0, LMAX
            Do M = -L, L
               iM = M+L+1
C               WRITE (6,'(5I3,2X, 3F10.6,2E15.4)') L,M, ij,1,ij,
C     &                    CntrX(ij), CntrY(ij), CntrZ(ij),
C     &                    SphMoms(ij,iM,L+1), Dens(ij)
               WRITE (98) L,M, I,J,ij,
     &                    CntrX(ij), CntrY(ij), CntrZ(ij),
     &                    SphMoms(ij,iM,L+1), Dens(ij)
            End Do
         End Do
      End Do
      End Do

*     Mark end of file with negative angular momentum
      WRITE (98) -1,0, 0,0,0, 0d0,0d0,0d0, 0d0,0d0
      CLOSE(98,STATUS='KEEP')
*
*---- Now call multipole code to update the Fock matrix with the
*     long-range multipole-computed Coulomb matrix elements.
*
#ifndef  _NO_F90_COMPILER_
C     CALL mm_call_get_J_matrix(TwoHam,ndim,nBasTot,LMAX)
#endif
*
*     Coulomb contributions of TwoHam should now be complete!
*
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(Dens)
         Call Unused_integer(ndim)
      End If
#endif
c Avoid unused argument warnings
      If (.False.) Call Unused_real(TwoHam)
      Return
      End
