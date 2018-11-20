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
      Subroutine pseudospin(M,dim,z,iDir,iOpt,iprint)
c dim - dimension (input)
c moment(l,dim,dim) (input)
c z - pseuDospin eigenfunctions (output)
      Implicit None
#include "stdalloc.fh"
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: dim, iprint
      Integer, intent(in)           :: iDir, iOpt
      Complex(kind=wp), intent(in)  :: M(3,dim,dim)
      Complex(kind=wp), intent(out) :: z(dim,dim)
      ! local variables:
      Integer                       :: info, i
      Real(kind=wp), allocatable    :: w(:)
      Complex(kind=wp), allocatable :: z1(:,:)
      Real(kind=wp)                 :: dznrm2_
      External                      :: dznrm2_
      Logical                       :: dbg
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Call qEnter('pseudospin')

      Call mma_allocate(W,dim,'W')
      Call mma_allocate(Z1,dim,dim,'Z1')
      dbg=iprint.ge.3
      Call dcopy_(dim,0.0_wp,0,W,1)
      Call zcopy_(dim*dim,(0.0_wp,0.0_wp),0,Z,1)
      Call zcopy_(dim*dim,(0.0_wp,0.0_wp),0,Z1,1)
      info=0
      Call diag_c2(M(iDir,1:dim,1:dim),dim,info,w,z1)
      If(dbg) Then
        Do i=1,dim
          Write(6,'(A,i3,A,F24.14)') 'i=',i,' eigenvalue=',w(i)
        End Do
      End If
      If (info.ne.0) Then
        Write(6,'(5x,a)') 'PSEUDO::  diagonalization of the '//
     &                    'zeeman hamiltonian failed.'
        Go To 199
      End If
      If(dbg) Write(6,*) "PSEUDO:  norm of  M is:",
     &                                            dznrm2_(3*dim*dim,M,1)
      If(dbg) Write(6,*) "PSEUDO:  norm of Z1 is:",dznrm2_(dim*dim,Z1,1)
      If(iDir==3) Then
         If(iOpt==1) Then
            Call spin_phase(M,dim,z1,z)
         Else
            Call zcopy_(dim*dim,z1,1,z,1)
            Write(6,*) 'PSEUDOSPIN:  iOpt = ',iOpt
            Call WarningMessage(2,'PSEUDOSPIN: iOpt '//
     &                            'is not understood.')
         End If
      Else
         Call zcopy_(dim*dim,z1,1,z,1)
      End If

 199  Continue
      Call mma_deallocate(W)
      Call mma_deallocate(Z1)
      Call qExit('pseudospin')
      Return
      End




