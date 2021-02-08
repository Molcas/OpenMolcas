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
      Subroutine fetch_data_RunFile_init( nss, nstate )
      Implicit None

      Integer :: nss, nstate, njob, mxjob
      Logical :: FOUND



      ! check the presence of the RUNFILE
      FOUND=.false.
      Call f_inquire('RUNFILE',FOUND)
      If(FOUND.eqv..false.) Then
         Write(6,'(5X,A)') 'The RUNFILE was not found in the $WorkDir'
         Write(6,'(5X,A)') 'Are you running this calculation in the'
         Write(6,'(5X,A)') 'same $WorkDir where RASSI calculation'
         Write(6,'(5X,A)') 'was executed?'
         Write(6,'(5X,A)') 'Check your calculation again, and If'
         Write(6,'(5X,A)') 'necessary, submit a BUG report.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      ! check the presence of the necessary information on RUNFILE:
      FOUND= .FALSE.
      Call qpg_iscalar('NSS_SINGLE',FOUND)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The NSS Value was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_iscalar('NJOB_SINGLE',FOUND)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The NJOB Value was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_iscalar('MXJOB_SINGLE',FOUND)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The MXJOB Value was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_iscalar('NSTATE_SINGLE',FOUND)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The NSTATE Value was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      ! fetch scalar data from RUNFILE
      Call get_iscalar('NSS_SINGLE',NSS)
      Call get_iscalar('NJOB_SINGLE',NJOB)
      Call get_iscalar('MXJOB_SINGLE',MXJOB)
      Call get_iscalar('NSTATE_SINGLE',NSTATE)

      ! check the presence of saved arrays on RUNFILE:

      FOUND=.FALSE.
      Call qpg_iArray('MLTP_SINGLE',FOUND,MXJOB)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The MLTP array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_iArray('JBNUM_SINGLE',FOUND,NSTATE)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The JBNUM array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_iArray('LROOT_SINGLE',FOUND,NSTATE)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The LROOT array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_dArray('ESO_SINGLE',FOUND,NSS)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The ESO array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_dArray('UMATR_SINGLE',FOUND,NSS**2)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The UMATR array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_dArray('UMATI_SINGLE',FOUND,NSS**2)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The UMATI array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Please report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_dArray('ANGM_SINGLE',FOUND,3*NSTATE*NSTATE)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The ANGMOM array was not found on RUNFILE'
         Write(6,'(5X,A)') '1. Check If ANGM keyword was used for '//
     &                     'SEWARD.'
         Write(6,'(5X,A)') '2. Check If MEES keyword was used for '//
     &                     'RASSI.'
         Write(6,'(5X,A)') '3. Check If PROP keyword was used for '//
     &                     'RASSI:'
         Write(6,'(9X,A)') 'PROP'
         Write(6,'(9X,A)') '3'
         Write(6,'(9X,A)') ' ''ANGMOM'' 1'
         Write(6,'(9X,A)') ' ''ANGMOM'' 2'
         Write(6,'(9X,A)') ' ''ANGMOM'' 3'
         Write(6,'(5X,A)') 'If MEES, ANGMOM and PROP keywords  were '//
     &                     'used and you still see this problem,'
         Write(6,'(5X,A)') 'please, report a BUG.'
         Call xFlush(6)
         Call Quit_OnUserError()
      End If

      FOUND=.FALSE.
      Call qpg_dArray('DIP1_SINGLE',FOUND,3*NSTATE*NSTATE)
      If (FOUND.EQV..FALSE.) Then
         Write(6,'(5X,A)') 'The DIPMOM array was not found on RUNFILE'
         Write(6,'(5X,A)') 'Absorption intensities will not be computed'
      End If
      Return
      End Subroutine fetch_data_RunFile_init




      Subroutine fetch_data_RunFile_all( nss, nstate,
     &                                   multiplicity, eso, esfs,
     &                                   U, MM, MS, ML, DM, angmom,
     &                                   eDmom, amfi, HSO,
     &                                   eso_au, esfs_au )
      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer :: nss, nstate
      Integer :: multiplicity(nstate)
      Real(kind=8) :: eso(nss), esfs(nstate), angmom(3,nstate,nstate),
     &                 eDmom(3,nstate,nstate), amfi(3,nstate,nstate),
     &                 eso_au(nss), esfs_au(nstate)
      Complex(kind=8) :: MM(3,nss,nss), MS(3,nss,nss), ML(3,nss,nss)
      Complex(kind=8) :: DM(3,nss,nss)
      Complex(kind=8) :: U(nss,nss), HSO(nss,nss)
      ! local variables:
      Integer              :: njob, mxjob, iss, ibas(nstate,-50:50)
      Integer              :: i, j, i1, j1, ist, jst, mult, multI, multJ
      Integer              :: l, ipar, info
      Real(kind=8)         :: g_e, au2cm, thr_deg, diff
      ! allocatable local arrays:
      Integer, allocatable :: mltplt(:), jbnum(:), nstat(:) !,lroot(:)
      Real(kind=8), allocatable :: tmpR(:,:), tmpI(:,:), W(:)
      Complex(kind=8), allocatable :: tmp(:,:)
      Complex(kind=8)  :: Spin
      External         :: Spin
      Logical          :: found_edmom, found_amfi, found_hsor,
     &                    found_hsoi

      g_e=2.00231930437180_wp
      au2cm=219474.6313702_wp
      ! get basic sizes:
      njob=0
      mxjob=0
      Call get_iScalar('NJOB_SINGLE',njob)
      Call get_iScalar('MXJOB_SINGLE',mxjob)
      ! allocate temporary memory:
      Call mma_allocate(jbnum,nstate,'jbnum')
      Call mma_allocate(mltplt,mxjob,'mltplt')
      Call mma_allocate(nstat,mxjob,'nstat')
      ! get the information from RUNFILE:
      mltplt=0
      jbnum=0
      nstat=0
      Call get_iArray('MLTP_SINGLE',MLTPLT,MXJOB)
      Call get_iArray('JBNUM_SINGLE',JBNUM,NSTATE)
      Call get_iArray('NSTAT_SINGLE',NSTAT,MXJOB)
      ! computing the multiplicity of each state:
      multiplicity=0
      Do i=1,nstate
         multiplicity(i) = mltplt( jbnum(i) )
      End Do
      Call mma_deallocate(jbnum)
      Call mma_deallocate(mltplt)
      Call mma_deallocate(nstat)

      ! fetch the spin-orbit energies:
      eso=0.0_wp
      eso_au=0.0_wp
      Call get_dArray('ESO_SINGLE',eso,nss)
      Call get_dArray('ESO_LOW',eso_au,nss)

      ! fetch the spin-free energies:
      esfs=0.0_wp
      esfs_au=0.0_wp
      Call get_dArray('ESFS_SINGLE',esfs,nstate)
      Call get_dArray('ESFS_SINGLEAU',esfs_au,nstate)

      ! fetch the U matrix:
      Call mma_allocate(tmpR,nss,nss,'tmpR')
      Call mma_allocate(tmpI,nss,nss,'tmpI')


      tmpR=0.0_wp
      tmpI=0.0_wp
      Call get_dArray('UMATR_SINGLE',tmpR,nss*nss)
      Call get_dArray('UMATI_SINGLE',tmpI,nss*nss)
      U=(0.0_wp,0.0_wp)
      Do i=1,nss
         Do j=1,nss
            U(i,j) = cmplx( tmpR(i,j), tmpI(i,j), wp )
         End Do
      End Do


      ! fetch the angular momentum integrals:
      angmom=0.0_wp
      Call get_dArray('ANGM_SINGLE',angmom,3*nstate*nstate)

      ! fetch the electric dipole moment integrals:
      edmom=0.0_wp
      found_EDMOM=.false.
      Call qpg_dArray('DIP1_SINGLE',FOUND_EDMOM,3*NSTATE*NSTATE)
      If (found_edmom)
     &  Call get_dArray('DIP1_SINGLE',edmom,3*nstate*nstate)


      ! fetch the amfi integrals:
      amfi=0.0_wp
      found_AMFI=.false.
      Call qpg_dArray('AMFI_SINGLE',FOUND_AMFI,3*NSTATE*NSTATE)
      If (found_amfi)
     &  Call get_dArray('AMFI_SINGLE',amfi,3*nstate*nstate)

      ! fetch the spin-orbit hamiltonian
      FOUND_HSOR=.FALSE.
      FOUND_HSOI=.FALSE.
      Call qpg_dArray('HAMSOR_SINGLE',FOUND_HSOR,NSS*NSS)
      Call qpg_dArray('HAMSOI_SINGLE',FOUND_HSOI,NSS*NSS)
      If (FOUND_HSOR.AND.FOUND_HSOI) Then
         tmpR=0.0_wp
         tmpI=0.0_wp
         Call get_dArray('HAMSOR_SINGLE',tmpR,nss*nss)
         Call get_dArray('HAMSOI_SINGLE',tmpI,nss*nss)
         Call zcopy_(nss*nss,[(0.0_wp,0.0_wp)],0,HSO,1)
         Do i=1,nss
            Do j=1,nss
               HSO(i,j) = cmplx( tmpR(i,j), tmpI(i,j), wp )
            End Do
         End Do
!-----------------------------------------------------------------------
!       if HSO is found, proceed to diagonalize it
         Call mma_allocate(W,nss,'W')
         Call dcopy_(nss,[0.0_wp],0,W,1)
         Call zcopy_(nss*nss,[(0.0_wp,0.0_wp)],0,U,1)
         info=0
         Call diag_c2(hso,nss,info,W,U)
         ! correct for numerical degeneracies:
         thr_deg=0.2D-13 ! a.u. = 0.2D-13*au2cm = 4.38949263E-09 cm-1
         Do i=1,nss-1
!           wtmp=W(i)
           Do j=i+1,nss
             diff=ABS( w(i)-w(j) )
             If( diff < thr_deg ) Then
               w(j)=w(i)
             End If
           End Do
         End Do

         Do i=1,nss
            ESO(i) = (W(i)-W(1))*au2cm
         End Do

         Call mma_deallocate(W)
!-----------------------------------------------------------------------
      End If
      Call mma_deallocate(tmpR)
      Call mma_deallocate(tmpI)


c-----
      ! generate a local indexing table:
      iss=0
      ibas=0
      ipar=mod(multiplicity(1),2)
      Do Ist=1,nstate
         Mult=Multiplicity(Ist)
         Do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            If( (Ipar==0) .AND. (I==0)) Go To 310
               Iss=Iss+1
               Ibas(Ist,I)=Iss
  310       Continue
         End Do ! i
      End Do ! ist

c----- expand the spin free basis to the spin-orbit basis:
      Call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,MM,1)
      Call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,ML,1)
      Call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,MS,1)
      Do Ist=1,nstate
         Mult=Multiplicity(Ist)
         Do I=-(Mult-Ipar)/2,(Mult-Ipar)/2
            If( (Ipar==0) .AND. (I==0) ) Go To 301
            Do J=-(Mult-Ipar)/2,(Mult-Ipar)/2
               If( (Ipar==0) .AND. (J==0) ) Go To 302
               Do l=1,3
                  i1=Ibas(Ist,I)
                  j1=Ibas(Ist,J)
                  MM(l,i1,j1)=-Spin(l,Mult,I,J)*g_e
                  MS(l,i1,j1)= Spin(l,Mult,I,J)
               End Do ! l
  302          Continue
            End Do ! J
  301       Continue
         End Do ! I
      End Do ! Ist


      Do Ist=1,nstate
         MultI=Multiplicity(Ist)
         Do Jst=1,nstate
            MultJ=Multiplicity(Jst)
            If(MultI==MultJ) Then
               Do I=-(MultI-Ipar)/2,(MultI-Ipar)/2
               If( (Ipar==0) .AND. (I.eq.0) ) Go To 303
                  Do l=1,3
                     i1=Ibas(Ist,I)
                     j1=Ibas(Jst,I)
                     MM(l,i1,j1)=MM(l,i1,j1)
     &                              -CMPLX(0.0_wp,Angmom(l,Ist,Jst),wp)
                     ML(l,i1,j1)=ML(l,i1,j1)
     &                              +CMPLX(0.0_wp,Angmom(l,Ist,Jst),wp)
                     DM(l,i1,j1)=DM(l,i1,j1)
     &                              +CMPLX(eDmom(l,Ist,Jst),0.0_wp,wp)
                  End Do   ! l
  303             Continue
               End Do   ! I
            End If
         End Do   ! Jst
      End Do   ! Ist

      ! calculate the matrix elements of the spin and magnetic moment
      ! in the spin-orbit basis:
      Call mma_allocate(tmp,nss,nss,'tmp')
      Do L=1,3
         TMP=(0.0_wp,0.0_wp)
         ! spin moment
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             MS(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             MS(L,:,:), nss )
         ! orbital moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             ML(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             ML(L,:,:), nss )
         ! magnetic moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             MM(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             MM(L,:,:), nss )

         If(found_EDMOM) Then
         ! electric dipole moment
         TMP=(0.0_wp,0.0_wp)
         CALL ZGEMM_('C', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                     U, nss,
     &             DM(L,:,:), nss,           (0.0_wp,0.0_wp),
     &                   TMP, nss )
         CALL ZGEMM_('N', 'N', nss, nss, nss, (1.0_wp,0.0_wp),
     &                   TMP, nss,
     &                     U, nss,           (0.0_wp,0.0_wp),
     &             DM(L,:,:), nss )
         End If
      End Do !L
      Call mma_deallocate(tmp)


      Return
      End Subroutine fetch_data_RunFile_all
