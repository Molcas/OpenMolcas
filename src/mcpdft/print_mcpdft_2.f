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
* Copyright (C) Giovanni Li Manni                                      *
************************************************************************
      SUBROUTINE print_MCPDFT_2(CASDFT_E,E_nuc,E_cor,E_cas,E_ot,
     & jroot,Ref_Ener)
******************************************************************
* Purpose:
* This routine is called from RASSCF when performing MC-PDFT jobs,
* for printing out some densities information and some functionals
* as computed in DRVNQ and children routines.
*
* Author:
* G. Li Manni (GLM)
* S Dong, 2018 (added print outs related to scaling)
******************************************************************
      use KSDFT_Info, only: CoefR, CoefX, Funcaa, Funcbb, Funccc
      use nq_Info, only: Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I
      Use hybridpdft, only: E_NoHyb
      use mspdft, only: mspdftmethod, do_rotate
      use printlevel, only: usual
      use mcpdft_output, only: iPrGlb
      use mcpdft_input, only: mcpdft_options

      Implicit Real*8 (A-H,O-Z)
      Real*8 CASDFT_E,E_nuc,E_cor,E_cas,E_ot
      Real*8 CASDFT_E_1,E_ot_1,Funcaa1,Funcbb1,Funccc1
      Dimension Ref_Ener(*)
      integer jroot, left
      character*6 Fmt2
      character*120 Line
#include "WrkSpc.fh"

      if (iPrGlb >= usual) then

      left=6
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
      Line=''
      IF(Do_Rotate) Then
      Write(Line(left-2:),'(2A,1X,I4.4)') MSpdftMethod,
     & ' INTERMEDIATE STATE',jroot
      ELSE
      Write(Line(left-2:),'(A,1X,I4.4)') 'MC-PDFT RESULTS, STATE',jroot
      ENDIF
      write(6,*)
      Call CollapseOutput(1,Line)
      Write(6,Fmt2//'A)') repeat('-',len_trim(Line)-3)
      write(6,*)

      write(6,'(6X,A,40X,F18.8)') 'MCSCF reference energy',
     &                           Ref_Ener(jroot)
      write(6,*)
      write(6,'(6X,A,45X,F10.3)') 'Integrated total density:',Dens_I
      write(6,'(6X,A,12X,F10.3)') 'Integrated alpha density '//
     &           'before functional transformation:', Dens_a1
      write(6,'(6X,A,12X,F10.3)') 'Integrated  beta density '//
     &           'before functional transformation:', Dens_b1
      write(6,'(6X,A,12X,F10.3)') 'Integrated alpha density '//
     &           ' after functional transformation:', Dens_a2
      write(6,'(6X,A,12X,F10.3)') 'Integrated  beta density '//
     &           ' after functional transformation:', Dens_b2
      write(6,'(6X,A)') 'NOTE:'
      write(6,'(6X,2A)') 'Densities after transformation are ',
     &           'intermediate quantities'
      write(6,'(6X,2A)') 'and should not be interpreted as ',
     &           'real spin densities'
      write(6,*)
      write(6,'(6X,A,32X,F18.6)') 'Exchange energy scaling factor',
     &          CoefX
      write(6,'(6X,A,29X,F18.6)') 'Correlation energy scaling factor',
     &          CoefR
      write(6,*)
      write(6,'(6X,A,30X,F18.6)') 'Integrated alpha exchange energy',
     &          Funcaa
      write(6,'(6X,A,30X,F18.6)') 'Integrated beta  exchange energy',
     &          Funcbb
      write(6,'(6X,A,30X,F18.6)') 'Integrated  correlation   energy',
     &          Funccc
      write(6,*)

      write(6,'(6X,A,38X,F18.8)') 'Nuclear Repulsion energy',E_nuc
      write(6,'(6X,A,51X,F18.8)') 'Core energy',E_cor
      write(6,'(6X,A,36X,F18.8)') 'CASSCF contribution energy',E_cas
      write(6,'(6X,A,49X,F18.8)') 'On-top energy',E_ot
      write(6,*)
      IF(mcpdft_options%otfnal%is_hybrid()) Then
       write(6,'(6X,A)') 'Information for hybrid PDFT:'
       write(6,'(6X,A,37X,F6.2)')
     & 'Wave function percentage (Lambda*100)',
     &  mcpdft_options%otfnal%lambda*1.0d2
       write(6,'(6X,A,42X,F18.8)')
     & 'Wave function energy',
     &  mcpdft_options%otfnal%lambda*Ref_Ener(jRoot)
       write(6,'(6X,A,51X,F18.8)')
     & 'PDFT energy',(1-mcpdft_options%otfnal%lambda)*E_NoHyb
       write(6,*)
      END IF

      IF(Do_Rotate) Then
      write(6,'(6X,A,2X,I4,13X,F18.8)')
     &'Total MC-PDFT energy for intermediate state', jroot,CASDFT_E
      ELSE
      call PrintResult(6,Fmt2//'A,2X,I4,A,26X,F18.8)',
     & 'Total MC-PDFT energy for state',jroot,'',[CASDFT_E],1)
      END IF
      if ((CoefX*CoefR.ne.0.0).and.(CoefX.ne.1.0.or.CoefR.ne.1.0)) Then
         Funcaa1 = Funcaa/CoefX
         Funcbb1 = Funcbb/CoefX
         Funccc1 = Funccc/CoefR
         E_ot_1 = E_ot-Funcaa-Funcbb-Funccc+Funcaa1+Funcbb1+Funccc1
         CASDFT_E_1 = CASDFT_E-E_ot+E_ot_1
         write(6,*)
         write(6,*)
         write(6,'(6X,A,19X,F18.6)') 'Integrated alpha exchange '//
     &          'energy (unscaled)',
     &          Funcaa1
         write(6,'(6X,A,19X,F18.6)') 'Integrated beta  exchange '//
     &          'energy (unscaled)',
     &          Funcbb1
         write(6,'(6X,A,19X,F18.6)') 'Integrated  correlation   '//
     &          'energy (unscaled)',
     &          Funccc1
!         write(6,*)
         write(6,'(6X,A,38X,F18.8)') 'On-top energy (unscaled)',E_ot_1
!         write(6,*)
         write(6,'(6X,A,31X,F18.8)') 'Total MC-PDFT energy '//
     &         '(unscaled)',
     &         CASDFT_E_1
      end if

      Call CollapseOutput(0,Line)
      write(6,*)
      end if

      Call Add_Info('dens_tt',[Dens_I],1,6)
      Call Add_Info('dens_a1',[Dens_a1],1,6)
      Call Add_Info('dens_b1',[Dens_b1],1,6)
      Call Add_Info('dens_a2',[Dens_a2],1,6)
      Call Add_Info('dens_b2',[Dens_b2],1,6)
      Call Add_Info('exch_f',[CoefX],1,6)
      Call Add_Info('corr_f',[CoefR],1,6)
      Call Add_Info('excha_a',[Funcaa],1,6)
      Call Add_Info('excha_b',[Funcbb],1,6)
      Call Add_Info('corr_e', [Funccc],1,6)
      Call Add_Info('CASDFTE',[CASDFT_E],1,8)

      RETURN
      END
