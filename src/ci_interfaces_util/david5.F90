Subroutine Mk_H_Psi(SGS,EXS,CIS,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb, &
                    nTU,TU,nTUVX,TUVX)
use Lucia_Interface, only: Lucia_Util
use lucia_data, only: Sigma_on_disk
use citrans, only: citrans_csf2sd, citrans_sd2csf, citrans_sort
use sguga, only: SGStruct, EXStruct, CIStruct
use rasscf_global, only: DoFaro
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use Constants, only: Zero
use faroald, only: my_norb, sigma_update, htu, gtuvx
use definitions, only: wp, iwp, u6
Implicit None

type(SGStruct), intent(in) :: SGS
type(EXStruct), intent(in) :: EXS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(in):: STSYM, nCSF
real(kind=wp), intent(in) :: CI_Vec(nCSF)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
integer(kind=iwp), intent(in) :: ntemp,ndeta,ndetb
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
integer(kind=iwp), intent(in):: nTU, nTUVX
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

integer(kind=iwp) :: itu, ituvx, it, iu, iv, ixmax, ix, iprlev
real(kind=wp), pointer:: Faroald_PSI(:,:), Faroald_SGM(:,:)

if (DOFARO) then

  htu(:,:) = Zero
  gtuvx(:,:,:,:) = Zero
  itu = 0
  ituvx = 0
  do it=1,my_norb
    do iu=1,it
      itu = itu+1
      htu(iu,it) = TU(itu)
      htu(it,iu) = TU(itu)
      do iv=1,it
        ixmax = iv
        if (it == iv) ixmax = iu
        do ix=1,ixmax
          ituvx = ituvx+1
          GTUVX(IT,IU,IV,IX) = TUVX(ITUVX)
          GTUVX(IU,IT,IV,IX) = TUVX(ITUVX)
          GTUVX(IT,IU,IX,IV) = TUVX(ITUVX)
          GTUVX(IU,IT,IX,IV) = TUVX(ITUVX)
          GTUVX(IV,IX,IT,IU) = TUVX(ITUVX)
          GTUVX(IX,IV,IT,IU) = TUVX(ITUVX)
          GTUVX(IV,IX,IU,IT) = TUVX(ITUVX)
          GTUVX(IX,IV,IU,IT) = TUVX(ITUVX)
        end do
      end do
    end do
  end do

  Faroald_Psi(1:nDetA,1:nDetB) => ctemp(:)
  Faroald_SGM(1:nDetA,1:nDetB) => sigtemp(:)

  call SG_REORD(SGS,EXS,STSYM,0,CIS%nCSF(STSYM),CI_Vec,ctemp)
  call CITRANS_SORT('C',ctemp,Sigma_Vec)
  Faroald_PSI(:,:) = Zero
  call CITRANS_CSF2SD(Sigma_Vec,Faroald_PSI)
  Faroald_SGM(:,:) = Zero
  call SIGMA_UPDATE(HTU,GTUVX,Faroald_SGM,Faroald_PSI)
  call CITRANS_SD2CSF(Faroald_SGM,Sigma_Vec)
  call CITRANS_SORT('O',Sigma_Vec,ctemp)
  call SG_Reord(SGS,EXS,STSYM,1,CIS%nCSF(STSYM),ctemp,Sigma_Vec)

  Faroald_Psi => Null()
  Faroald_SGM => Null()

else

  ! Convert the CI-vector from CSF to Det. basis
  ! sigtemp is scratch, converted vector is stored in ctemp.
  ! Note that ctemp is of the size of the nDet. basis

  ctemp(1:nCSF) = CI_Vec(1:nCSF)
  sigtemp(:) = Zero
  call csdtvc(ctemp,sigtemp,1,STSym,1)

  ! Calling Lucia to determine the sigma vector
  call Lucia_Util('Sigma',                   &
                  CI_Vector=ctemp(:),        &
                  Sigma_Vector=sigtemp(:),   &
                  nTU=Size(TU),TU=TU,        &
                  nTUVX=Size(TUVX),TUVX=TUVX)

  ! Set mark so densi_master knows that the Sigma-vector exists on disk.
  Sigma_on_disk = .true.
  ! Convert the Sigma vector from Det. to CSF basis. Converted vector is
  ! stored in Sigma_vec.
  call CSDTVC(Sigma_Vec,sigtemp,2,stSym,1)

end if

End Subroutine Mk_H_Psi
