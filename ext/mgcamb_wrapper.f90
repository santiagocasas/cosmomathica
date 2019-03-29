module cambwrapper

    use iso_c_binding

    use CAMB
    use MGCAMB
    use LambdaGeneral
    use Lensing
    use AMLUtils
    use Transfer
    use constants
    use Bispectrum
    use CAMBmain
    use ModelParams
    use Precision
    use ModelData
    use GaugeInterface
    use InitialPower
    use NonLinear
    use Reionization
    use Recombination
    
    implicit none

contains
    
    function int2bool(i)
        implicit none
        integer, intent(in) :: i
        logical :: int2bool
        
        if (i == 0) then
            int2bool = .false.
        else
            int2bool = .true.
        endif

    end function int2bool


    subroutine add1Darray(array, output, offset, intarray, intoffset)
        implicit none
        integer, intent(inout) :: offset, intoffset
        double precision, intent(in) :: array(:)
        double precision, intent(out) :: output(:)
        integer, intent(out) :: intarray(:)
        integer i

        if (offset + size(array) > size(output)) stop "camb_wrapper: array too small"

        do i=1,size(array,1)
            output(offset+i) = array(i)
        end do

        intarray(intoffset+1) = 1
        intarray(intoffset+2) = size(array)
        offset = offset + size(array)
        intoffset = intoffset + 2

    end subroutine add1Darray


    subroutine add2Darray(array, output, offset, intarray, intoffset)
        implicit none
        integer, intent(inout) :: offset, intoffset
        double precision, intent(in) :: array(:,:)
        double precision, intent(out) :: output(:)
        integer, intent(out) :: intarray(:)
        integer i,j

        if (offset + size(array) > size(output)) stop "camb_wrapper: array too small"

        do i=1,size(array,1)
            do j=1,size(array,2)
                output(offset+j+(i-1)*size(array,2)) = array(i,j)
            end do
        end do

        intarray(intoffset+1) = 2
        intarray(intoffset+2) = size(array,1)
        intarray(intoffset+3) = size(array,2)
        offset = offset + size(array)
        intoffset = intoffset + 3

    end subroutine add2Darray


    subroutine add3Darray(array, output, offset, intarray, intoffset)
        implicit none
        integer, intent(inout) :: offset, intoffset
        double precision, intent(in) :: array(:,:,:)
        double precision, intent(out) :: output(:)
        integer, intent(out) :: intarray(:)
        integer i,j,k

        if (offset + size(array) > size(output)) stop "camb_wrapper: array too small"

        do i=1,size(array,1)
            do j=1,size(array,2)
                do k=1,size(array,3)
                    output(offset+k+(j-1)*size(array,3)+(i-1)*size(array,2)*size(array,3)) = array(i,j,k)
                end do
            end do
        end do

        intarray(intoffset+1) = 3
        intarray(intoffset+2) = size(array,1)
        intarray(intoffset+3) = size(array,2)
        intarray(intoffset+4) = size(array,3)
        offset = offset + size(array)
        intoffset = intoffset + 4

    end subroutine add3Darray


    subroutine runcamb(floats, floats_len, ints, ints_len, floats_out, floats_out_len, ints_out, ints_out_len) bind(c)
        use LambdaGeneral
        use NonLinear 
        use Transfer
        implicit none

        integer(c_int), intent(in)    :: floats_len, ints_len
        integer(c_int), intent(inout) :: ints_out_len, floats_out_len
        real(c_double), intent(in)    :: floats(floats_len)
        integer(c_int), intent(in)    :: ints(ints_len)
        real(c_double), intent(out)   :: floats_out(floats_out_len)
        integer(c_int), intent(out)   :: ints_out(ints_out_len)
        double precision :: omegak
        integer error, fi, ii, eigenstates, fitemp, i, float_offset, int_offset
        integer testkper

        Type(CAMBparams) P
        Type(MatterPowerData) PK_data
        Type (CAMBdata)  :: OutData
        Type (MGCAMB_parameter_cache) :: mgcamb_par_cache
        ! We have two arrays containing all the parameters: one full of doubles,
        ! one full of ints. We need to fill the structure P with it. See
        ! camb/inidriver.f90 for how it works. We avoid passing structures
        ! between C and Fortran because of compatibility issues, so we pass flat
        ! arrays.


        ! call CAMB_SetDefParams(P)
        ! Compare the following to the function CAMB_SetDefParams.
        ! Obviously the order of the following matters a lot - don't change it
        ! unless you really know what you're doing.

        fi = 1
        ii = 1


        


        !### Standard CAMB parameters
        P%omegac  = floats(fi); fi=fi+1
        P%omegab  = floats(fi); fi=fi+1
        P%H0      = floats(fi); fi=fi+1
        P%omegan  = floats(fi); fi=fi+1
        omegak    = floats(fi); fi=fi+1
        P%omegav  = 1 - omegak - P%omegab-P%omegac - P%omegan
       
        !use_tabulated_w = .false. 
        w0DE     = floats(fi); fi=fi+1
        waDE    = floats(fi); fi=fi+1

        P%TCMB    = floats(fi); fi=fi+1
        P%YHe     = floats(fi); fi=fi+1

        !## Fill mgcamb_par_cache      
         mgcamb_par_cache%omegab = P%omegab
         mgcamb_par_cache%omegac = P%omegac
         mgcamb_par_cache%omegav = P%omegav
         mgcamb_par_cache%h0     = P%H0
         mgcamb_par_cache%h0_Mpc = P%H0 * (1.d3/c)
         ! mgcamb_par_cache%output_root = outroot
 


        !###### Floats, parameters for above parametrizations
        !# Choose at which time to turn on MG
        GRtrans = floats(fi); fi=fi+1
        !#BZ parameters:
        B1   = floats(fi); fi=fi+1
        lambda1_2 = floats(fi); fi=fi+1
        B2 = floats(fi); fi=fi+1
        lambda2_2 = floats(fi); fi=fi+1
        ss = floats(fi); fi=fi+1
        !#Planck parameters
        E11 = floats(fi); fi=fi+1
        E22 =  floats(fi); fi=fi+1
       !###### Part 3.1.2. - mu, Sigma functions
        !# musigma_par = 1 : DES parametrization
        mu0 =  floats(fi); fi=fi+1
        sigma0 =  floats(fi); fi=fi+1
       !###### Part 3.1.3. - Q,R functions
        !#Bean parameters :
        !#(Q,R)
        MGQfix= floats(fi); fi=fi+1
        MGRfix= floats(fi); fi=fi+1       
        !#(Q0,R0,s)
        Qnot=floats(fi); fi=fi+1
        Rnot=floats(fi); fi=fi+1
        sss=floats(fi); fi=fi+1
        !##### Part 3.2.1 - Linder Gamma
        ! Linder's gamma :
        Linder_gamma =floats(fi); fi=fi+1
        !##### Part 3.3.1 - QSA f(R) model
        !B0 =floats(fi); fi=fi+1   !! not declared in mgcamb.f90
        !##### Part 3.3.2 - QSA Symmetron model
        beta_star = floats(fi); fi=fi+1
        a_star = floats(fi); fi=fi+1
        xi_star = floats(fi); fi=fi+1
        !##### Part 3.3.3 - QSA Dilaton model
        beta0 = floats(fi); fi=fi+1
        xi0 = floats(fi); fi=fi+1
        DilS = floats(fi); fi=fi+1
        DilR = floats(fi); fi=fi+1
        !!A2 = floats(fi); fi=fi+1!! not declared in mgcamb.f90
        !##### Part 3.3.4 - QSA Hu-Sawicki f(R)
        F_R0 = floats(fi); fi=fi+1 
        FRn = floats(fi); fi=fi+1
        !< MGCAMB MOD END


        !> MGCAMB MOD START
         !###### Part 1. Choose the Modified Growth flag
         !
         !# MG_flag = 0 :  default GR
         !# MG_flag = 1 :  pure MG models
         !# MG_flag = 2 :  alternative MG models
         !# MG_flag = 3 :  QSA models
        MG_flag=ints(ii); ii=ii+1
        !###### Part 2.1 - Pure MG models
        !# pure_MG_flag = 1 : mu, gamma parametrization
        !# pure_MG_flag = 2 : mu, sigma parametrization
        !# pure_MG_flag = 3 : Q, R  parametrization
        pure_MG_flag = ints(ii); ii=ii+1
        !###### Part 2.2 - Alternative MG models
        !# alt_MG_flag = 1 : Linder Gamma parametrization ( introduced in arXiv:0507263 )
        alt_MG_flag = ints(ii); ii=ii+1
         !###### Part 2.3 - QSA models
         !# QSA_flag = 1 : f(R)
         !# QSA_flag = 2 : Symmetron
         !# QSA_flag = 3 : Dilaton
         !# QSA_flag = 4 : Hu-Sawicki f(R)
        QSA_flag= ints(ii); ii=ii+1
         !###### Part 3.1.1. - mu, gamma functions
         !# mugamma_par = 1 : BZ parametrization(arXiv:0809.3791 )
         !# mugamma_par = 2 : Planck parametrization
        mugamma_par = ints(ii); ii=ii+1
         !###### Part 3.1.2. - mu, Sigma functions
         !# musigma_par = 1 : DES parametrization
        musigma_par = ints(ii); ii=ii+1
          !###### Part 3.1.3. - Q,R functions
          !# QR_par = 1 : (Q,R)(arXiv:1002.4197 )
          !# QR_par = 2 : (Q0,R0,s)(arXiv:1002.4197 )
        QR_par = ints(ii); ii=ii+1
          !##DE_model = 0 : LCDM
          !##DE_model = 1 : wCDM                
          !##DE_model = 2 : (w0,wa)CDM          
          !##DE_model = 3 : user defined                
        DE_model = ints(ii); ii=ii+1
     !> MGCAMB MOD START: reading models and params
      !     call MGCAMB_read_model_params( mgcamb_par_cache )
     !< MGCAMB MOD END

        !! ### Neutrino Parameters!!
        P%Num_Nu_massless     = floats(fi); fi=fi+1
        P%Num_Nu_massive      = ints(ii); ii=ii+1
        P%Nu_mass_eigenstates = ints(ii); ii=ii+1
        eigenstates = P%Nu_mass_eigenstates
        P%share_delta_neff = .true.  
        !temporary fix for standard neutrinos
        !P%Nu_mass_degeneracies(1:eigenstates) = floats(fi:fi+eigenstates)
        fi = fi+eigenstates
        P%Nu_mass_fractions(1:eigenstates) = floats(fi:fi+eigenstates)
        fi = fi+eigenstates

        P%Scalar_initial_condition = ints(ii); ii=ii+1
        P%NonLinear                = ints(ii); ii=ii+1
        halofit_version            = ints(ii); ii=ii+1
        ! call SetDefPowerParams(P%InitPower)
        ! Compare with SetDefPowerParams
        ! These are arrays with length nn
         P%InitPower%nn     = ints(ii); ii=ii+1 !number of initial power spectra
         fitemp = ints(ii-1)
         P%InitPower%an(1:fitemp)     = floats(fi:fi+fitemp); fi=fi+fitemp !scalar spectral index
         P%InitPower%n_run(1:fitemp)  = floats(fi:fi+fitemp); fi=fi+fitemp !running of scalar spectral index
         P%InitPower%ant(1:fitemp)    = floats(fi:fi+fitemp); fi=fi+fitemp !tensor spectra index
         P%InitPower%rat = 1d0
         P%InitPower%rat(1:fitemp)    = floats(fi:fi+fitemp); fi=fi+fitemp
         P%InitPower%ScalarPowerAmp(1:fitemp) = floats(fi:fi+fitemp); fi=fi+fitemp
         P%InitPower%k_0_scalar = floats(fi); fi=fi+1
         P%InitPower%k_0_tensor = floats(fi); fi=fi+1

        ! call Recombination_SetDefParams(P%Recomb)
             P%Recomb%RECFAST_fudge = 1.14d0
             P%Recomb%RECFAST_fudge_He =.86d0
             P%Recomb%RECFAST_Heswitch = 6
             P%Recomb%RECFAST_Hswitch  = .true.
        ! call Reionization_SetDefParams(P%Reion)
         P%Reion%Reionization      = int2bool(ints(ii)); ii=ii+1
         P%Reion%use_optical_depth = int2bool(ints(ii)); ii=ii+1
         P%Reion%optical_depth     = floats(fi); fi=fi+1
         P%Reion%redshift          = floats(fi); fi=fi+1
         P%Reion%fraction          = floats(fi); fi=fi+1
         P%Reion%delta_redshift    = floats(fi); fi=fi+1
        
         AccuracyBoost             = floats(fi); fi=fi+1;
         !TODO bispectrum?

        P%Transfer%high_precision   = int2bool(ints(ii)); ii=ii+1
        transfer_interp_matterpower = int2bool(ints(ii)); ii=ii+1

        P%Want_CMB     = int2bool(ints(ii)); ii=ii+1
        P%PK_WantTransfer = int2bool(ints(ii)); ii=ii+1
        P%WantTransfer = P%PK_WantTransfer 
        P%WantCls      = int2bool(ints(ii)); ii=ii+1
        P%WantScalars = int2bool(ints(ii)); ii=ii+1
        P%WantVectors = int2bool(ints(ii)); ii=ii+1
        P%WantTensors = int2bool(ints(ii)); ii=ii+1
        P%want_zstar  = int2bool(ints(ii)); ii=ii+1
        P%want_zdrag  = int2bool(ints(ii)); ii=ii+1

        P%OutputNormalization = ints(ii); ii=ii+1

        P%Max_l            = ints(ii); ii=ii+1
        P%Max_eta_k        = floats(fi); fi=fi+1
        P%Max_l_tensor     = ints(ii); ii=ii+1
        P%Max_eta_k_tensor = floats(fi); fi=fi+1
        P%Transfer%kmax          = floats(fi); fi=fi+1
        testkper  = ints(ii); ii=ii+1
        P%Transfer%k_per_logint  = testkper
        !P%Transfer%k_per_logint  = ints(ii); ii=ii+1
        P%Transfer%PK_num_redshifts = ints(ii); ii=ii+1
        ! This is an array with length num_redshifts
        P%Transfer%PK_redshifts     = floats(fi:fi+ints(ii-1)); fi=fi+ints(ii-1) 

        P%AccuratePolarization = int2bool(ints(ii)); ii=ii+1
        P%AccurateReionization = int2bool(ints(ii)); ii=ii+1
        P%AccurateBB           = int2bool(ints(ii)); ii=ii+1

        P%DoLensing = int2bool(ints(ii)); ii=ii+1
        P%OnlyTransfers   = int2bool(ints(ii)); ii=ii+1
        P%DerivedParameters = int2bool(ints(ii)); ii=ii+1
        P%MassiveNuMethod = ints(ii); ii=ii+1

        CP%DerivedParameters = .true.
        !TODO DoTensorNeutrinos? ThreadNum?

            write(42,*) "Wrong number of parameters: ", fi-1,ii-1
            write(42,*) "Expected: ", floats_len, ints_len
        
       if (ii-1 /= ints_len .or. fi-1 /= floats_len) then
            write(42,*) "Wrong number of parameters: ", fi-1,ii-1
            write(42,*) "Expected: ", floats_len, ints_len
            stop
        endif

        !call cmbmain
        error = 0
        call CAMB_GetTransfers(P, OutData, error)
        !call Transfer_SetForNonlinearLensing(P%Transfer) 
        !this routine resets the k_max_per_logint to 0
        ! cannot use NonLinear = "both" anymore  
        call Transfer_SortAndIndexRedshifts(P%Transfer)

        error = 0
        if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'
        call  CAMBParams_Set(P, error, .false.)
        if (error>0) write(*,*) "Error: ",error,trim(global_error_message)

        !call InitTransfer

        call CAMB_GetResults(P)
        if (global_error_flag/=0) then
            write (*,*) 'Error result '//trim(global_error_message)
            stop
        endif
        open(23, file='./testprints.txt', status='REPLACE', ACTION="READWRITE")

        write(23,*)  'File containing test prints'
        write(23,*)  'num_q_trans=',MT%num_q_trans
        write(23,*)  'kmax=',CP%Transfer%kmax
        write(23,*)  'passed P : kmax=',P%Transfer%kmax
        write(23,*)  'CP:  k_per_logint=',CP%Transfer%k_per_logint
        write(23,*)  'passed P : k_per_logint=',P%Transfer%k_per_logint
        write(23,*)  'really passed testkper: k_per_logint=',testkper
        close(23)

        ! the global variables MT, Cl_scalar, Cl_vector, Cl_tensor now contain
        ! meaningful data

        ! Return the results

        ints_out(1) = error
        int_offset = 1

!#define add3d(array) call add3darray(array, floats_out, float_offset, ints_out, int_offset)
!#define add2d(array) call add2darray(array, floats_out, float_offset, ints_out, int_offset)
!#define add1d(array) call add1darray(array, floats_out, float_offset, ints_out, int_offset)

        ! Derived parameters first
        float_offset = 0
        call add1darray( (/ ThermoDerivedParams( derived_Age ),&
            ThermoDerivedParams( derived_zstar ),&
            ThermoDerivedParams( derived_rstar ),&
            ThermoDerivedParams( derived_thetastar ),&
            ThermoDerivedParams( derived_zdrag ),&
            ThermoDerivedParams( derived_rdrag ),&
            ThermoDerivedParams( derived_kD ),&
            ThermoDerivedParams( derived_thetaD ),&
            ThermoDerivedParams( derived_zEQ ),&
            ThermoDerivedParams( derived_thetaEQ ) /),&
        floats_out, float_offset, ints_out, int_offset)

        if (P%WantCls) then 
           if(P%WantScalars) call add3darray(Cl_scalar, floats_out, float_offset, ints_out, int_offset)   !add1d(Cl_scalar)   !add3d
           if(P%WantVectors) call add3darray(Cl_vector, floats_out, float_offset, ints_out, int_offset)   !add3d(Cl_vector)
           if(P%WantTensors) call add3darray(Cl_tensor, floats_out, float_offset, ints_out, int_offset)   !add3d(Cl_tensor)
       endif 

        if (P%WantTransfer) then
            do i=1,P%InitPower%nn     
                call Transfer_GetMatterPowerData(MT, PK_data, i)
                call add1darray(PK_data%log_kh, floats_out, float_offset, ints_out, int_offset) !add1d(PK_data%log_kh)
                call add2darray(PK_data%matpower, floats_out, float_offset, ints_out, int_offset) !add2d(PK_data%matpower)
                call MatterPowerdata_MakeNonlinear(PK_data)
                call add2darray(PK_data%matpower, floats_out, float_offset, ints_out, int_offset) !add2d(PK_data%matpower)
            end do
            call add3darray(dble(MT%TransferData), floats_out, float_offset, ints_out, int_offset) !add3d(dble(MT%TransferData))
            call add1darray(PK_data%redshifts, floats_out, float_offset, ints_out, int_offset) !add1d(PK_data%redshifts)
            call Transfer_Get_sigma8(MT, 8d0)
            call add2darray(MT%sigma_8, floats_out, float_offset, ints_out, int_offset) !add2d(MT%sigma_8)
        endif


        ! Background: TODO: CAMB never fills these arrays. But the functions
        ! exist, so we have to do it by hand.
        if (associated(BackgroundOutputs%z_outputs)) then
            call add1darray(BackgroundOutputs%z_outputs, floats_out, float_offset, ints_out, int_offset) !add1d(BackgroundOutputs%z_outputs)
            call add1darray(BackgroundOutputs%H, floats_out, float_offset, ints_out, int_offset) !add1d(BackgroundOutputs%H)
            call add1darray(BackgroundOutputs%DA, floats_out, float_offset, ints_out, int_offset) !add1d(BackgroundOutputs%DA)
            call add1darray(BackgroundOutputs%rs_by_D_v, floats_out, float_offset, ints_out, int_offset) !add1d(BackgroundOutputs%rs_by_D_v)
        endif

        ints_out_len = int_offset
        floats_out_len = float_offset
        
        call CAMB_cleanup
    end subroutine runcamb


end module cambwrapper
