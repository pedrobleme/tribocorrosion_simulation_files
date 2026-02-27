!> @brief Module to centralize default Abaqus parameters.
!!
!! This module replaces the need to use 'INCLUDE vaba_param.inc' and works correctly with IMPLICIT NONE.
!--------------------------------------------------------------------------------------------------------------
      MODULE AbaqusParameters
      IMPLICIT NONE

    ! --- Job control parameters (lOp) ---
      INTEGER, PARAMETER :: J_INT_START_ANALYSIS    = 0
      INTEGER, PARAMETER :: J_INT_START_STEP        = 1
      INTEGER, PARAMETER :: J_INT_SETUP_INCREMENT   = 2
      INTEGER, PARAMETER :: J_INT_START_INCREMENT   = 3
      INTEGER, PARAMETER :: J_INT_END_INCREMENT     = 4
      INTEGER, PARAMETER :: J_INT_END_STEP          = 5
      INTEGER, PARAMETER :: J_INT_END_ANALYSIS      = 6

    ! --- Indices for the i_Array (Job information) ---
      INTEGER, PARAMETER :: I_INT_NTOTALNODES       = 1
      INTEGER, PARAMETER :: I_INT_NTOTALELEMENTS    = 2
      INTEGER, PARAMETER :: I_INT_KSTEP             = 3
      INTEGER, PARAMETER :: I_INT_KINC              = 4
      INTEGER, PARAMETER :: I_INT_ISTATUS           = 5
      INTEGER, PARAMETER :: I_INT_LWRITERESTART     = 6

    ! --- Indices for the r_Array (Time information) ---
      INTEGER, PARAMETER :: I_FLT_TOTALTIME         = 1
      INTEGER, PARAMETER :: I_FLT_STEPTIME          = 2
      INTEGER, PARAMETER :: I_FLT_DTIME             = 3

      END MODULE AbaqusParameters
!--------------------------------------------------------------------------------------------------------------
!> @brief Declares global variables, parameters and shared arrays.
!!
!! This module centralizes all variables that need to be accessed by different subroutines of the simulation.
!--------------------------------------------------------------------------------------------------------------
      MODULE SharedVariables
      IMPLICIT NONE

    ! --- Numerical Precision ---
      INTEGER, PARAMETER :: dp = KIND(1.0d0)

    ! --- Simulation Parameters ---
      INTEGER             :: n_elems             = 0    ! Total number of elements in the mesh.
      INTEGER, PARAMETER  :: max_influence_elems = 500  ! Maximum number of elements in the influence radius.
      INTEGER, PARAMETER :: max_face_neighbors = 6

    ! --- Physical and Model Parameters ---
      REAL(KIND=dp), PARAMETER :: lr    = 15.7      ! [mm] Intrinsic length for non-local model.
      REAL(KIND=dp), PARAMETER :: nt    = 10950.0  ! [days] Total corrosion time.
      REAL(KIND=dp), PARAMETER :: sigth = 121      ! [MPa] Material yield stress threshold.

    ! --- Indices for the elem_properties Matrix ---
    ! This matrix stores all properties read from the file and topology information per element.
      INTEGER, PARAMETER :: PROP_IDX_ELEM_ID           = 1   ! Column 1: Element number (ID)
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_1        = 2   ! Column 2: Neighbor ID on face 1
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_2        = 3   ! Column 3: Neighbor ID on face 2
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_3        = 4   ! Column 4: Neighbor ID on face 3
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_4        = 5   ! Column 5: Neighbor ID on face 4
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_5        = 6   ! Column 6: Neighbor ID on face 5
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_6        = 7   ! Column 7: Neighbor ID on face 6
      INTEGER, PARAMETER :: PROP_IDX_SURFACE_FLAG      = 8   ! Column 8: Flag indicating whether element is a surface (1.0) or not (0.0)
      INTEGER, PARAMETER :: PROP_IDX_PITTING           = 9   ! Column 9: Normalized pitting value
      INTEGER, PARAMETER :: PROP_IDX_VOLUME            = 10  ! Column 10: Element volume
      INTEGER, PARAMETER :: PROP_IDX_COORD_X           = 11  ! Column 11: Centroid X coordinate
      INTEGER, PARAMETER :: PROP_IDX_COORD_Y           = 12  ! Column 12: Centroid Y coordinate
      INTEGER, PARAMETER :: PROP_IDX_COORD_Z           = 13  ! Column 13: Centroid Z coordinate
      INTEGER, PARAMETER :: PROP_IDX_PITTING_ABSOLUTE  = 14  ! Column 14: Absolute pitting value

    !-----------------------------------------------------------------------
    ! --- Flags and Analysis Control States ---
    !-----------------------------------------------------------------------
      LOGICAL :: failure_occurred_this_increment = .FALSE.

    ! --- Corrosion Start Marker ---
      INTEGER, PARAMETER :: flag_corrosion_started = 10

    ! --- Element Deletion Control ---
      INTEGER, PARAMETER :: delete_element_flag = 0      ! Flag to trigger element deletion via stateNew(k,3).

    ! --- Element corrosion status (element_corrosion_status) ---
      INTEGER, PARAMETER :: flag_corrosion_active   = 0    ! Element is subject to corrosion.
      INTEGER, PARAMETER :: flag_corrosion_inactive = -10  ! Element is no longer corroding.

    ! --- Strain-based deletion status (strain_deletion_status) ---
      INTEGER, PARAMETER :: flag_deletion_by_strain_enabled  = 0    ! Strain deletion enabled.
      INTEGER, PARAMETER :: flag_deletion_by_strain_disabled = -20  ! Strain deletion disabled.

    ! --- Property Update Lock (property_update_lock) ---
      INTEGER, PARAMETER :: key_prop_update_unlocked = 0    ! Element properties may be updated.
      INTEGER, PARAMETER :: key_prop_update_locked   = -5   ! Element properties are locked.

    ! --- Trigger for General Property Transfer ---
      INTEGER, PARAMETER :: flag_general_property_transfer        = -8
      INTEGER, PARAMETER :: flag_general_property_transfer_locked = 0

    ! --- Flag for Global Corrosion Stop (global_corrosion_stop_flag) ---
      INTEGER, PARAMETER :: flag_global_stop = -15

    ! --- Execution control flags ---
      INTEGER :: prop_update_trigger         ! Trigger for the transfer_properties routine.
      INTEGER :: global_corrosion_stop_flag  ! Controls global stop of the corrosion process.

    ! --- Material and Simulation Properties ---
      REAL(KIND=dp)                                 :: mCO, timetot, lmCO
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: elem_properties

    ! --- Shared Matrices and Vectors ---
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: element_corrosion_status
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: strain_deletion_status
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: property_update_lock
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: counter
      INTEGER,       DIMENSION(:,:), ALLOCATABLE    :: influ
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: dist
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: temp_pitting_write
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: temp_surface_flag_write
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: nonlocal_property
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: temp_nonlocal_write
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: damage_current
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: von_mises_stress
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: non_local_stress

    ! --- Variables for Tribocorrosion ---
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: tribological_damage_accumulated
      REAL(KIND=dp)                               :: accumulated_pure_tribo
      REAL(KIND=dp)                               :: enhanced_tribo_damage_inc
      REAL(KIND=dp)                               :: combined_damage, total_damage

      END MODULE SharedVariables
!-------------------------------------------------------------------------------------------------------------------------------------
!> @brief Main VUMAT entry point called by Abaqus.
!!
! This subroutine acts as a wrapper. Its role is to unpack arguments from the 'jblock' vector
! (such as element number and integration point) and pass them in an organized way to the main worker subroutine `vumatXtrArg`.
!-------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE vumat (
! Read only -
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )
!
      INCLUDE 'vaba_param.inc'
!     
!      USE subs
      DIMENSION jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)

      character*80 cmname
      PARAMETER (i_umt_nblock = 1,
     *     i_umt_npt    = 2,
     *     i_umt_layer  = 3,
     *     i_umt_kspt   = 4,
     *     i_umt_noel   = 5 )

      CALL  vumatXtrArg ( jblock(i_umt_nblock),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(i_umt_noel), jblock(i_umt_npt),
     *     jblock(i_umt_layer), jblock(i_umt_kspt))

      RETURN
      END SUBROUTINE vumat
!-------------------------------------------------------------------------------------------------------------------------------------
!> @brief Main implementation of the constitutive model logic.
!!
!! This subroutine contains all calculations for the damage model (corrosion, wear, tribocorrosion)
!! and plasticity. It receives a clearer argument list than `vumat` to ease readability and maintenance.
!-------------------------------------------------------------------------------------------------------------------------------------  
      SUBROUTINE vumatXtrArg (
  ! Input arguments (Read only)
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relS  pinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
  ! Output arguments (Write only)
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
  ! Extra arguments (Read only)
     *     nElement, nMatPoint, nLayer, nSecPoint )
     
      USE SharedVariables
      INCLUDE 'vaba_param.inc'
  !-----------------------------------------------------------------------
  ! --- 1. DECLARATION OF VARIABLES AND ARGUMENTS ---
  !-----------------------------------------------------------------------
    ! --- Subroutine Arguments ---
      DIMENSION props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
    ! Documentation of extra arguments:
    !  nElement: Array of internal element numbers
    !  nMatPoint: Integration point number
    !  nLayer   : Layer number for composite shells and layered solids
    !  nSecPoint: Section point number within the current layer
      DIMENSION nElement(nblock), nLayer(nblock),nMatPoint(nblock)

          character*80 cmname
      ! --- Local Calculation Variables ---
        INTEGER :: k, nvalue, OK

      ! Working variables for stresses and strains
        REAL(KIND=dp) :: s11, s22, s33, s12, s13, s23, trace
        REAL(KIND=dp) :: smean, sigdif, vmises, factor
        REAL(KIND=dp) :: peeqOld_k, deqps, eqps, plasticWorkInc
        REAL(KIND=dp) :: yieldOld, yieldNew, hard
        REAL(kind=dp) :: e, xnu, twomu, alamda, thremu
        REAL(kind=dp) ::  masslimit, deltat ! Verify necessity of these variables

      ! Working variables for energy
        REAL(KIND=dp) :: stressPower

      ! Working variables for corrosion and damage
        REAL(KIND=dp) :: tnt, ku, corrosion_damage
        REAL(KIND=dp) :: dCOld_k, dSCld_k
        REAL(KIND=dp) :: sigmaJK, epsilonJK, epsilonf
        REAL(KIND=dp) :: initial_volume
      
      ! Working variables for tribology
        REAL(KIND=dp) :: calculated_slip_inc, contact_area
        REAL(KIND=dp) :: force_of_friction, energy_dissipated_increment
        REAL(KIND=dp) :: wear_volume_increment, tribological_damage
        REAL(KIND=dp) :: tribological_damage_inc

      ! Local arrays
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: damage_new
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: relative_slip_increment
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: contact_pressure
        REAL(KIND=dp), DIMENSION(n_elems)        :: elem_length
        REAL(KIND=dp), DIMENSION(n_elems)        :: corroded_volume 
                

  !-----------------------------------------------------------------------
  ! --- 2. ASSIGNMENT OF PHYSICAL AND MODEL PARAMETERS ---
  !-----------------------------------------------------------------------
      ! --- Mathematical and numerical constants ---
      REAL(KIND=dp), PARAMETER :: zero         = 0.0_dp
      REAL(KIND=dp), PARAMETER :: one          = 1.0_dp
      REAL(KIND=dp), PARAMETER :: two          = 2.0_dp
      REAL(KIND=dp), PARAMETER :: half         = 0.5_dp
      REAL(KIND=dp), PARAMETER :: third        = 1.0_dp / 3.0_dp
      REAL(KIND=dp), PARAMETER :: one_point_five = 1.5_dp

      ! --- General Damage Model Parameters ---
      REAL(KIND=dp), PARAMETER :: dmax = 0.999_dp         ! Maximum damage before element deletion [dimensionless]
      REAL(KIND=dp), PARAMETER :: beta = 0.8_dp           ! Damage propagation factor for neighbors [dimensionless]

      ! --- Corrosion Model Parameters ---
      REAL(KIND=dp), PARAMETER :: corrosion_rate       = 0.471_dp ! Corrosion rate [mm/year]
      REAL(KIND=dp), PARAMETER :: corrosion_delay_time = 0.0_dp    ! Delay time before corrosion starts [days]

      ! --- Wear Model Parameters (Tribology) ---
      REAL(KIND=dp), PARAMETER :: slips_per_day      = 1393548_dp    ! [cycles/day] Sliding frequency
      REAL(KIND=dp), PARAMETER :: slip_amplitude     = 1.0_dp     ! [mm] Sliding amplitude
      REAL(KIND=dp), PARAMETER :: test_pressure      = 1000.0_dp    ! [MPa] Reference contact pressure
      REAL(KIND=dp), PARAMETER :: friction_coefficient = 0.496_dp   ! Friction coefficient [dimensionless]
      REAL(KIND=dp), PARAMETER :: wear_coefficient   = 1.72e-13_dp  ! Wear coefficient [mm^3/(N*mm)]
      REAL(KIND=dp), PARAMETER :: effective_pressure_exponent = 50.0_dp     ! Exponent for contact pressure in the wear model [dimensionless]

      ! --- Interaction parameters (tribocorrosion) ---
      REAL(KIND=dp), PARAMETER :: coupling_factor_alpha = 0.5_dp ! Synergy factor [dimensionless]

      ! --- Physical properties and control parameters ---
      REAL(KIND=dp), PARAMETER :: mass        = 7.86e-6_dp ! [ton/mm^3] Material density
      REAL(KIND=dp), PARAMETER :: corr        = 1.e20_dp    ! Factor for corroded-mass criterion 
      

!-------------------------------------------------------------------------------------------------------------------------------------
  ! --- 3. DECLARATION OF INDICES FOR STATE VARIABLES (SDVs) ---
  ! SDVs (Solution-Dependent State Variables) store material state at each integration point.
  !-------------------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: SDV_PEEQ                 = 1  ! Equivalent plastic strain
      INTEGER, PARAMETER :: SDV_TOTAL_DAMAGE         = 2  ! Total accumulated damage
      INTEGER, PARAMETER :: SDV_DELETE_FLAG          = 3  ! Element deletion flag
      INTEGER, PARAMETER :: SDV_PITTING_FACTOR       = 4  ! Pitting factor
      INTEGER, PARAMETER :: SDV_CORROSION_STATUS     = 5  ! Corrosion status (active/inactive)
      INTEGER, PARAMETER :: SDV_NONLOCAL_PROPERTY    = 6  ! Non-local property
      INTEGER, PARAMETER :: SDV_STRESS_S11           = 7  ! Stress component S11
      INTEGER, PARAMETER :: SDV_STRESS_S22           = 8  ! Stress component S22
      INTEGER, PARAMETER :: SDV_STRESS_S33           = 9  ! Stress component S33
      INTEGER, PARAMETER :: SDV_STRESS_S12           = 10 ! Stress component S12
      INTEGER, PARAMETER :: SDV_STRESS_S13           = 11 ! Stress component S13
      INTEGER, PARAMETER :: SDV_STRESS_S23           = 12 ! Stress component S23
      INTEGER, PARAMETER :: SDV_VON_MISES            = 13 ! Von Mises stress
      INTEGER, PARAMETER :: SDV_DCOld                = 14 ! (Not used, kept for compatibility)
      INTEGER, PARAMETER :: SDV_DSCLd                = 15 ! (Not used, kept for compatibility)
      INTEGER, PARAMETER :: SDV_NONLOCAL_STRESS      = 16 ! Non-local stress
      INTEGER, PARAMETER :: SDV_SURFACE_FLAG         = 17 ! Surface flag
      INTEGER, PARAMETER :: SDV_TOTAL_CORR_DAMAGE    = 18 ! Total corrosion damage
      INTEGER, PARAMETER :: SDV_CURRENT_TOTAL_DAMAGE = 19 ! Current total damage
      INTEGER, PARAMETER :: SDV_JC_DAMAGE_INC        = 20 ! Johnson-Cook damage increment
      INTEGER, PARAMETER :: SDV_CORROSION_DEPTH      = 21 ! Corrosion depth
      INTEGER, PARAMETER :: SDV_CORROSION_START_TIME = 22 ! Corrosion start time
      INTEGER, PARAMETER :: SDV_CORROSION_INIT_FLAG  = 24 ! Corrosion initialization flag
      INTEGER, PARAMETER :: SDV_TRIBO_DAMAGE_ACCUM   = 26 ! Accumulated tribological damage
      INTEGER, PARAMETER :: SDV_CORRODED_VOLUME      = 30 ! Corroded volume
      
  !-----------------------------------------------------------------------
  ! --- 4. DYNAMIC ARRAY ALLOCATION ---
  !-----------------------------------------------------------------------
      ALLOCATE (damage_new(n_elems))
      ALLOCATE (relative_slip_increment(n_elems))
      ALLOCATE (contact_pressure(n_elems))

  !-----------------------------------------------------------------------
  ! --- 5. OBTAINING ELASTIC PROPERTIES ---
  ! props(1): Young's modulus, props(2): Poisson's ratio
  ! props(3..) - Yield and hardening data
  !-----------------------------------------------------------------------
      e         = props(1)
      xnu       = props(2)
      twomu     = e / (one + xnu)
      alamda    = xnu * twomu / (one - (two * xnu))
      thremu    = one_point_five * twomu
      nvalue    = (nprops / 2) - 1
      masslimit = mass * corr

  !-----------------------------------------------------------------------
  ! --- 6. MAIN MATERIAL LOGIC ---
  !-----------------------------------------------------------------------
      IF (totalTime .eq. zero) THEN
        ! --- Initial Elastic Guess (First Increment) ---
        ! Essential for initial stability in explicit simulations.
        DO k = 1, nblock
            trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
            stressNew(k,1) = stressOld(k,1) + twomu * strainInc(k,1) + alamda * trace
            stressNew(k,2) = stressOld(k,2) + twomu * strainInc(k,2) + alamda * trace
            stressNew(k,3) = stressOld(k,3) + twomu * strainInc(k,3) + alamda * trace
            stressNew(k,4) = stressOld(k,4) + twomu * strainInc(k,4)
            IF (nshr > 1) THEN
                stressNew(k,5) = stressOld(k,5) + twomu * strainInc(k,5)
                stressNew(k,6) = stressOld(k,6) + twomu * strainInc(k,6)
            END IF
        END DO
        lmCO = 0.0

      ELSE
        ! --- Start of calculations for subsequent increments ---
        deltat = nt * dt

        ! --- Main loop over integration points ---
        DO k = 1, nblock

              ! --- Loop variable initialization ---
            corrosion_damage      = 0.0_dp
            tribological_damage   = 0.0_dp
            tribological_damage_inc = 0.0_dp
            deqps                 = 0.0_dp
            s13                   = 0.0_dp
            s23                   = 0.0_dp

            ! Retrieve previous increment values of state variables
            peeqOld_k = stateOld(k, SDV_PEEQ)
            dCOld_k   = stateOld(k, SDV_DCOld)
            dSCld_k   = stateOld(k, SDV_DSCLd)
            damage_new(nElement(k)) = stateOld(k, SDV_TOTAL_DAMAGE)
      

  !===============================================================
  ! --- 7. DAMAGE CALCULATION (CORROSION, WEAR, TRIBOCORROSION) ---
  !===============================================================
      
        IF ((totalTime < one) .AND. (element_corrosion_status(nElement(k)) /= flag_corrosion_inactive)) THEN
        
!  Modification for galvanic corrosion / tribology
        IF (totalTime > corrosion_delay_time) then 
! Mises stress for contact simulation (test)
      IF (nshr == 1) THEN
      vmises = SQRT(one_point_five * (stressOld(k,1)*stressOld(k,1) + 
     1 stressOld(k,2)*stressOld(k,2) + stressOld(k,3)*stressOld(k,3) + two*stressOld(k,4)*stressOld(k,4)))
      ELSE
      vmises = SQRT(one_point_five * (stressOld(k,1)*stressOld(k,1) + 
     1 stressOld(k,2)*stressOld(k,2) + stressOld(k,3)*stressOld(k,3) + two*stressOld(k,4)*stressOld(k,4)
     2   + two*stressOld(k,6)*stressOld(k,6) + two*stressOld(k,5)*stressOld(k,5)))
      END IF        
              
        !-----------------------------------------------------------
        ! --- 7.1. PURE CORROSION DAMAGE ---
        !-----------------------------------------------------------
          ! tnt: adjusted time considering delay (corrosion_delay_time)
          tnt = (totalTime-corrosion_delay_time+dt) * nt                 
          ku = 1.011*((tnt-stateNew(k,SDV_CORROSION_START_TIME))/365.)**-0.05   
          ! ku = 0.0_dp !Disable corrosion        
         
        ! If this is the first corrosion increment, store the start time
        IF (stateOld(k, SDV_CORROSION_INIT_FLAG) /= flag_corrosion_started) THEN
          stateNew(k, SDV_CORROSION_INIT_FLAG) = flag_corrosion_started
          stateNew(k, SDV_CORROSION_START_TIME) = tnt
        END IF
                   
         
                
        ! Compute corrosion depth
        stateNew(k, SDV_CORROSION_DEPTH) = ABS(corrosion_rate *ku* (tnt - stateNew(k, SDV_CORROSION_START_TIME)) / 365.0)
     1  * elem_properties(nElement(k),PROP_IDX_SURFACE_FLAG)

        IF (stateNew(k, SDV_CORROSION_DEPTH) < stateOld(k, SDV_CORROSION_DEPTH)) THEN
            stateNew(k, SDV_CORROSION_DEPTH) = stateOld(k, SDV_CORROSION_DEPTH)
        END IF
        
        

        ! Compute corroded volume and corresponding damage
        elem_length(nElement(k)) = elem_properties(nElement(k), PROP_IDX_VOLUME) ** (1.0/3.0)             
        corroded_volume(nElement(k)) = ABS(elem_length(nElement(k)) * elem_length(nElement(k))
     1  * stateNew(k, SDV_CORROSION_DEPTH)) * nonlocal_property(nElement(k))
         
         
        initial_volume = elem_properties(nElement(k), PROP_IDX_VOLUME)

        IF (initial_volume > 0.0_dp) THEN
          IF (corroded_volume(nElement(k)) >= initial_volume) THEN
            corrosion_damage = dmax
            stateNew(k, SDV_CORROSION_DEPTH) = 0.0_dp
          ELSE
            corrosion_damage = (ABS(corroded_volume(nElement(k))) ) / initial_volume
          END IF
        ELSE
        corrosion_damage = 0.0_dp
        END IF

        IF (vmises>10.0) THEN
    !-----------------------------------------------------------
    ! --- 7.2. PURE WEAR DAMAGE ---
    !-----------------------------------------------------------
    ! Conceptual energy-based wear model
        calculated_slip_inc = slips_per_day * slip_amplitude   * deltat 
        relative_slip_increment(nElement(k)) = calculated_slip_inc
        contact_pressure(nElement(k)) = sqrt(stressOld(k,1)**2.0 + stressOld(k,3)**2.0) 

        IF (relative_slip_increment(nElement(k)) > zero) THEN
        contact_area = elem_properties(nElement(k), PROP_IDX_VOLUME)**(2.0/3.0)
        force_of_friction = effective_pressure_exponent * contact_pressure(nElement(k)) * friction_coefficient * contact_area
        energy_dissipated_increment = force_of_friction * relative_slip_increment(nElement(k))
        wear_volume_increment = wear_coefficient * energy_dissipated_increment

              ! Convert wear volume into incremental damage
              IF (elem_properties(nElement(k), PROP_IDX_VOLUME) > zero) THEN
                tribological_damage_inc = wear_volume_increment / elem_properties(nElement(k), PROP_IDX_VOLUME)            
              ELSE
                tribological_damage = zero
              END IF
        END IF
! End of tribological/galvanic condition
      END IF
    !-----------------------------------------------------------
    ! --- 7.3. COMBINED DAMAGE (TRIBOCORROSION) ---
    !-----------------------------------------------------------
      ! Retrieve previously accumulated tribological damage
      accumulated_pure_tribo = stateOld(k, SDV_TRIBO_DAMAGE_ACCUM)

      ! Synergy: corrosion accelerates the new wear increment
      enhanced_tribo_damage_inc = tribological_damage_inc * (one + coupling_factor_alpha * corrosion_damage)
                
      ! Accumulate wear damage (including synergy)
      accumulated_pure_tribo = accumulated_pure_tribo + enhanced_tribo_damage_inc

      ! Total damage is the sum of pure corrosion + accumulated wear (with synergy)
      combined_damage = corrosion_damage + accumulated_pure_tribo
      total_damage = MIN(combined_damage, dmax)
      damage_new(nElement(k)) = total_damage
      timetot = totalTime - corrosion_delay_time
        
        
  
! Time for corrosion - contact occurrence
      END IF ! End of corrosion/tribology condition
  !-----------------------------------------------------------
  ! --- 8. CHECK MASS LOSS AND ELEMENT FAILURE DUE TO DAMAGE ---
  !-----------------------------------------------------------
    ! Stop condition by corroded mass
      IF (mCO >= masslimit) THEN
        !write(*,*) mCO, timetot
        !write(*,*) 'End of corrosion: 100% of element volume affected'
        CALL block_all(ok)
      END IF

    ! Element failure condition due to maximum damage
        IF (damage_new(nElement(k)) >= dmax) THEN
          stateNew(k, SDV_DELETE_FLAG) = delete_element_flag
          damage_new(nElement(k))      = dmax
          element_corrosion_status(nElement(k)) = flag_corrosion_inactive
          prop_update_trigger = flag_general_property_transfer
          CALL process_failed_element_neighbors(nElement(k), beta, OK)
          failure_occurred_this_increment = .TRUE.
          stateNew(k, SDV_CORROSION_DEPTH) = 0.0
        END IF
      END IF ! End of damage calculation block

    ! Ensure damage never decreases
      IF (damage_new(nElement(k)) < stateOld(k, SDV_TOTAL_DAMAGE)) THEN
          damage_new(nElement(k)) = stateOld(k, SDV_TOTAL_DAMAGE)
      END IF

  !===============================================================
  ! --- 9. PLASTICITY CALCULATION (VON MISES CRITERION) ---
  !===============================================================
      CALL vuhard(yieldOld, hard, peeqOld_k, props(3), nvalue)

    ! Elastic trial stress
      trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
      s11 = stateOld(k, SDV_STRESS_S11) + twomu * strainInc(k,1) + alamda * trace
      s22 = stateOld(k, SDV_STRESS_S22) + twomu * strainInc(k,2) + alamda * trace
      s33 = stateOld(k, SDV_STRESS_S33) + twomu * strainInc(k,3) + alamda * trace
      s12 = stateOld(k, SDV_STRESS_S12) + twomu * strainInc(k,4)
      IF (nshr > 1) THEN
        s13 = stateOld(k, SDV_STRESS_S13) + twomu * strainInc(k,5)
        s23 = stateOld(k, SDV_STRESS_S23) + twomu * strainInc(k,6)
      END IF

    ! Deviatoric stress and Von Mises stress
      smean = third * (s11 + s22 + s33)
      s11 = s11 - smean
      s22 = s22 - smean
      s33 = s33 - smean

      IF (nshr == 1) THEN
      vmises = SQRT(one_point_five * (s11*s11 + s22*s22 + s33*s33 + two*s12*s12))
      ELSE
      vmises = SQRT(one_point_five * (s11*s11 + s22*s22 + s33*s33 + two*s12*s12
     1   + two*s13*s13 + two*s23*s23))
      END IF

    ! Equivalent plastic strain increment
      sigdif = vmises - yieldOld
      IF (sigdif > zero) THEN
      deqps = sigdif / (thremu + hard)
      END IF 
  !===============================================================
  ! --- 10. UPDATE OF MECHANICAL DAMAGE AND STRESS ---
  !===============================================================
  ! --- Stress Update (Return Mapping with Damage) ---
      yieldNew = yieldOld + hard * deqps
      IF ((yieldNew + thremu * deqps) > 0.0_dp) THEN
          factor = yieldNew / (yieldNew + thremu * deqps)
      ELSE
          factor = 1.0_dp
      END IF

      stressNew(k,1) = (s11 * factor + smean) * (one - damage_new(nElement(k)))
      stressNew(k,2) = (s22 * factor + smean) * (one - damage_new(nElement(k)))
      stressNew(k,3) = (s33 * factor + smean) * (one - damage_new(nElement(k)))
      stressNew(k,4) = (s12 * factor) * (one - damage_new(nElement(k)))
      IF (nshr > 1) THEN
          stressNew(k,5) = (s13 * factor) * (one - damage_new(nElement(k)))
          stressNew(k,6) = (s23 * factor) * (one - damage_new(nElement(k)))
      END IF

  ! --- Removal due to excessive plastic deformation ---
      eqps = stateOld(k, SDV_PEEQ) + deqps
      IF ((eqps > 0.137_dp) .AND. (strain_deletion_status(nElement(k)) /= flag_deletion_by_strain_disabled)) THEN
        stateNew(k, SDV_DELETE_FLAG) = delete_element_flag
        CALL process_failed_element_neighbors(nElement(k), beta, OK)
        CALL compute_nonlocal_property(ok)
        prop_update_trigger = flag_general_property_transfer
      END IF

    !===============================================================
    ! --- 11. UPDATE OF STATE VARIABLES (stateNew) ---
    !===============================================================
      von_mises_stress(nElement(k)) = vmises
      damage_current(nElement(k))   = damage_new(nElement(k))
      
      

      stateNew(k, SDV_PEEQ)                 = eqps
      stateNew(k, SDV_TOTAL_DAMAGE)         = damage_new(nElement(k))
      stateNew(k, SDV_PITTING_FACTOR)       = elem_properties(nElement(k), PROP_IDX_PITTING)
      stateNew(k, SDV_CORROSION_STATUS)     = element_corrosion_status(nElement(k))
      stateNew(k, SDV_NONLOCAL_PROPERTY)    = nonlocal_property(nElement(k))
      stateNew(k, SDV_STRESS_S11)           = s11 * factor + smean
      stateNew(k, SDV_STRESS_S22)           = s22 * factor + smean
      stateNew(k, SDV_STRESS_S33)           = s33 * factor + smean
      stateNew(k, SDV_STRESS_S12)           = s12 * factor
      stateNew(k, SDV_STRESS_S13)           = s13 * factor
      stateNew(k, SDV_STRESS_S23)           = s23 * factor
      stateNew(k, SDV_VON_MISES)            = von_mises_stress(nElement(k))
      stateNew(k, SDV_DCOld)                = dCOld_k
      stateNew(k, SDV_DSCLd)                = dSCld_k
      stateNew(k, SDV_NONLOCAL_STRESS)      = non_local_stress(nElement(k))
      stateNew(k, SDV_SURFACE_FLAG)         = elem_properties(nElement(k), PROP_IDX_SURFACE_FLAG)
      stateNew(k, SDV_TOTAL_CORR_DAMAGE)    = corrosion_damage
      stateNew(k, SDV_CURRENT_TOTAL_DAMAGE) = damage_current(nElement(k))
      stateNew(k, SDV_JC_DAMAGE_INC)        = tribological_damage
      stateNew(k, SDV_TRIBO_DAMAGE_ACCUM)   = accumulated_pure_tribo
      stateNew(k, SDV_CORRODED_VOLUME)      = corroded_volume(nElement(k))      

  !===============================================================
  ! --- 12. UPDATE OF ENERGIES ---
  !===============================================================
      IF (nshr == 1) THEN
        stressPower = half * (
     1 (stateOld(k, SDV_STRESS_S11) + stressNew(k,1)) * strainInc(k,1) + +
     2 (stateOld(k, SDV_STRESS_S22) + stressNew(k,2)) * strainInc(k,2) +
     3 (stateOld(k, SDV_STRESS_S33) + stressNew(k,3)) * strainInc(k,3) +
     4 (stateOld(k, SDV_STRESS_S12) + stressNew(k,4)) * strainInc(k,4) )
      ELSE
        stressPower = half * (
     1  (stateOld(k, SDV_STRESS_S11) + stressNew(k,1)) * strainInc(k,1) +
     2  (stateOld(k, SDV_STRESS_S22) + stressNew(k,2)) * strainInc(k,2) +
     3  (stateOld(k, SDV_STRESS_S33) + stressNew(k,3)) * strainInc(k,3) +
     4  (stateOld(k, SDV_STRESS_S12) + stressNew(k,4)) * strainInc(k,4) +
     5  (stateOld(k, SDV_STRESS_S13) + stressNew(k,5)) * strainInc(k,5) +
     6  (stateOld(k, SDV_STRESS_S23) + stressNew(k,6)) * strainInc(k,6) )
      END IF
      enerInternNew(k) = enerInternOld(k) + stressPower / density(k)

      plasticWorkInc = half * (yieldOld + yieldNew) * deqps
      enerInelasNew(k) = enerInelasOld(k) + plasticWorkInc / density(k)

        END DO ! End of main loop over integration points
      END IF ! End of main material logic
                 
      ! --- Deallocate local arrays ---
      DEALLOCATE (damage_new)
      DEALLOCATE (relative_slip_increment)
      DEALLOCATE (contact_pressure)
      RETURN           
            
      END SUBROUTINE vumatXtrArg
!-----------------------------------------------------------------------      
!> @brief Computes the yield stress and isotropic hardening.
!!
!! Performs a linear interpolation on the stress vs. equivalent plastic strain curve
!! provided in the material property table.
!-----------------------------------------------------------------------
      SUBROUTINE vuhard(syield, hard, eqplas, table, nvalue)
      USE SharedVariables
      USE AbaqusParameters

      ! --- Argumentos ---
      INTEGER,       INTENT(IN)  :: nvalue
      REAL(KIND=dp), INTENT(IN)  :: eqplas, table(2,*)
      REAL(KIND=dp), INTENT(OUT) :: syield, hard

      ! --- Local variables ---
      INTEGER :: k1
      REAL(kind=dp) :: eqpl1, eqpl0, deqpl, syiel0, syiel1, dsyiel
      REAL(kind=dp), PARAMETER :: zero = 0.0_dp

    ! Initialize output values assuming right extrapolation
    ! (uses the last value in the table if strain exceeds all).
      syield = table(1, nvalue)
      hard   = zero

    ! If table has more than one point, search for correct interval.
      IF (nvalue > 1) THEN
        DO k1 = 1, nvalue - 1
          eqpl1 = table(2, k1 + 1)

          ! If plastic strain is within interval [eqpl0, eqpl1]
          IF (eqplas < eqpl1) THEN
              eqpl0  = table(2, k1)
              syiel0 = table(1, k1)
              syiel1 = table(1, k1 + 1)

              deqpl  = eqpl1 - eqpl0
              dsyiel = syiel1 - syiel0

            ! Computes hardening modulus (derivative)
            IF (deqpl > 0.0_dp) THEN
                hard = dsyiel / deqpl
            ELSE
                hard = 0.0_dp ! Avoid division by zero
            END IF

            ! Linear interpolation for yield stress
            syield = syiel0 + (eqplas - eqpl0) * hard
            EXIT ! Exit loop since correct interval was found
          END IF
        END DO
      END IF
      RETURN
      END SUBROUTINE vuhard

!----------------------------------------------------------------------------------------------------------
!> @brief Entry point for data exchange and control with Abaqus.
!!
!! This subroutine is called by Abaqus at different analysis points (start, end increment, etc.)
!! to enable variable initialization, data reading, and post-processing routine execution.
!----------------------------------------------------------------------------------------------------------   
      SUBROUTINE vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

      USE AbaqusParameters
      USE SharedVariables      

      ! --- Arguments ---
      INTEGER,       INTENT(IN) :: lOp, niArray, nrArray, i_Array(niArray)
      REAL(KIND=dp), INTENT(IN) :: r_Array(nrArray)

      ! --- Local variables ---
      INTEGER :: k, p, l, ok
      INTEGER :: kStep, kInc, iStatus, nElems
      INTEGER :: line_count, read_status
      INTEGER, PARAMETER :: unit_num = 100
      CHARACTER(LEN=512) :: line_buffer  

      ! --- Initialization of control variables ---
      kStep   = i_Array(I_INT_KSTEP)
      kInc    = i_Array(I_INT_KINC)
      iStatus = i_Array(I_INT_ISTATUS)
      nElems  = i_Array(I_INT_NTOTALELEMENTS)     

      !=======================================================================
      ! --- ANALYSIS EVENT SELECTION ---
      !=======================================================================
      SELECT CASE (lOp)     

      !-----------------------------------------------------------------------
      ! --- Event: Start of analysis ---
      !-----------------------------------------------------------------------
      CASE (J_INT_START_ANALYSIS)
        ! --- Read properties file in two passes ---
        ! 1. Open file and count lines to determine n_elems
      OPEN(unit_num, file='element_properties.txt', ! Best practice: use absolute path or ensure file is in working directory
     1           status='old', action='read', iostat=read_status)           

      ! Error handling on file open
      IF (read_status /= 0) THEN
          PRINT *, '============================================================'
          PRINT *, '>>> FATAL ERROR IN VEXTERNALDB SUBROUTINE <<<'
          PRINT *, '>>> Could not open properties file.'
          PRINT *, '>>> Verify file path and existence.'
          PRINT *, '============================================================'
          STOP 'Execution stopped due to file open error.'
      END IF

        line_count = 0

      DO
          READ(unit_num, '(A)', iostat=read_status) line_buffer
          IF (read_status /= 0) EXIT ! Exit on error or end of file
          IF (TRIM(line_buffer) /= '') THEN
              line_count = line_count + 1
          END IF
        END DO
        REWIND(unit_num)
        
      ! 2. Set n_elems, allocate memory, and read data
        n_elems = line_count
        PRINT *, 'Diagnostics: Properties file read. Number of elements: ', n_elems
        CALL allocate_arrays(ok)

        DO k = 1, n_elems
            READ(unit_num, *) (elem_properties(k,l), l=1,14)
        END DO
        CLOSE(unit_num)

        ! --- Initialization of properties and data structures ---
        CALL block_surface(ok)
        DO p = 1, n_elems
            IF (elem_properties(p,PROP_IDX_SURFACE_FLAG) == 1.0_dp) THEN
                temp_pitting_write(p) = elem_properties(p,PROP_IDX_PITTING)
                temp_surface_flag_write(p) = elem_properties(p,PROP_IDX_SURFACE_FLAG)
            ELSE
                elem_properties(p,PROP_IDX_PITTING) = 0.0_dp
                temp_pitting_write(p) = elem_properties(p,PROP_IDX_PITTING)
            END IF
        END DO
        
        prop_update_trigger = flag_general_property_transfer
        CALL build_influence_map(ok)
        CALL compute_nonlocal_property(ok)
        CALL transfer_properties(ok)

      !-----------------------------------------------------------------------
      ! --- Event: End of each increment ---
      !-----------------------------------------------------------------------
      CASE (J_INT_END_INCREMENT)
        ! If element failure occurred, update neighbor properties
        IF (failure_occurred_this_increment) THEN
            CALL transfer_properties(ok)
            CALL compute_nonlocal_property(ok)
            failure_occurred_this_increment = .FALSE. ! Reset flag
        END IF

        ! Routines executed each increment
        CALL compute_mass_loss(ok)
        CALL compute_nonlocal_stress(ok)

      !-----------------------------------------------------------------------
      ! --- Event: End of analysis ---
      !-----------------------------------------------------------------------
      CASE (J_INT_END_ANALYSIS)
        CALL deallocate_arrays(ok)
        CLOSE(unit_num)

      END SELECT
      RETURN
      END SUBROUTINE vexternaldb
      
!----------------------------------------------------------------------------------------------------------------------------
!> @brief Updates a property for neighbors of a failed element.
!!
!! Propagates the effect of failure of element ('nElement') to its direct neighbors. Operates on a specific column of matrix
!! 'elem_properties' and stores result in a temporary array, avoiding code duplication.
!----------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE update_neighbor_property(nElement, beta, prop_column_index, target_write_array, ok)
        USE SharedVariables
        IMPLICIT NONE

      ! --- Argumentos ---
        INTEGER,       INTENT(IN)    :: nElement
        REAL(KIND=dp), INTENT(IN)    :: beta
        INTEGER,       INTENT(IN)    :: prop_column_index
        REAL(KIND=dp), INTENT(INOUT) :: target_write_array(n_elems)
        INTEGER,       INTENT(OUT)   :: ok

      ! --- Local variables ---
        INTEGER       :: j, neighbor_element
        REAL(KIND=dp) :: neighbor_prop_old, prop_from_failed

      ! --- Main logic ---
      ! Retrieve the property from the failed element for propagation to neighbors.
        prop_from_failed = elem_properties(nElement, prop_column_index)

      ! Iterate over the 6 neighbors of the failed element (assumes structured mesh).
        DO j = 1, max_face_neighbors
          ! Columns 2 to 7 of elem_properties contain neighbor IDs.
          neighbor_element = elem_properties(nElement, j + 1)

          ! Proceed only if neighbor exists (ID not zero).
          IF (neighbor_element /= 0) THEN
              ! If neighbor is not "locked", calculate and update property.
            IF (property_update_lock(neighbor_element) /= key_prop_update_locked) THEN
                neighbor_prop_old = elem_properties(neighbor_element, prop_column_index)

                ! MAX logic ensures neighbor property only increases or stays same,
                ! never decreases due to another neighbor's failure.
                target_write_array(neighbor_element) = MAX(neighbor_prop_old, prop_from_failed * beta)
            END IF
          END IF
        END DO

      ! Zero the property in the failed element and lock it for
      ! future neighbor updates.
        target_write_array(nElement) = 0.0_dp
        property_update_lock(nElement) = key_prop_update_locked

        ok = 1
        RETURN
      END SUBROUTINE update_neighbor_property
!-----------------------------------------------------------------------------------------------------
!> @brief Orchestrates all updates for neighbors of a failed element.
!!
!! When an element fails in VUMAT, this routine ensures failure consequences
!! (pitting update and surface flag) are propagated correctly to neighbors.
!-----------------------------------------------------------------------------------------------------
      SUBROUTINE process_failed_element_neighbors(nElement, beta, ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Arguments ---
        INTEGER,       INTENT(IN)  :: nElement
        REAL(KIND=dp), INTENT(IN)  :: beta
        INTEGER,       INTENT(OUT) :: ok
        
        ! This routine encapsulates the two actions always occurring together at failure.

        ! 1. Update surface flag property (column 8).
        CALL update_neighbor_property(nElement, beta, PROP_IDX_SURFACE_FLAG, temp_surface_flag_write, ok)

        ! 2. Update pitting factor property (column 9).
        CALL update_neighbor_property(nElement, beta, PROP_IDX_PITTING, temp_pitting_write, ok)
        
        ! Return OK. Error handling logic can be improved.
        ok = 1  
        RETURN
      END SUBROUTINE process_failed_element_neighbors
! -------------------------------------------------------------------------
!> @brief Computes the non-local (homogenized) property for surface elements.
!!
!! Implements a weighted-average model to compute a non-local property,
!! considering the influence of neighbor elements within radius 'lr'.
! -------------------------------------------------------------------------
      SUBROUTINE compute_nonlocal_property(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Arguments ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Local variables ---
        INTEGER       :: i, p, neighbor_element
        REAL(KIND=dp) :: denominator_sum, numerator_sum, weight_p
        REAL(KIND=dp) :: dist_sq, lr_sq

        ! --- Main logic ---
        ! Precompute the square of the influence radius for performance.
        lr_sq = lr**2.0_dp

        DO i = 1, n_elems
            ! Proceed only if element 'i' is a surface element (pitting property not zero).
            IF (elem_properties(i, PROP_IDX_PITTING) /= 0.0_dp) THEN

            ! Initialize sums with the contribution of the element itself.
              denominator_sum = elem_properties(i, PROP_IDX_VOLUME) ! Denominator: weighted sum of volumes
              numerator_sum   = elem_properties(i, PROP_IDX_PITTING) * elem_properties(i, PROP_IDX_VOLUME) ! Numerator: weighted sum of (prop * vol)

            ! Iterate over neighbors to accumulate their contributions.
              DO p = 1, counter(i)
              neighbor_element = influ(i, p)

                ! Consider only neighbors that are also surface elements.
                IF (elem_properties(neighbor_element, PROP_IDX_PITTING) /= 0.0_dp) THEN
                    dist_sq = dist(i, p)**2.0_dp
                        
                    ! Quadratic weight function
                    weight_p = (1.0_dp - (dist_sq / lr_sq))**2.0_dp

                    denominator_sum = denominator_sum + (weight_p * elem_properties(neighbor_element, PROP_IDX_VOLUME))
                    numerator_sum = numerator_sum + (weight_p * elem_properties(neighbor_element, PROP_IDX_PITTING)
     1                                                        * elem_properties(neighbor_element, PROP_IDX_VOLUME))
                END IF
              END DO

                ! Computes final non-local property value (weighted average).
                IF (denominator_sum > 1.0e-20_dp) THEN
                    temp_nonlocal_write(i) = numerator_sum / denominator_sum
                ELSE
                    temp_nonlocal_write(i) = 0.0_dp
                END IF
            ELSE
            ! For internal elements, non-local property is zero.
                temp_nonlocal_write(i) = 0.0_dp
            END IF
        END DO

        ok = 1
        RETURN
      END SUBROUTINE compute_nonlocal_property
!-----------------------------------------------------------------------------------------------------------------
!> @brief Transfers computed properties from temporary arrays to globals.
!!
!! This routine acts as a gate, copying data from temporary vectors (e.g. temp_nonlocal_write)
!! into persistent state variables (e.g. nonlocal_property). Execution is controlled by a global trigger.
!-----------------------------------------------------------------------------------------------------------------
      SUBROUTINE transfer_properties(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Arguments ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Local variables ---
          INTEGER :: i

            ! --- Main logic ---
            ! Transfer only occurs if global trigger is active.
            IF (prop_update_trigger == flag_general_property_transfer) THEN

              ! Copy data from temporary arrays into global arrays.
              DO i = 1, n_elems
                nonlocal_property(i) = temp_nonlocal_write(i)
                elem_properties(i,PROP_IDX_PITTING) = temp_pitting_write(i)
                elem_properties(i,PROP_IDX_SURFACE_FLAG) = temp_surface_flag_write(i)
              END DO

              ! Disable trigger to avoid unnecessary re-execution.
              ! Will be reactivated when a new element failure occurs.
              prop_update_trigger = flag_general_property_transfer_locked
          END IF

          ok = 1
          RETURN
      END SUBROUTINE transfer_properties   
!------------------------------------------------------------------------------------------------------------------
!> @brief Builds the influence map for each element in the model.
!!
!! For every element 'i', this routine finds all other elements 'j' whose Euclidean distance
!! is less than the intrinsic length 'lr'. It stores neighbor IDs, distances and neighbor counts.
!------------------------------------------------------------------------------------------------------------------    
      SUBROUTINE build_influence_map(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Arguments ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Local variables ---
        INTEGER       :: t, u, i, j, m
        REAL(KIND=dp) :: distance

        ! --- 1. Initialize influence vectors and matrices ---
        DO t = 1, n_elems
            counter(t) = 0
            DO u = 1, max_influence_elems
                dist(t,u)  = 0.0_dp
                influ(t,u) = 0
            END DO
        END DO

        ! --- 2. Compute distance between each pair of elements (i, j) ---
        ! This is an O(N^2) loop and may be slow for very large meshes.
        DO i = 1, n_elems
            m = 1 ! Restart neighbor counter for element 'i'
            DO j = 1, n_elems
                IF (j /= i) THEN

                ! Euclidean distance calculation
                distance = SQRT((elem_properties(j,PROP_IDX_COORD_X) - elem_properties(i,PROP_IDX_COORD_X))**2.0 +
     1                          (elem_properties(j,PROP_IDX_COORD_Y) - elem_properties(i, PROP_IDX_COORD_Y))**2.0 +
     2                          (elem_properties(j,PROP_IDX_COORD_Z) - elem_properties(i,PROP_IDX_COORD_Z))**2.0 )

                ! If 'j' is within influence radius, store it as a neighbor.
                  IF (distance < lr) THEN
                          ! Check if number of neighbors does not exceed allocated limit.
                          IF (m <= max_influence_elems) THEN
                              influ(i,m)   = j
                              dist(i,m)    = distance
                              counter(i)   = m
                              m = m + 1
                          ELSE                        
                          END IF
                  END IF
                END IF
            END DO
        END DO

        ok = 1
        RETURN
      END SUBROUTINE build_influence_map

!--------------------------------------------------------------------------------------------------------------
!> @brief Safely allocate memory for shared variables.
!!
!! Ensures memory is only allocated if not already allocated, avoiding
!! "array is already allocated" errors.
!--------------------------------------------------------------------------------------------------------------
      SUBROUTINE allocate_arrays(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Arguments ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Local variables ---
          INTEGER :: stat ! Variable to check allocation status

          ! --- Safe allocation of all global arrays ---
          IF (.NOT. ALLOCATED(elem_properties))          ALLOCATE(elem_properties(n_elems, 14), stat=stat)
          IF (.NOT. ALLOCATED(nonlocal_property))        ALLOCATE(nonlocal_property(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(property_update_lock))     ALLOCATE(property_update_lock(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(influ))                    ALLOCATE(influ(n_elems, max_influence_elems), stat=stat)
          IF (.NOT. ALLOCATED(dist))                     ALLOCATE(dist(n_elems, max_influence_elems), stat=stat)
          IF (.NOT. ALLOCATED(temp_pitting_write))       ALLOCATE(temp_pitting_write(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(temp_surface_flag_write))  ALLOCATE(temp_surface_flag_write(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(counter))                  ALLOCATE(counter(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(temp_nonlocal_write))      ALLOCATE(temp_nonlocal_write(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(damage_current))           ALLOCATE(damage_current(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(element_corrosion_status)) ALLOCATE(element_corrosion_status(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(strain_deletion_status))   ALLOCATE(strain_deletion_status(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(von_mises_stress))         ALLOCATE(von_mises_stress(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(non_local_stress))         ALLOCATE(non_local_stress(n_elems), stat=stat)

            ! Check allocation status
             IF (stat /= 0) THEN
               PRINT *, 'Error allocating memory!'
               STOP
             END IF

          ok = 1
          RETURN
      END SUBROUTINE allocate_arrays

!--------------------------------------------------------------------------------------------------------------
!> @brief Deallocate shared variables memory at the end of the analysis.
!!
!! Ensures all dynamically allocated memory is released to avoid memory leaks.
!--------------------------------------------------------------------------------------------------------------           
      SUBROUTINE deallocate_arrays(ok)
        USE SharedVariables
        IMPLICIT NONE
        INTEGER, intent(out) :: ok

        ! --- Desalocação Segura de todos os Arrays Globais ---
        IF (allocated(elem_properties))          DEALLOCATE(elem_properties)
        IF (allocated(nonlocal_property))        DEALLOCATE(nonlocal_property)
        IF (allocated(property_update_lock))     DEALLOCATE(property_update_lock)
        IF (allocated(influ))                    DEALLOCATE(influ) 
        IF (allocated(dist))                     DEALLOCATE(dist)
        IF (allocated(temp_pitting_write))       DEALLOCATE(temp_pitting_write)
        IF (allocated(temp_surface_flag_write))  DEALLOCATE(temp_surface_flag_write)
        IF (allocated(counter))                  DEALLOCATE(counter)
        IF (allocated(temp_nonlocal_write))      DEALLOCATE(temp_nonlocal_write)
        IF (allocated(damage_current))           DEALLOCATE(damage_current)
        IF (allocated(element_corrosion_status)) DEALLOCATE(element_corrosion_status)
        IF (allocated(strain_deletion_status))   DEALLOCATE(strain_deletion_status)
        IF (allocated(von_mises_stress))         DEALLOCATE(von_mises_stress)
        IF (allocated(non_local_stress))         DEALLOCATE(non_local_stress)

        ok=1
      END SUBROUTINE deallocate_arrays    
!---------------------------------------------------------------------------------------------------------------------------------
!> @brief Computes total mass loss due to corrosion each increment.
!!
!! Sums the mass loss contribution of all elements based on current damage and prints the total if change is significant.
!---------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE compute_mass_loss(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Arguments ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Local variables ---
        REAL(KIND=dp) :: rho, deltam, dmCO
        INTEGER       :: i

        ! --- Local parameters ---
        ! Note: Density (rho) is defined locally here.
        rho    = 7.86e-09 ! [ton/mm^3] Material density
        deltam = 0.0_dp   ! Threshold for mass change to output

        ! --- Main logic ---
        ! Zero corroded mass of current increment before summing.
        mCO = 0.0_dp

        ! Sum mass loss contribution from each element.
        ! (See suggestions regarding index 10)
        DO i = 1, n_elems
            mCO = mCO + (damage_current(i) * rho * elem_properties(i,PROP_IDX_VOLUME))
        END DO

        ! Calculates mass change in increment.
        dmCO = mCO - lmCO

        ! If corrosion not stopped and mass change is significant,
        ! output status and update reference mass.
        IF((global_corrosion_stop_flag /= flag_global_stop) .AND. (dmCO >= deltam)) THEN
            WRITE(*,*) 'Total Corroded Mass:', mCO, ' Time:', timetot
            lmCO = mCO ! Update mass value from previous step.
        END IF

        ok = 1
        RETURN
      END SUBROUTINE compute_mass_loss
!--------------------------------------------------------------------------------------------------------
!> @brief Blocks corrosion in specific regions of the geometry.
!!
!! Called at the start of the analysis to initialize corrosion flags based on element Z positions.
!--------------------------------------------------------------------------------------------------------
      SUBROUTINE block_surface(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Arguments ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Local variables ---
          INTEGER       :: k
          REAL(KIND=dp) :: coordz1, coordz2, coordz3, coordz4, coordx1, coordx2,
     1 coordx3, coordx4, coordy1, coordy2, r2, somaQ

            ! --- Main logic ---
            ! Geometric control to define the non-corroding zone.
            ! (Values and indices taken from Abaqus model)
            coordz1 = 7.4 ! Value from the Abaqus model
            coordz2 = 9.6 ! Value from the Abaqus model
            coordz3 = 0.0 ! Value from the Abaqus model
            coordz4 = 0.4 ! Value from the Abaqus model
      
            coordx1 = -11.0 ! Value from the Abaqus model
            coordx2 = 9.0 ! Value from the Abaqus model
            coordx3 = -0.6 ! Value from the Abaqus model
            coordx4 = -1.0 ! Value from the Abaqus model
      
          coordy1 = -153.0 ! Value from the Abaqus model
          coordy2 = -793.5 ! Value from the Abaqus model
      
      
          r2 = 277.5**2.0

          DO k = 1, n_elems
              ! If element centroid is outside the region of interest...
              somaQ = elem_properties(k,11)**2.0 + elem_properties(k,13)**2.0
              IF ((somaQ.ge.r2).or.(elem_properties(k,13).lt.coordz1).or.
     1 (elem_properties(k,11).gt.coordx1).or.(elem_properties(k,12).gt.coordy1).or.
     2 (elem_properties(k,12).lt.coordy2)) THEN
                  ! ...deactivate all flags relevant to corrosion and damage.
                  property_update_lock(k)                         = key_prop_update_locked
                  elem_properties(k,PROP_IDX_SURFACE_FLAG)        = 0.0_dp
                  element_corrosion_status(k)                     = flag_corrosion_inactive
                  damage_current(k)                               = 0.0_dp
                  strain_deletion_status(k)                       = flag_deletion_by_strain_disabled
              ELSE
                  ! ...otherwise, activate flags to allow corrosion and damage.
                  property_update_lock(k)     = key_prop_update_unlocked
                  element_corrosion_status(k) = flag_corrosion_active
                  damage_current(k)           = 0.0_dp
                  strain_deletion_status(k)   = flag_deletion_by_strain_enabled
              END IF
          END DO

          ok = 1
          RETURN
      END SUBROUTINE block_surface
!--------------------------------------------------------------------------------------------------------------------------
!> @brief Globally stops the corrosion process.
!!
!! Called when a stop condition is reached (e.g., maximum corroded mass). Deactivates corrosion for ALL elements.
!--------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE block_all(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Arguments ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Local variables ---
        INTEGER :: k

        ! --- Main logic ---
        DO k = 1, n_elems
            ! Lock neighbor property updates.
            property_update_lock(k) = key_prop_update_locked
            ! Deactivate corrosion calculation for the element.
            element_corrosion_status(k) = flag_corrosion_inactive
        END DO

        ! Activate global flag to signal complete corrosion stop.
        global_corrosion_stop_flag = flag_global_stop

        ok = 1
        RETURN
      END SUBROUTINE block_all
!-----------------------------------------------------------------------
!> @brief Computes non-local stress for surface elements.
!!
!! Implements a non-local weighted-average stress model,
!! analogous to `compute_nonlocal_property`. Result is stored only if
!! it exceeds material yield stress (`sigth`).
!-----------------------------------------------------------------------
      SUBROUTINE compute_nonlocal_stress(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Arguments ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Local variables ---
          INTEGER       :: i, p, neighbor_element
          REAL(KIND=dp) :: denominator_sum, numerator_sum, weight_p
          REAL(KIND=dp) :: dist_sq, lr_sq

          lr_sq = lr**2.0_dp

          DO i = 1, n_elems
              ! Consider only non-blocked surface elements.
            IF ((elem_properties(i, PROP_IDX_PITTING) /= 0.0_dp) .AND.
     1                (element_corrosion_status(i) /= flag_corrosion_inactive)) THEN             

                  ! Initialize sums with contribution from element itself.
                  denominator_sum = elem_properties(i, PROP_IDX_VOLUME)
                  numerator_sum   = von_mises_stress(i) * elem_properties(i, PROP_IDX_VOLUME)

                  ! Iterate over neighbors to compute their contributions.
                  DO p = 1, counter(i)
                      neighbor_element = influ(i, p)

                      IF (element_corrosion_status(neighbor_element) /= flag_corrosion_inactive) THEN
                          dist_sq  = dist(i, p)**2.0_dp
                          weight_p = (1.0_dp - (dist_sq / lr_sq))**2.0_dp

                          denominator_sum = denominator_sum + (weight_p * elem_properties(neighbor_element, PROP_IDX_VOLUME))
                          numerator_sum = numerator_sum + (weight_p * von_mises_stress(neighbor_element) *
     1                elem_properties(neighbor_element, PROP_IDX_VOLUME))
                      END IF
                  END DO

                  ! Computes final non-local stress value (weighted average).
                  IF (denominator_sum > 1.0e-20_dp) THEN
                      non_local_stress(i) = numerator_sum / denominator_sum
                  ELSE
                      non_local_stress(i) = 0.0_dp
                  END IF

                  ! Apply yield stress threshold.
                  IF (non_local_stress(i) <= sigth) THEN
                      non_local_stress(i) = 0.0_dp
                  END IF
              ELSE
                  ! For internal or inactive elements, non-local stress is zero.
                  non_local_stress(i) = 0.0_dp
              END IF
          END DO

          ok = 1
          RETURN
      END SUBROUTINE compute_nonlocal_stress
