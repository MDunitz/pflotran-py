module Reaction_Sandbox_AWInhibitAcetate_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  PetscInt, parameter :: AWINHIBITACETATE_THRESHOLD_INHIBITION = 1
  PetscInt, parameter :: AWINHIBITACETATE_SMOOTHSTEP_INHIBITION = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibitacetate_type
    
    ! Water activity inhibition parameters
    PetscReal :: aw_threshold
    PetscInt :: inhibition_type

    ! Reaction parameters for acetoclastic methanogenesis: Acetate- + H2O -> CH4(aq) + HCO3- + Tracer
    PetscReal :: rate_constant
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    ! Species indices
    PetscInt :: i_acetate
    PetscInt :: i_h2o
    PetscInt :: i_ch4
    PetscInt :: i_hco3
    PetscInt :: i_tracer

  contains
    procedure, public :: ReadInput  => AWInhibitAcetateRead
    procedure, public :: Setup      => AWInhibitAcetateSetup
    procedure, public :: Evaluate   => AWInhibitAcetateEvaluate
    procedure, public :: Destroy    => AWInhibitAcetateDestroy
  end type reaction_sandbox_awinhibitacetate_type

  public :: AWInhibitAcetateCreate

contains

! ************************************************************************** !

function AWInhibitAcetateCreate()
  !
  ! Allocates AWInhibitAcetate reaction object for water activity inhibited acetoclastic methanogenesis
  !

  implicit none

  class(reaction_sandbox_awinhibitacetate_type), pointer :: AWInhibitAcetateCreate

  allocate(AWInhibitAcetateCreate)
  
  ! Default water activity threshold of 0.5
  AWInhibitAcetateCreate%aw_threshold = 0.5d0
  AWInhibitAcetateCreate%inhibition_type = AWINHIBITACETATE_THRESHOLD_INHIBITION
  
  ! Reaction parameters
  AWInhibitAcetateCreate%rate_constant = UNINITIALIZED_DOUBLE
  AWInhibitAcetateCreate%activation_energy = UNINITIALIZED_DOUBLE
  AWInhibitAcetateCreate%reference_temperature = UNINITIALIZED_DOUBLE
  
  ! Species indices
  AWInhibitAcetateCreate%i_acetate = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_h2o = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_ch4 = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_hco3 = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_tracer = UNINITIALIZED_INTEGER

  nullify(AWInhibitAcetateCreate%next)

end function AWInhibitAcetateCreate

! ************************************************************************** !

subroutine AWInhibitAcetateRead(this,input,option)
  !
  ! Reads input deck for water activity inhibited acetoclastic methanogenesis parameters
  !

  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_awinhibitacetate_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  error_string = 'AWINHIBITACETATE'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       trim(error_string))
    call StringToUpper(word)

    select case(trim(word))
      case('WATER_ACTIVITY_THRESHOLD')
        call InputReadDouble(input,option,this%aw_threshold)
        call InputErrorMsg(input,option,'water_activity_threshold',error_string)
        if (this%aw_threshold < 0.d0 .or. this%aw_threshold > 1.d0) then
          option%io_buffer = 'WATER_ACTIVITY_THRESHOLD must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant',error_string)
        call InputReadAndConvertUnits(input,this%rate_constant,'mol/m^3-sec',&
                        trim(error_string)//',rate_constant',option)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = AWINHIBITACETATE_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = AWINHIBITACETATE_SMOOTHSTEP_INHIBITION
          case default
            error_string = trim(error_string) // ',INHIBITION_TYPE'
            call InputKeywordUnrecognized(input,word,error_string ,option)
        end select

      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,this%activation_energy)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%activation_energy,'j/mol',&
                          trim(error_string)//',activation energy',option)

      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,this%reference_temperature)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%reference_temperature,'C',&
                          trim(error_string)//',reference temperature',option)

      case default
        call InputKeywordUnrecognized(input,word,error_string ,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine AWInhibitAcetateRead

! ************************************************************************** !

subroutine AWInhibitAcetateSetup(this,reaction,option)
  !
  ! Sets up the water activity inhibited acetoclastic methanogenesis reaction
  !

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_awinhibitacetate_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  ! Check that rate constant is provided
  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided for AWInhibitAcetate reaction'
    call PrintErrMsg(option)
  endif

  ! Get species indices for acetoclastic methanogenesis: Acetate- + H2O -> CH4(aq) + HCO3- + Tracer
  word = 'Acetate-'
  this%i_acetate = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'H2O'
  this%i_h2o = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'CH4(aq)'
  this%i_ch4 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'HCO3-'
  this%i_hco3 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'Tracer'
  this%i_tracer = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  if (Initialized(this%activation_energy) .and. &
      UnInitialized(this%reference_temperature)) then
    option%io_buffer = 'A REFERENCE_TEMPERATURE must be provided when an &
      &ACTIVATION_ENERGY is defined in the AWInhibitAcetate Reaction Sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitAcetateSetup

! ************************************************************************** !

subroutine AWInhibitAcetateEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Evaluates the water activity inhibited acetoclastic methanogenesis reaction
  ! Reaction: Acetate- + H2O -> CH4(aq) + HCO3- + Tracer
  ! Inhibited below water activity threshold
  !

  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  class(reaction_sandbox_awinhibitacetate_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: molality_to_molarity
  PetscReal :: water_activity, aw_inhibition, tempreal
  PetscReal :: rate_constant
  PetscReal :: reaction_rate

  ! Concentrations
  PetscReal :: C_acetate, C_h2o, C_ch4, C_hco3, C_tracer

  ! Calculate water volume in liters
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! m^3 -> L

  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3

  ! Get water activity from rt_auxvar
  water_activity = exp(rt_auxvar%ln_act_h2o)

  ! Apply temperature correction if specified
  rate_constant = this%rate_constant
  if (Initialized(this%activation_energy)) then
    rate_constant = rate_constant * Arrhenius(this%activation_energy, &
                                            global_auxvar%temp, &
                                            this%reference_temperature)
  endif

  ! Get aqueous concentrations (convert from molality to molarity)
  C_acetate = rt_auxvar%pri_molal(this%i_acetate) * &
              rt_auxvar%pri_act_coef(this%i_acetate) * molality_to_molarity
  C_h2o = rt_auxvar%pri_molal(this%i_h2o) * &
          rt_auxvar%pri_act_coef(this%i_h2o) * molality_to_molarity
  C_ch4 = rt_auxvar%pri_molal(this%i_ch4) * &
          rt_auxvar%pri_act_coef(this%i_ch4) * molality_to_molarity
  C_hco3 = rt_auxvar%pri_molal(this%i_hco3) * &
           rt_auxvar%pri_act_coef(this%i_hco3) * molality_to_molarity
  C_tracer = rt_auxvar%pri_molal(this%i_tracer) * &
             rt_auxvar%pri_act_coef(this%i_tracer) * molality_to_molarity

  ! Calculate water activity inhibition
  select case(this%inhibition_type)
    case(AWINHIBITACETATE_SMOOTHSTEP_INHIBITION)
      ! Smooth transition around threshold
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.05d0, aw_inhibition, tempreal)
      ! For smoothstep, we want inhibition when aw < threshold, so invert
      aw_inhibition = 1.d0 - aw_inhibition
    case(AWINHIBITACETATE_THRESHOLD_INHIBITION)
      ! Sharp threshold
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  ! Calculate reaction rate: Acetate- + H2O -> CH4(aq) + HCO3- + Tracer
  ! Rate law assumes first order in Acetate- (H2O concentration is essentially constant)
  reaction_rate = rate_constant * C_acetate * aw_inhibition
  
  ! Convert to mol/sec
  reaction_rate = reaction_rate * L_water

  ! Update residuals (negative stoichiometry for reactants, positive for products)
  ! 1 Acetate- consumed
  Residual(this%i_acetate) = Residual(this%i_acetate) + reaction_rate
  
  ! 1 H2O consumed (though this is often negligible compared to solvent water)
  Residual(this%i_h2o) = Residual(this%i_h2o) + reaction_rate
  
  ! 1 CH4(aq) produced
  Residual(this%i_ch4) = Residual(this%i_ch4) - reaction_rate
  
  ! 1 HCO3- produced
  Residual(this%i_hco3) = Residual(this%i_hco3) - reaction_rate
  
  ! 1 Tracer produced
  Residual(this%i_tracer) = Residual(this%i_tracer) - reaction_rate

  if (compute_derivative) then
    option%io_buffer = 'REACTION_SANDBOX AWINHIBITACETATE must be run with &
      &NUMERICAL_JACOBIAN listed in the NUMERICAL_METHODS TRANSPORT &
      &NEWTON_SOLVER block as analytical derivatives are not calculated &
      &in the sandbox evaluate routine.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitAcetateEvaluate

! ************************************************************************** !

subroutine AWInhibitAcetateDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !

  implicit none
  class(reaction_sandbox_awinhibitacetate_type) :: this

  ! Nothing to deallocate in this simple version

end subroutine AWInhibitAcetateDestroy

end module Reaction_Sandbox_AWInhibitAcetate_class