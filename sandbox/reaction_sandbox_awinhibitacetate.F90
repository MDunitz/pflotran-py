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
  PetscInt, parameter :: AWINHIBITACETATE_ONE_MINUS_AW_INHIBITION = 3

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibitacetate_type

    PetscReal :: aw_threshold
    PetscInt :: inhibition_type
    PetscReal :: fixed_water_activity

    ! Network-matching Monod kinetics for acetoclastic methanogenesis:
    !   Acetate- + H2O -> CH4(aq) + HCO3- + Tracer
    PetscReal :: rate_constant
    PetscReal :: half_saturation_acetate
    PetscReal :: threshold_acetate
    PetscReal :: o2_inhibition
    PetscReal :: fe_inhibition
    PetscReal :: h_inhibition_above
    PetscReal :: h_inhibition_below
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    PetscInt :: i_acetate
    PetscInt :: i_ch4
    PetscInt :: i_hco3
    PetscInt :: i_tracer
    PetscInt :: i_o2
    PetscInt :: i_fe
    PetscInt :: i_h

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

  implicit none

  class(reaction_sandbox_awinhibitacetate_type), pointer :: AWInhibitAcetateCreate

  allocate(AWInhibitAcetateCreate)

  AWInhibitAcetateCreate%aw_threshold = 0.95d0
  AWInhibitAcetateCreate%inhibition_type = AWINHIBITACETATE_SMOOTHSTEP_INHIBITION
  AWInhibitAcetateCreate%fixed_water_activity = UNINITIALIZED_DOUBLE

  AWInhibitAcetateCreate%rate_constant = UNINITIALIZED_DOUBLE
  AWInhibitAcetateCreate%half_saturation_acetate = 4.0d-2
  AWInhibitAcetateCreate%threshold_acetate = 1.1d-15
  AWInhibitAcetateCreate%o2_inhibition = 1.0d-6
  AWInhibitAcetateCreate%fe_inhibition = 1.0d-9
  AWInhibitAcetateCreate%h_inhibition_above = 2.88d-5
  AWInhibitAcetateCreate%h_inhibition_below = 2.88d-7
  AWInhibitAcetateCreate%activation_energy = UNINITIALIZED_DOUBLE
  AWInhibitAcetateCreate%reference_temperature = UNINITIALIZED_DOUBLE

  AWInhibitAcetateCreate%i_acetate = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_ch4 = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_hco3 = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_tracer = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_o2 = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_fe = UNINITIALIZED_INTEGER
  AWInhibitAcetateCreate%i_h = UNINITIALIZED_INTEGER

  nullify(AWInhibitAcetateCreate%next)

end function AWInhibitAcetateCreate

! ************************************************************************** !

subroutine AWInhibitAcetateRead(this,input,option)

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
    call InputErrorMsg(input,option,'keyword',trim(error_string))
    call StringToUpper(word)

    select case(trim(word))
      case('WATER_ACTIVITY_THRESHOLD')
        call InputReadDouble(input,option,this%aw_threshold)
        call InputErrorMsg(input,option,'water_activity_threshold',error_string)
        if (this%aw_threshold < 0.d0 .or. this%aw_threshold > 1.d0) then
          option%io_buffer = 'WATER_ACTIVITY_THRESHOLD must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('FIXED_WATER_ACTIVITY')
        call InputReadDouble(input,option,this%fixed_water_activity)
        call InputErrorMsg(input,option,'fixed_water_activity',error_string)
        if (this%fixed_water_activity < 0.d0 .or. &
            this%fixed_water_activity > 1.d0) then
          option%io_buffer = 'FIXED_WATER_ACTIVITY must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant',error_string)
        call InputReadAndConvertUnits(input,this%rate_constant,'mol/L-sec',&
                        trim(error_string)//',rate_constant',option)

      case('HALF_SATURATION_ACETATE')
        call InputReadDouble(input,option,this%half_saturation_acetate)
        call InputErrorMsg(input,option,'half_saturation_acetate',error_string)

      case('THRESHOLD_ACETATE')
        call InputReadDouble(input,option,this%threshold_acetate)
        call InputErrorMsg(input,option,'threshold_acetate',error_string)

      case('O2_INHIBITION')
        call InputReadDouble(input,option,this%o2_inhibition)
        call InputErrorMsg(input,option,'o2_inhibition',error_string)

      case('FE_INHIBITION')
        call InputReadDouble(input,option,this%fe_inhibition)
        call InputErrorMsg(input,option,'fe_inhibition',error_string)

      case('H_INHIBITION_ABOVE')
        call InputReadDouble(input,option,this%h_inhibition_above)
        call InputErrorMsg(input,option,'h_inhibition_above',error_string)

      case('H_INHIBITION_BELOW')
        call InputReadDouble(input,option,this%h_inhibition_below)
        call InputErrorMsg(input,option,'h_inhibition_below',error_string)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = AWINHIBITACETATE_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = AWINHIBITACETATE_SMOOTHSTEP_INHIBITION
          case('ONE_MINUS_AW')
            this%inhibition_type = AWINHIBITACETATE_ONE_MINUS_AW_INHIBITION
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

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_awinhibitacetate_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided for AWInhibitAcetate'
    call PrintErrMsg(option)
  endif

  word = 'Acetate-'
  this%i_acetate = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'CH4(aq)'
  this%i_ch4 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'HCO3-'
  this%i_hco3 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Tracer'
  this%i_tracer = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%i_o2 = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Fe+++'
  this%i_fe = ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'H+'
  this%i_h = ReactionAuxGetPriSpecIDFromName(word,reaction,option)

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
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water, molality_to_molarity
  PetscReal :: water_activity, aw_inhibition, tempreal
  PetscReal :: rate_constant, reaction_rate
  PetscReal :: C_acetate, C_o2, C_fe, C_h
  PetscReal :: monod_ac, inhib_o2, inhib_fe, inhib_h_above, inhib_h_below

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  if (Initialized(this%fixed_water_activity)) then
    water_activity = this%fixed_water_activity
  else
    water_activity = exp(rt_auxvar%ln_act_h2o)
  endif

  rate_constant = this%rate_constant
  if (Initialized(this%activation_energy)) then
    rate_constant = rate_constant * Arrhenius(this%activation_energy, &
                                            global_auxvar%temp, &
                                            this%reference_temperature)
  endif

  C_acetate = rt_auxvar%pri_molal(this%i_acetate) * &
              rt_auxvar%pri_act_coef(this%i_acetate) * molality_to_molarity
  C_o2 = rt_auxvar%pri_molal(this%i_o2) * &
         rt_auxvar%pri_act_coef(this%i_o2) * molality_to_molarity
  C_fe = rt_auxvar%pri_molal(this%i_fe) * &
         rt_auxvar%pri_act_coef(this%i_fe) * molality_to_molarity
  C_h = rt_auxvar%pri_molal(this%i_h) * &
        rt_auxvar%pri_act_coef(this%i_h) * molality_to_molarity

  if (C_acetate < this%threshold_acetate) return

  monod_ac = C_acetate / (this%half_saturation_acetate + C_acetate)
  inhib_o2 = this%o2_inhibition / (this%o2_inhibition + C_o2)
  inhib_fe = this%fe_inhibition / (this%fe_inhibition + C_fe)
  inhib_h_above = this%h_inhibition_above / (this%h_inhibition_above + C_h)
  ! Monod inhibition BELOW threshold: C / (Ki + C)
  inhib_h_below = C_h / (this%h_inhibition_below + C_h)

  select case(this%inhibition_type)
    case(AWINHIBITACETATE_SMOOTHSTEP_INHIBITION)
      ! Positive threshold => INHIBIT_BELOW: factor -> 1 as a_w rises (wet on).
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.20d0, aw_inhibition, tempreal)
    case(AWINHIBITACETATE_ONE_MINUS_AW_INHIBITION)
      if (this%aw_threshold >= 1.d0) then
        aw_inhibition = 1.d0
      else if (water_activity <= this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = (water_activity - this%aw_threshold) / &
                        (1.d0 - this%aw_threshold)
      endif
    case(AWINHIBITACETATE_THRESHOLD_INHIBITION)
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  reaction_rate = rate_constant * monod_ac * inhib_o2 * inhib_fe * &
                  inhib_h_above * inhib_h_below * aw_inhibition
  reaction_rate = reaction_rate * L_water

  Residual(this%i_acetate) = Residual(this%i_acetate) + reaction_rate
  Residual(this%i_ch4) = Residual(this%i_ch4) - reaction_rate
  Residual(this%i_hco3) = Residual(this%i_hco3) - reaction_rate
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

  implicit none
  class(reaction_sandbox_awinhibitacetate_type) :: this

end subroutine AWInhibitAcetateDestroy

end module Reaction_Sandbox_AWInhibitAcetate_class
