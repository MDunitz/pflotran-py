module Reaction_Sandbox_AWInhibit_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  PetscInt, parameter :: AWINHIBIT_THRESHOLD_INHIBITION = 1
  PetscInt, parameter :: AWINHIBIT_SMOOTHSTEP_INHIBITION = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_awinhibit_type
    
    ! Water activity inhibition parameters
    PetscReal :: aw_threshold
    PetscInt :: inhibition_type

    ! Reaction parameters for methanogenesis: 4 H2(aq) + HCO3- + H+ -> CH4(aq) + 3 H2O
    PetscReal :: rate_constant
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    ! Species indices
    PetscInt :: i_h2
    PetscInt :: i_hco3
    PetscInt :: i_h
    PetscInt :: i_ch4

  contains
    procedure, public :: ReadInput  => AWInhibitRead
    procedure, public :: Setup      => AWInhibitSetup
    procedure, public :: Evaluate   => AWInhibitEvaluate
    procedure, public :: Destroy    => AWInhibitDestroy
  end type reaction_sandbox_awinhibit_type

  public :: AWInhibitCreate

contains

! ************************************************************************** !

function AWInhibitCreate()
  !
  ! Allocates AWInhibit reaction object for water activity inhibited methanogenesis
  !

  implicit none

  class(reaction_sandbox_awinhibit_type), pointer :: AWInhibitCreate

  allocate(AWInhibitCreate)
  
  ! Default water activity threshold of 0.5
  AWInhibitCreate%aw_threshold = 0.5d0
  AWInhibitCreate%inhibition_type = AWINHIBIT_THRESHOLD_INHIBITION
  
  ! Reaction parameters
  AWInhibitCreate%rate_constant = UNINITIALIZED_DOUBLE
  AWInhibitCreate%activation_energy = UNINITIALIZED_DOUBLE
  AWInhibitCreate%reference_temperature = UNINITIALIZED_DOUBLE
  
  ! Species indices
  AWInhibitCreate%i_h2 = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_hco3 = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_h = UNINITIALIZED_INTEGER
  AWInhibitCreate%i_ch4 = UNINITIALIZED_INTEGER

  nullify(AWInhibitCreate%next)

end function AWInhibitCreate

! ************************************************************************** !

subroutine AWInhibitRead(this,input,option)
  !
  ! Reads input deck for water activity inhibited methanogenesis parameters
  !

  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_awinhibit_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  error_string = 'CHEMISTRY,RXN_SANDBOX,AWINHIBIT'

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
            this%inhibition_type = AWINHIBIT_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = AWINHIBIT_SMOOTHSTEP_INHIBITION
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

end subroutine AWInhibitRead

! ************************************************************************** !

subroutine AWInhibitSetup(this,reaction,option)
  !
  ! Sets up the water activity inhibited methanogenesis reaction
  !

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_awinhibit_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  ! Check that rate constant is provided
  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided for AWInhibit reaction'
    call PrintErrMsg(option)
  endif

  ! Get species indices for methanogenesis reaction: 4 H2(aq) + HCO3- + H+ -> CH4(aq) + 3 H2O
  word = 'H2(aq)'
  this%i_h2 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'HCO3-'
  this%i_hco3 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'H+'
  this%i_h = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  word = 'CH4(aq)'
  this%i_ch4 = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  if (Initialized(this%activation_energy) .and. &
      UnInitialized(this%reference_temperature)) then
    option%io_buffer = 'A REFERENCE_TEMPERATURE must be provided when an &
      &ACTIVATION_ENERGY is defined in the AWInhibit Reaction Sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitSetup

! ************************************************************************** !

subroutine AWInhibitEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Evaluates the water activity inhibited methanogenesis reaction
  ! Reaction: 4 H2(aq) + HCO3- + H+ -> CH4(aq) + 3 H2O
  ! Inhibited below water activity threshold
  !

  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  class(reaction_sandbox_awinhibit_type) :: this
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
  PetscReal :: C_h2, C_hco3, C_h, C_ch4

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
  C_h2 = rt_auxvar%pri_molal(this%i_h2) * &
         rt_auxvar%pri_act_coef(this%i_h2) * molality_to_molarity
  C_hco3 = rt_auxvar%pri_molal(this%i_hco3) * &
           rt_auxvar%pri_act_coef(this%i_hco3) * molality_to_molarity
  C_h = rt_auxvar%pri_molal(this%i_h) * &
        rt_auxvar%pri_act_coef(this%i_h) * molality_to_molarity
  C_ch4 = rt_auxvar%pri_molal(this%i_ch4) * &
          rt_auxvar%pri_act_coef(this%i_ch4) * molality_to_molarity

  ! Calculate water activity inhibition
  select case(this%inhibition_type)
    case(AWINHIBIT_SMOOTHSTEP_INHIBITION)
      ! Smooth transition around threshold
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.05d0, aw_inhibition, tempreal)
      ! For smoothstep, we want inhibition when aw < threshold, so invert
      aw_inhibition = 1.d0 - aw_inhibition
    case(AWINHIBIT_THRESHOLD_INHIBITION)
      ! Sharp threshold
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  ! Calculate reaction rate: 4 H2(aq) + HCO3- + H+ -> CH4(aq) + 3 H2O
  ! Rate law assumes first order in each reactant
  reaction_rate = rate_constant * (C_h2**4) * C_hco3 * C_h * aw_inhibition
  
  ! Apply water activity inhibition
  reaction_rate = reaction_rate * L_water  ! Convert to mol/sec

  ! Update residuals (negative stoichiometry for reactants, positive for products)
  ! 4 H2(aq) consumed
  Residual(this%i_h2) = Residual(this%i_h2) + 4.d0 * reaction_rate
  
  ! 1 HCO3- consumed
  Residual(this%i_hco3) = Residual(this%i_hco3) + reaction_rate
  
  ! 1 H+ consumed
  Residual(this%i_h) = Residual(this%i_h) + reaction_rate
  
  ! 1 CH4(aq) produced
  Residual(this%i_ch4) = Residual(this%i_ch4) - reaction_rate

  ! Note: H2O is not tracked as it's the solvent

  if (compute_derivative) then
    option%io_buffer = 'REACTION_SANDBOX AWINHIBIT must be run with &
      &NUMERICAL_JACOBIAN listed in the NUMERICAL_METHODS TRANSPORT &
      &NEWTON_SOLVER block as analytical derivatives are not calculated &
      &in the sandbox evaluate routine.'
    call PrintErrMsg(option)
  endif

end subroutine AWInhibitEvaluate

! ************************************************************************** !

subroutine AWInhibitDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !

  implicit none
  class(reaction_sandbox_awinhibit_type) :: this

  ! Nothing to deallocate in this simple version

end subroutine AWInhibitDestroy

end module Reaction_Sandbox_AWInhibit_class
