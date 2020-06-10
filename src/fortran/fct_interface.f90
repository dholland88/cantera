! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module fct
interface

   subroutine cantera_error(proc, msg)
     character*(*), intent(in) :: proc
     character*(*), intent(in) :: msg
   end subroutine cantera_error

   integer function phase_nelements(n)
     integer, intent(in) :: n
   end function phase_nelements

    integer function phase_nspecies(n)
        integer, intent(in) :: n
    end function phase_nspecies

    double precision function phase_temperature(n)
        integer, intent(in) :: n
    end function phase_temperature

    integer function phase_settemperature(n, t)
        integer, intent(in) :: n
        double precision, intent(in) :: t
    end function phase_settemperature

    double precision function phase_density(n)
        integer, intent(in) :: n
    end function phase_density

    integer function phase_setdensity(n, rho)
        integer, intent(in) :: n
        double precision, intent(in) :: rho
    end function phase_setdensity

    double precision function phase_molardensity(n)
        integer, intent(in) :: n
    end function phase_molardensity

    double precision function phase_meanmolecularweight(n)
        integer, intent(in) :: n
    end function phase_meanmolecularweight

    integer function phase_elementindex(n, nm)
        integer, intent(in) :: n
        character*(*), intent(in) :: nm
    end function phase_elementindex

    integer function phase_speciesindex(n, nm)
        integer, intent(in) :: n
        character*(*), intent(in) :: nm
    end function phase_speciesindex

    integer function phase_getmolefractions(n, x)
        integer, intent(in) :: n
        double precision, intent(out) :: x(*)
    end function phase_getmolefractions

    double precision function phase_molefraction(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function phase_molefraction

    integer function phase_getmassfractions(n, y)
        integer, intent(in) :: n
        double precision, intent(out) :: y(*)
    end function phase_getmassfractions

    double precision function phase_massfraction(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function phase_massfraction

    integer function phase_setmolefractions(n, x, norm)
        integer, intent(in) :: n
        double precision, intent(in) :: x(*)
        integer, intent(in) :: norm
    end function phase_setmolefractions

    integer function phase_setmolefractionsbyname(n, x)
        integer, intent(in) :: n
        character*(*), intent(in) :: x
    end function phase_setmolefractionsbyname

    integer function phase_setmassfractions(n, y, norm)
        integer, intent(in) :: n
        double precision, intent(in) :: y(*)
        integer, intent(in) :: norm
    end function phase_setmassfractions

    integer function phase_setmassfractionsbyname(n, y)
        integer, intent(in) :: n
        character*(*), intent(in) :: y
    end function phase_setmassfractionsbyname

    integer function phase_getatomicweights(n, atw)
        integer, intent(in) :: n
        double precision, intent(out) :: atw(*)
    end function phase_getatomicweights

    integer function phase_getmolecularweights(n, mw)
        integer, intent(in) :: n
        double precision, intent(out) :: mw(*)
    end function phase_getmolecularweights

    integer function phase_getspeciesname(n, k, nm)
        integer, intent(in) :: n
        integer, intent(in) :: k
        character*(*), intent(out) :: nm
    end function phase_getspeciesname

    integer function phase_getelementname(n, m, nm)
        integer, intent(in) :: n
        integer, intent(in) :: m
        character*(*), intent(out) :: nm
    end function phase_getelementname

    double precision function phase_natoms(n, k, m)
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer, intent(in) :: m
    end function phase_natoms

    integer function th_newfromfile(filename, id)
        character*(*), intent(in) :: filename
        character*(*), intent(in) :: id
    end function th_newfromfile

    integer function newthermofromxml(mxml)
        integer, intent(in) :: mxml
    end function newthermofromxml

    integer function th_geteostype(n, buf)
        integer, intent(in) :: n
        character*(*), intent(out) :: buf
    end function th_geteostype

    double precision function th_enthalpy_mole(n)
        integer, intent(in) :: n
    end function th_enthalpy_mole

    double precision function th_intenergy_mole(n)
        integer, intent(in) :: n
    end function th_intenergy_mole

    double precision function th_entropy_mole(n)
        integer, intent(in) :: n
    end function th_entropy_mole

    double precision function th_gibbs_mole(n)
        integer, intent(in) :: n
    end function th_gibbs_mole

    double precision function th_cp_mole(n)
        integer, intent(in) :: n
    end function th_cp_mole

    double precision function th_cv_mole(n)
        integer, intent(in) :: n
    end function th_cv_mole

    double precision function th_pressure(n)
        integer, intent(in) :: n
    end function th_pressure

    double precision function th_enthalpy_mass(n)
        integer, intent(in) :: n
    end function th_enthalpy_mass

    double precision function th_intEnergy_mass(n)
        integer, intent(in) :: n
    end function th_intEnergy_mass

    double precision function th_entropy_mass(n)
        integer, intent(in) :: n
    end function th_entropy_mass

    double precision function th_gibbs_mass(n)
        integer, intent(in) :: n
    end function th_gibbs_mass

    double precision function th_cp_mass(n)
        integer, intent(in) :: n
    end function th_cp_mass

    double precision function th_cv_mass(n)
        integer, intent(in) :: n
    end function th_cv_mass

    integer function th_chempotentials(n, murt)
        integer, intent(in) :: n
        double precision, intent(out) :: murt(*)
    end function th_chempotentials

    integer function th_setpressure(n, p)
        integer, intent(in) :: n
        double precision, intent(in) :: p
    end function th_setpressure

    integer function th_set_hp(n, v1, v2)
        integer, intent(in) :: n
        double precision, intent(in) :: v1
        double precision, intent(in) :: v2
    end function th_set_hp

    integer function th_set_uv(n, v1, v2)
        integer, intent(in) :: n
        double precision, intent(in) :: v1
        double precision, intent(in) :: v2
    end function th_set_uv

    integer function th_set_sv(n, v1, v2)
        integer, intent(in) :: n
        double precision, intent(in) :: v1
        double precision, intent(in) :: v2
    end function th_set_sv

    integer function th_set_sp(n, v1, v2)
        integer, intent(in) :: n
        double precision, intent(in) :: v1
        double precision, intent(in) :: v2
    end function th_set_sp

    integer function th_equil(n, XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level)
        integer, intent(in) :: n
        character*(*), intent(in) :: XY
        character*(*), intent(in) :: solver
        double precision, intent(in) :: rtol
        integer, intent(in) :: max_steps
        integer, intent(in) :: max_iter
        integer, intent(in) :: estimate_equil
        integer, intent(in) :: log_level
    end function th_equil

    double precision function th_refpressure(n)
        integer, intent(in) :: n
    end function th_refpressure

    double precision function th_mintemp(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function th_mintemp

    double precision function th_maxtemp(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function th_maxtemp

    integer function th_getenthalpies_rt(n, h_rt)
        integer, intent(in) :: n
        double precision, intent(out) :: h_rt(*)
    end function th_getenthalpies_rt

    integer function th_getentropies_r(n, s_r)
        integer, intent(in) :: n
        double precision, intent(out) :: s_r(*)
    end function th_getentropies_r

    integer function th_getcp_r(n, lenm, cp_r)
        integer, intent(in) :: n
        integer, intent(out) :: lenm
        double precision, intent(out) :: cp_r(*)
    end function th_getcp_r

    integer function th_getpartialmolarintenergies_r(n, ie)
        integer, intent(in) :: n
        double precision, intent(out) :: ie(*)
    end function th_getpartialmolarintenergies_r

    integer function kin_newfromfile(filename, id, reactingPhase, neighbor1, neighbor2, neighbor3, neighbor4)
        character*(*), intent(in) :: filename
        character*(*), intent(in) :: id
        integer, intent(in) :: reactingPhase
        integer, intent(in) :: neighbor1
        integer, intent(in) :: neighbor2
        integer, intent(in) :: neighbor3
        integer, intent(in) :: neighbor4
    end function kin_newfromfile

    integer function newkineticsfromxml(mxml, iphase, neighbor1, neighbor2, neighbor3, neighbor4)
        integer, intent(in) :: mxml
        integer, intent(in) :: iphase
        integer, intent(in) :: neighbor1
        integer, intent(in) :: neighbor2
        integer, intent(in) :: neighbor3
        integer, intent(in) :: neighbor4
    end function newkineticsfromxml

    integer function kin_gettype(n, buf)
        integer, intent(in) :: n
        character*(*), intent(out) :: buf
    end function kin_gettype

    integer function kin_start(n, p)
        integer, intent(in) :: n
        integer, intent(in) :: p
    end function kin_start

    integer function kin_speciesindex(n, nm, ph)
        integer, intent(in) :: n
        character*(*), intent(in) :: nm
        character*(*), intent(in) :: ph
    end function kin_speciesindex

    integer function kin_ntotalspecies(n)
        integer, intent(in) :: n
    end function kin_ntotalspecies

    integer function kin_nreactions(n)
        integer, intent(in) :: n
    end function kin_nreactions

    integer function kin_nphases(n)
        integer, intent(in) :: n
    end function kin_nphases

    integer function kin_phaseIndex(n, phase)
        integer, intent(in) :: n
        character*(*), intent(in) :: phase
    end function kin_phaseindex

    double precision function kin_reactantstoichcoeff(n, k, i)
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer, intent(in) :: i
    end function kin_reactantstoichcoeff

    double precision function kin_productstoichcoeff(n, k, i)
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer, intent(in) :: i
    end function kin_productstoichcoeff

    integer function kin_reactiontype(n, i)
        integer, intent(in) :: n
        integer, intent(in) :: i
    end function kin_reactiontype

    integer function kin_getfwdratesofprogress(n, fwdROP)
        integer, intent(in) :: n
        double precision, intent(out) :: fwdROP(*)
    end function kin_getfwdratesofprogress

    integer function kin_getrevratesofprogress(n, revROP)
        integer, intent(in) :: n
        double precision, intent(out) :: revROP(*)
    end function kin_getrevratesofprogress

    integer function kin_isreversible(n, i)
        integer, intent(in) :: n
        integer, intent(in) :: i
    end function kin_isreversible

    integer function kin_getnetratesofprogress(n, netROP)
        integer, intent(in) :: n
        double precision, intent(out) :: netROP(*)
    end function kin_getnetratesofprogress

    integer function kin_getcreationrates(n, cdot)
        integer, intent(in) :: n
        double precision, intent(out) :: cdot(*)
    end function kin_getcreationrates

    integer function kin_getdestructionrates(n, ddot)
        integer, intent(in) :: n
        double precision, intent(out) :: ddot(*)
    end function kin_getdestructionrates

    integer function kin_getnetproductionrates(n, wdot)
        integer, intent(in) :: n
        double precision, intent(out) :: wdot(*)
    end function kin_getnetproductionrates

    double precision function kin_multiplier(n, i)
        integer, intent(in) :: n
        integer, intent(in) :: i
    end function kin_multiplier

    integer function kin_getequilibriumconstants(n, kc)
        integer, intent(in) :: n
        double precision, intent(out) :: kc(*)
    end function kin_getequilibriumconstants

    integer function kin_getreactionstring(n, i, buf)
        integer, intent(in) :: n
        integer, intent(in) :: i
        character*(*), intent(out) :: buf
    end function kin_getreactionstring

    integer function kin_setmultiplier(n, i, v)
        integer, intent(in) :: n
        integer, intent(in) :: i
        double precision, intent(out) :: v
    end function kin_setmultiplier

    integer function kin_advancecoverages(n, tstep)
        integer, intent(in) :: n
        double precision, intent(in) :: tstep
    end function kin_advancecoverages

    integer function newtransport(model, ith, loglevel)
        character*(*), intent(in) :: model
        integer, intent(in) :: ith
        integer, intent(in) :: loglevel
    end function newtransport

    integer function trans_newdefault(ith, loglevel)
        integer, intent(in) :: ith
        integer, intent(in) :: loglevel
    end function trans_newdefault

    double precision function trans_viscosity(n)
        integer, intent(in) :: n
    end function trans_viscosity

    double precision function trans_electricalConductivity(n)
        integer, intent(in) :: n
    end function trans_electricalConductivity

    double precision function trans_thermalConductivity(n)
        integer, intent(in) :: n
    end function trans_thermalConductivity

    integer function trans_getThermalDiffCoeffs(n, dt)
        integer, intent(in) :: n
        double precision, intent(out) :: dt(*)
    end function trans_getThermalDiffCoeffs

    integer function trans_getMixDiffCoeffs(n, d)
        integer, intent(in) :: n
        double precision, intent(out) :: d(*)
    end function trans_getMixDiffCoeffs

    integer function trans_getMixDiffCoeffsMass(n, d)
        integer, intent(in) :: n
        double precision, intent(out) :: d(*)
    end function trans_getMixDiffCoeffsMass

    integer function trans_getMixDiffCoeffsMole(n, d)
        integer, intent(in) :: n
        double precision, intent(out) :: d(*)
    end function trans_getMixDiffCoeffsMole

    integer function trans_getBinDiffCoeffs(n, ld, d)
        integer, intent(in) :: n
        integer, intent(in) :: ld
        double precision, intent(out) :: d(*)
    end function trans_getBinDiffCoeffs

    integer function trans_getMultiDiffCoeffs(n, ld, d)
        integer, intent(in) :: n
        integer, intent(in) :: ld
        double precision, intent(out) :: d(*)
    end function trans_getMultiDiffCoeffs

    integer function trans_setParameters(n, type, k, d)
        integer, intent(in) :: n
        integer, intent(in) :: type
        integer, intent(in) :: k
        double precision, intent(in) :: d(*)
    end function trans_setParameters

    integer function ctphase_report(nth, buf, show_thermo)
        integer, intent(in) :: nth
        character*(*), intent(out) :: buf
        integer, intent(in) :: show_thermo
    end function ctphase_report

    integer function ctgetCanteraError(buf)
        character*(*), intent(out) :: buf
    end function ctgetCanteraError

    integer function ctaddCanteraDirectory(buflen, buf)
        integer, intent(in) :: buflen
        character*(*), intent(in) :: buf
    end function ctaddCanteraDirectory

    integer function ctbuildSolutionFromXML(src, ixml, id, ith, ikin)
        character*(*), intent(in) :: src
        integer, intent(in) :: ixml
        character*(*), intent(in) :: id
        integer, intent(in) :: ith
        integer, intent(in) :: ikin
    end function ctbuildSolutionFromXML

    integer function mphase_addmultiphase()
    end function mphase_addmultiphase

    integer function mphase_addphase(n, p, x)
        integer, intent(in) :: n
        integer, intent(in) :: p
        double precision, intent(in) :: x
    end function mphase_addphase

    integer function mphase_elementindex(n, nm)
        integer, intent(in) :: n
        character*(*), intent(in) :: nm
    end function mphase_elementindex

    integer function mphase_elementname(n, k, nm)
        integer, intent(in) :: n
        integer, intent(in) :: k
        character*(*), intent(out) :: nm
    end function mphase_elementname

    integer function mphase_speciesname(n, k, nm)
        integer, intent(in) :: n
        integer, intent(in) :: k
        character*(*), intent(out) :: nm
    end function mphase_speciesname

    double precision function mphase_natoms(n, k, m)
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer, intent(in) :: m
    end function mphase_natoms

    integer function mphase_getmolefractions(n, x)
        integer, intent(in) :: n
        double precision, intent(out) :: x(*)
    end function mphase_getmolefractions

    integer function mphase_init(n)
        integer, intent(in) :: n
    end function mphase_init

    integer function mphase_phasename(n, p, nm)
        integer, intent(in) :: n
        integer, intent(in) :: p
        character*(*), intent(out) :: nm
    end function mphase_phasename

    integer function mphase_phaseindex(n, nm)
        integer, intent(in) :: n
        character*(*), intent(in) :: nm
    end function mphase_phaseindex

    double precision function mphase_phasemoles(n, p)
        integer, intent(in) :: n
        integer, intent(in) :: p
    end function mphase_phaseMoles

    integer function mphase_setphasemoles(n, p, x)
        integer, intent(in) :: n
        integer, intent(in) :: p
        double precision, intent(in) :: x
    end function mphase_setphasemoles

    integer function mphase_speciesmoles(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function mphase_speciesmoles

    integer function mphase_speciesindex(n, k, p)
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer, intent(in) :: p
    end function mphase_speciesindex

    integer function mphase_speciesindexbyname(n, nm, pname)
        integer, intent(in) :: n
        character*(*), intent(in) :: nm
        character*(*), intent(in) :: pname
    end function mphase_speciesindexbyname

    double precision function mphase_mintemp(n)
        integer, intent(in) :: n
    end function mphase_mintemp

    double precision function mphase_maxtemp(n)
        integer, intent(in) :: n
    end function mphase_maxtemp

    double precision function mphase_charge(n)
        integer, intent(in) :: n
    end function mphase_charge

    double precision function mphase_phasecharge(n, p)
        integer, intent(in) :: n
        integer, intent(in) :: p
    end function mphase_phasecharge

    double precision function mphase_elementmoles(n, m)
        integer, intent(in) :: n
        integer, intent(in) :: m
    end function mphase_elementmoles

    integer function mphase_getchempotentials(n, mu)
        integer, intent(in) :: n
        double precision, intent(out) :: mu(*)
    end function mphase_getchempotentials

    double precision function mphase_temperature(n)
        integer, intent(in) :: n
    end function mphase_temperature

    integer function mphase_equil(n, XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level)
        integer, intent(in) :: n
        character*(*), intent(in) :: XY
        character*(*), intent(in) :: solver
        double precision, intent(in) :: rtol
        integer, intent(in) :: max_steps
        integer, intent(in) :: max_iter
        integer, intent(in) :: estimate_equil
        integer, intent(in) :: log_level
    end function mphase_equil

    integer function mphase_settemperature(n, t)
        integer, intent(in) :: n
        double precision, intent(in) :: t
    end function mphase_settemperature

    integer function mphase_set_tp(n, v1, v2)
        integer, intent(in) :: n
        double precision, intent(in) :: v1
        double precision, intent(in) :: v2
    end function mphase_set_tp

    integer function mphase_set_tpmoles(n, v1, v2, x)
        integer, intent(in) :: n
        double precision, intent(in) :: v1
        double precision, intent(in) :: v2
        double precision, intent(in) :: x(*)
    end function mphase_set_tpmoles

    double precision function mphase_volume(n)
        integer, intent(in) :: n
    end function mphase_volume

    double precision function mphase_pressure(n)
        integer, intent(in) :: n
    end function mphase_pressure

    integer function mphase_setpressure(n, p)
        integer, intent(in) :: n
        double precision, intent(in) :: p
    end function mphase_setpressure

    double precision function mphase_enthalpy(n)
        integer, intent(in) :: n
    end function mphase_enthalpy

    double precision function mphase_intenergy(n)
        integer, intent(in) :: n
    end function mphase_intenergy

    double precision function mphase_entropy(n)
        integer, intent(in) :: n
    end function mphase_entropy

    double precision function mphase_gibbs(n)
        integer, intent(in) :: n
    end function mphase_gibbs

    double precision function mphase_cp(n)
        integer, intent(in) :: n
    end function mphase_cp

    integer function mphase_speciesphaseindex(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function mphase_speciesphaseindex

    double precision function mphase_molefraction(n, k)
        integer, intent(in) :: n
        integer, intent(in) :: k
    end function mphase_molefraction

    integer function mphase_setphasemolefractions(n, p, x)
        integer, intent(in) :: n
        integer, intent(in) :: p
        double precision, intent(in) :: x(*)
    end function mphase_setphasemolefractions

    integer function mphase_setmoles(n, x)
        integer, intent(in) :: n
        double precision, intent(in) :: x(*)
    end function mphase_setmoles

    integer function mphase_setmolesbyname(n, x)
        integer, intent(in) :: n
        character*(*), intent(in) :: x
    end function mphase_setmolesbyname

    double precision function mphase_getmoles(n, m)
    integer, intent(in) :: n
    double precision, intent(out) :: m(*)
end function mphase_getmoles

    integer function mphase_addspeciesmoles(n, k, x)
        integer, intent(in) :: n
        integer, intent(in) :: k
        double precision, intent(in) :: x
    end function mphase_addspeciesmoles

    integer function mphase_getelemabundances(n, v1)
        integer, intent(in) :: n
        double precision, intent(out) :: v1(*)
    end function mphase_getelemabundances

    integer function mphase_uploadmolefractions(n)
        integer, intent(in) :: n
    end function mphase_uploadmolefractions

    integer function mphase_report(n)
        integer, intent(in) :: n
    end function mphase_report

end interface
end module fct
