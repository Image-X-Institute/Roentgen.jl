#
#    Coordinate System Transformations
#
# Constains functions which convert between standard radiotherapy coordinate
# systems, such as IEC Fixed, IEC Gantry and IEC BeamLimitingDevice (BLD)
#

using CoordinateTransformations, Rotations, StaticArrays

export fixed_to_gantry, gantry_to_fixed
export gantry_to_bld, bld_to_gantry
export fixed_to_bld, bld_to_fixed
export patient_to_fixed

#--- Fixed -> Gantry ----------------------------------------------------------
"""
    fixed_to_gantry(ϕg)

Convert from IEC fixed to IEC gantry coordinates
""" 
fixed_to_gantry(ϕg) = LinearMap(RotY(-ϕg))

#--- Fixed <- Gantry ----------------------------------------------------------
"""
    gantry_to_fixed(ϕg)

Convert from IEC gantry to IEC fixed coordinates
"""
gantry_to_fixed(ϕg) = LinearMap(RotY(ϕg))

#--- Gantry -> BLD ------------------------------------------------------------
"""
    gantry_to_bld(θb, SAD)

Convert from IEC gantry to IEC BLD coordinates
"""
gantry_to_bld(θb, SAD::T) where T<:AbstractFloat = LinearMap(RotZ(-θb)) ∘ Translation(zero(T), zero(T), -SAD)

#--- Gantry <- BLD ------------------------------------------------------------
"""
    bld_to_gantry(θb, SAD)

Convert from IEC BLD to IEC gantry coordinates
"""
bld_to_gantry(θb, SAD::T) where T<:AbstractFloat = Translation(zero(T), zero(T), SAD) ∘ LinearMap(RotZ(θb))

#--- Fixed -> BLD -------------------------------------------------------------
"""
    fixed_to_bld(ϕg, θb, SAD)

Convert from IEC fixed to IEC BLD coordinates
"""
fixed_to_bld(ϕg, θb, SAD) = gantry_to_bld(θb, SAD) ∘ fixed_to_gantry(ϕg)

#--- Fixed <- BLD -------------------------------------------------------------
"""
    bld_to_fixed(ϕg, θb, SAD)

Convert from IEC BLD to IEC fixed coordinates
"""
bld_to_fixed(ϕg, θb, SAD) = gantry_to_fixed(ϕg) ∘ bld_to_gantry(θb, SAD)

#--- Patient -> Fixed ---------------------------------------------------------
"""
    patient_to_fixed(isocenter)

Convert from the patient-based coordinate system to IEC Fixed.

Isocenter is in the patient-based coordinate system, as specified in the DICOM
RP file.
"""
patient_to_fixed(isocenter) = LinearMap(@SMatrix [1  0  0
                                                  0  0  1
                                                  0 -1  0]) ∘ Translation(-isocenter)
