# Scaled Isoplane Kernel

The Scaled Isoplane Dose Calculation Kernel is a simplified pencil beam kernel algorithm.
It only uses 1 single pencil beam kernel for a given field size (typically 10 x 10 cm²).
It then scales the dose by using the measured percentage depth dose combined with a few correction factors.
More details can be found in the original paper.

## Data

Example kernel data is located in `examples/sample-data/dose-kernel/scaled-isoplane-kernel.json`.
This data specification is as follows:

* Energy: Beam energy (*e.g*. `"6MV"`)
* Kernel
  * SSD: Source Surface Distance in mm
  * Radius: Array of off-axis distance (radius) in mm
  * Kernel Value: Array of kernel values at specified radius
* Depth Dose
  * SSD: Source Surface Distance in mm
  * Field Size: Field size of measured depth dose
  * Depth: Array of depths below surface in mm
  * Dose: Dose at specified depths, normalised to 1

<!-- | File | Info |
|:-----|:-----| 
| - Energy | Beam energy (*e.g*. `"6MV"`) |
| - Kernel | |
|   - SSD | Source Surface Distance in mm |
|   - Radius | Array of off-axis distance (radius) in mm |
|   - Kernel Value | Array of kernel values at specified radius | -->

