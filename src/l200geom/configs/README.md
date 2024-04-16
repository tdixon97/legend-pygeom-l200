# HPGe string configuration

## Global HPGe string configuration

- `hpge_string` → HPGe string number
  - `radius_in_mm` → radial distance from the center of the cryostat to the
    string
  - `angle_in_deg` → azimutal position of the string with respect to the
    positive x-direction

## HPGe detector unit configuration

- `hpges` → HPGe detector name
  - `rodlength_in_mm` → length of the copper rods next to this detector. This is a
    "warm" length, i.e. it is multiplied by a factor < 1 to get the shorter rod
    length in the cryostat.
  - `baseplate` → size of the PEN plate below this detector (one value out of
    `small`, `medium`, `large`, `xlarge`)
