# Units and models in coupled systems

The coupled stack has two entities that must be understood clearly.

## Units

**Units** are the elements that are coupled together. They are the "slots" in the joint state space: unit 1, unit 2, unit 3, etc. The number of units is the length of the first element of `coupling`.

## Models

**Models** are specs that define a mode: they are given by the model-spec tuples such as `G`, `R`, `S`, `transitions`, `insertstep`, etc. Each model has an index (1, 2, 3, …) corresponding to the position in those tuples. **Two units can use the same model.**

## The unit–model map: `coupling[1]`

The first element of `coupling` is the **unit–model map** (often called `unit_model`). It defines which model each unit uses:

- **`[1, 2, 3]`**  
  Unit 1 uses model 1, unit 2 uses model 2, unit 3 uses model 3.

- **`[3, 2, 1]`**  
  Unit 1 uses model 3, unit 2 uses model 2, unit 3 uses model 1.

- **`[1, 2, 2]`**  
  Unit 1 uses model 1, unit 2 uses model 2, unit 3 uses model 2 (two units share the same model).

So `unit_model[i]` is the **model index** for **unit** `i`. When building matrix components (e.g. in `set_elements_*`), you must loop over **units** and use `unit_model[unit]` to get the **model** for that unit, then use that model’s specs (`G[model]`, `R[model]`, etc.) and that model’s rate offset in the flat rate vector.

## Flat rate vector order

The flattened rate vector is in **model order**, followed by coupling constants:

- Block 1: rates for model 1  
- Block 2: rates for model 2  
- …  
- Then: coupling strength parameters  

So rate offsets are **per model**, not per unit. When assigning element rate indices for unit `i` (which uses model `unit_model[i]`), the base index is the offset for that **model**. The component struct (e.g. `TCoupledFullComponents`) stores only its size and the elements; each element’s rate index already points into this flat layout. Code that builds the matrix (e.g. `make_mat_T`) takes the component and the flat rate vector and constructs the matrix.

## What needs what

- **Functions that build matrix components** (e.g. from `coupling`, model specs, rates, noiseparams) need both unit and model information: loop over units, use `unit_model[unit]` to get the model, then use that model’s specs and the model’s rate offset for indices.
- **The component struct itself** only needs its size (e.g. `N`) and the list of elements; it does not store `unit_model` or rate offsets.
- **Building the matrix** is then `make_mat_T(components, flat_rates)` where `flat_rates` is the full flat vector in model order followed by coupling constants.
