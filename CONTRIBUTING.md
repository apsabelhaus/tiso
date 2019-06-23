# Contributing

A quick note on the style convention in the tiso repository.
We use the MATLAB style guideline at:

https://www.mathworks.com/matlabcentral/fileexchange/46056-matlab-style-guidelines-2-0

Here are the basic conventions on naming:

#### Naming Conventions in tiso

1. Folder names: lower case with underscores. For example, /is_2d/.
2. Functions: either all lower case or Camel Case, with underscores used if needed to separate text for clarity or after an acronym. For example, getB_3d.m
3. Single letters, special examples: lower case. For example, 'dimension' = d, as in 2d or 3d.
4. File names: follow function names, camel case, unless previously used in matlab as all lowercase.
5. Variable names: follow function names, camel case, unless previously used in matlab as all lowercase.
6. Plural: no plural ever. For example, 'inputs' becomes 'input'. When necessary, words are shortened to make this easier to read: 'coordinates' becomes 'coord' so confusion doesn't exist with the 's' ending.
7. Struct names: when a variable is a struct, begin with a capital letter. For example, 'mandatory inputs' = MandInput.

### Style Conventions in tiso

1. Function, variable, file naming: Avoid words that are obvious in context. For example, 'node coordinates' is redundant since the only coordinates in the problem are for nodes (rigid body position vectors are states.)
2. Using 'disp': since matlab's disp function doesn't give a variable name, we commonly just use expressions without the ending semicolon for quick debugging.
3. 2D or 3D designation should be at the end of names of functions, etc.

#### Some abbreviations

Abbreviations used throughout:

- obj = objective function
- coord = coordinate(s)
- rxn = reaction
- fcn = function
- surf = surface
- pt = point (try to avoid the singular use of 'pt' without the 's', but sometimes necessary)
- pts = points
- discr = discretization
- rad = radius (in plotting functions), radians (sometimes with respect to state variables for bodies)
- traj = trajectory (-ies).

#### When to use 'is', 'iso', 'tiso'

In this library, we use the abbreviations:

1. is = inverse statics
2. iso = inverse statics optimization
3. tiso = tensegrity inverse statics optimization

These will be used depending on context and clarity.

- For example, when it's obvious we're working with tensegrity structures (i.e., internal to the library), use 'is' or 'iso.'
- External outputs should use 'tiso.'
- The 'o' ending versus simply 'is' used when it's unclear if the variable/function is related to the optimization program itself. For example, the root directories named 'is' because they address the inverse statics problem in general, while the 'parse' functions used 'iso' since they address the output of the optimization program.
- Avoid capitalizing these acronyms if possible, except where to do so increases readibility, for example in function names. The core optimization functions are named like 'rbISO_3d.m' for 'rigid body inverse statics optimization 3D.'

