# dimcheck: A physics dimension checker based on sympy

## Description
Dimcheck is a Python library that provides an interface for users to check, compare and manipulate dimensions of quantities. It is especially useful in scientific computing and physics, where ensuring correct dimensions is crucial.

## Dependencies
```
python3
sympy
numpy
```

## Installation

You can install the package using pip,
```
pip install dimcheck
```
## Basic usage

You can try the package as follows,
```python
# Import the package, the package requires some time to initialize the global instance "si". You can make the process faster by serialization in the later section.
>>> from dimcheck import si
Starting to parse the quantity definition file: 
/to/your/pip/lib/path/dimcheck/dimcheck/si.json
Successfully parsed!

# Get the dimension of quantity [E] where square bracket represents this is a quantity rather than basic unit like kg, m, etc. in the si unit system.
>>> si.dim("[E]")   
'kg*m**2/s**2'

# Get the dimension of combination quantities [m*v**2] which is same as [E].
>>> si.dim("[m*v**2]")  
'kg*m**2/s**2'

# Check whether two quantities have the same dimension.
>>> si.is_dc("[E]","[m*v**2]")  
True

# Derive the possible quantity of [G*m**2/r**2] which is the formula of Newton's law of universal gravitation. Here [G] is the gravitational Constant, [m] is the mass, [r] is the radius or length. You can use a space to replace the operator of the multiplication '*'.
>>> si.quant("[G m**2/r**2]")   
'[F]'

# Derive the possible quantity based on the Ohm's law. Here [V] is the voltage, [I] is the current, [R] is the resistance which have the same dimension as the von Klitzing constant [R_K].
>>> si.quant("[V/I]")
['[R]', '[R_K]']

# Restore the omitted quantity on the rhs. Here [hbar] is the reduced Planck constant, [k] is the wavenumber.
>>> si.omit_quant(lhs="[E]",rhs="[k]",omit_quant=["[hbar]","[v]"])  
'[E] = [k]*[hbar]*[v]'

# Derive the formula based on the given parameters. Here, we try to derive the pendulum formula.
>>> si.formula(lhs="[t]",parameters=["[l]","[g]"])   
'[t] = ([l]**0.5)*([g]**-0.5)'

# Display all quantities in the terminal.
>>> si.all_quant()
All Quantities:
[B]                 magnetic flux density, magnetic induction, B field           kg/(A*s**2)
[C]                 capacitance                                                  A**2*s**4/(kg*m**2)
[D]                 displacement field                                           A*s/m**2
...
None

# Save the definition of quantities in a ".csv" format file. By default, it will save to the "./Quantities.csv".
>>> si.save_all_quant()
True
```

## Convention
1. Unit: The basic units  in `"si"` unit system (which is the only unit system so far) are `{"kg", "m", "s", "A", "mol", "cd", "K"}`.

2. Quantity: In the package, this term is a combination of the number 1 and units , e.g. `[v]` is the unit velocity whose unit is `"m/s"`.

3. Expression: The combination of basic units e.g. `"kg*m**2/s**2"`.

4. In the package, it will follow the same convention to represent a quantity as Ref. "Fly by Night Physics". For a quantity or a combination of quantities you should wrap it with square brackets.
`"[F]","[sigma_3D]","[m*v**2]","[e**2/hbar]"`，
```python
>>> si.dim("[F]")
'kg*m/s**2'
>>> si.quant("[e**2/hbar]")
['[conductance]', '[sigma_2D]']
```

5. The basic units`{"kg", "m", "s", "A", "mol", "cd", "K"}` should not be wrapped with square brackets.
```python
>>> si.dim("m")
'm'
>>> si.quant("kg m/s")
'[p]'
```

6. You can replace `**` with `^` to represent the exponent. Also, the multiplication `*` can be replaced with a space ` `,
```python
>>> si.quant("[m v^2]")
'[E]'
```
There is also a few reserved keywords like `sqrt` and `cbrt` for convenience,
```python
>>> si.quant("[sqrt(hbar e B v^2)]")
'[E]'
```

## Advanced Usage

### Pretty printing
In the terminal mode, `Dimcheck` class (`si` is the instance of it) provides `is_pretty` property to make the results more easy to read.
```python
# Before setting "is_pretty" mode.
>>> si.is_pretty=False
>>> si.dim("[m a]")
'kg*m/s**2'

# After setting "is_pretty" mode.
>>> si.is_pretty=True
>>> si.dim("[m a]")
'kg*m/s²'

# Other methods
>>> si.is_dc("[sqrt(hbar e B v**2)]","[E]",is_print=True)
DC!
[sqrt(hbar e B v²)] == [E] (kg*m²/s²)
True

# Noted that "kg**(1/3)" will be shown as "kg¹ʴ³" since there is no other proper symbol found in the Unicode for the division that locates at the exponent.
>>> si.dim("[m]**(1/3)")
'kg¹ʴ³'
```

### Custom definition of the quantities
In `dimcheck`, you can define your own quantities and symbols by simply manipulate the `setting.json` file in the package directory. Here is the step to define your own unit system.

1. Find the position of the package:
```python
>>> from dimcheck import si
Starting to parse the quantity definition file: 
/to/your/pip/lib/path/dimcheck/dimcheck/si.json
Successfully parsed!
>>> si.setting_file
'/to/your/pip/lib/path/dimcheck/dimcheck/setting.json'
```

2. Copy the following setting to the `setting.json`. The `is_save` field is used to serialize the object.
```json
{
    "si_instance": true,
    "custom_instance": true,
    "custom_quant_def_file": "my_definition.json",
    "is_save": true,
    "output": "./"
}
```

3. Rename the file `'/to/your/pip/lib/path/dimcheck/dimcheck/custom.json'` to `'/to/your/pip/lib/path/dimcheck/dimcheck/my_definition.json'`

4. Add, delete or change the definition of symbols in the file. You need to be careful since the quantities should always be wrapped with a pair of square bracket `[]`. In the following example, the quantity `[epsilon_0]` is defined by the rhs `[Q]**2/([F]*[l]**2)`. Meanwhile, you can also give multiple aliases to the same quantity. Discription can be used to indicate the meaning of the quantity.
```json
{
    "quantity":"[epsilon_0]",
    "rhs":"[Q]**2/([F]*[l]**2)",
    "alias":["[eps_0]"],
    "discription":"vacuum dielectric constant"
}
```

5. After doing so, you can then import your definition by simply invoking
```python
# Import custom quantity definition. Here "cs" means custom
>>> from dimcheck import cs
Starting to parse the quantity definition file: 
/to/your/pip/lib/path/dimcheck/dimcheck/custom.json
Successfully parsed!

# Test your definition.
>>> cs.dim("[m v**2]")
'kg*m**2/s**2'
>>> cs.quant("[m v**2]")
'[E]'
```


### Serialization
Since sometimes loading `"si.json"` can be time-consuming, `dimcheck` provides the serialization to directly save the `Dimcheck` object (like `si`).
```python
# Serialize the instance
>>> si.save()
True

# The position of the file of quantity definition
>>> si.quant_def_file
'/to/your/pip/lib/path/dimcheck/dimcheck/si.json'

# The position of the serialized file
>>> si.serialized_file
'/to/your/pip/lib/path/dimcheck/dimcheck/si.pickle'

# Delete the serialized file
>>> si.clean()
True
```

### Save and display all quantities or expressions
`dimcheck` provides a few method to display and save the table of quantities and expressions. You can simply invoke the method,
```python
# Display all quantities in the terminal.
>>> si.all_quant()
All Quantities:
[B]                 magnetic flux density, magnetic induction, B field           kg/(A*s**2)
[C]                 capacitance                                                  A**2*s**4/(kg*m**2)
[D]                 displacement field                                           A*s/m**2
...

# Display all expressions with their aliases in the terminal.
>>> si.all_expr(is_inc_alias=True, is_sorted=True)
All Expression:
1                   [alpha] A['[alp]'] 
1/m                 [k] [nabla] A['[nbl]'] 
1/mol               [N_A]
...

# Make the output more easy to read
>>> si.is_pretty=True
>>> si.all_quant(is_sorted=False)
All Quantities:
[m]                  mass                           kg
[t]                  time                           s
...
[rho_m]              density of mass                kg/m³
[v]                  velocity                       m/s
[a]                  acceleration                   m/s²
[g]                  acceleration of gravity        m/s²
[F]                  force                          kg*m/s²

# Save quantities, by default the location will be "./Expressions.csv".
>>> si.save_all_expr()
True

# Save quantities to the path, by default the location will be "./Quantities.csv".
>>> si.save_all_quant(file_path="/path/to/output/Quantities.csv")
True
```
You can review above results from [Quantities.csv](./data/Quantities.csv) and [Expressions.csv](./data/Expressions.csv)
## Methods

### Property
| property    | Description           | Writable |
| :---:       | :----------------     | :---:|   
| is_pretty   | The output is setting to be more easy to read (based on some symbols in UTF-8) or not. | Yes|
| setting_file   | Returns the path of the global setting file.| No|   
| unit_system   | Returns the unit system.| No|  
| quant_def_file   | Returns the path of the file of quantity definition.| No|  
| serialized_file   | Returns the path of the serialized file.| No|  


### Dimension
| method    | Description           |  Parameters   |   Return   |
| :---:     | :----------------     |   :---:       |  :-----:  |
| dim       | To get the dimension of a quantity. | s: str| expr: str|
| dimension | Alias of dim.          | s: str        |expr: str|
| is_dc     | To check whether two quantities are same. | lhs: str, rhs: str, is_print: bool=False | dc: bool |
| quant | To print the possible quantity by deriving the combination of quantities. | s: str, is_print=False| quantity: str or quantities: list|

### Derive Formula
| method    | Description           |  Parameters   |   Return   |
| :---:     | :----------------     |   :---:       |  :-----:  |
| formula       | Print the formula of the target based on the parameters. | lhs : str, parameters : list| expr: str|
| omit_quant       | Restore the units of lhs based on the omitted quantities and rhs. | lhs : str, rhs : str, omit_quant: list| expr: str|




### Save, clean and display
| method    | Description           |  Parameters   |   Return |
| :---:     | :----------------     |   :---:       |:---: |
| all_expr | To display all expressions based on the quantity definition file. | is_inc_alias: bool=False, is_sorted: bool=True||
| all_quant | To display all quantities based on the quantity definition file. | is_inc_alias: bool=False, is_sorted: bool=True||
| clean | To clean the binary serialized file of the current object. | | success: bool|
| save  | To save the binary serialized file of the current object.  | | success: bool|
| save_all_expr | To save all expression in ".csv" format based on the quantity definition file. | file_path='./Expressions.csv', is_sorted: bool=True| success: bool|
| save_all_quant | To save all quantities in ".csv" format based on the quantity definition file. | file_path='./Quantities.csv', is_sorted: bool=True|success: bool|

