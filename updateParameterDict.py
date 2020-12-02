import ParameterDict
import sys

def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

dictionaries = [item for item in dir(ParameterDict) if not item.startswith("__")]

# Loop over command-line arguments
for parameterChange in sys.argv[2:]:
    # parse this argument
    [name, parameter, value] = ((parameterChange.replace(":"," ")).replace("="," ")).split()
    #if value=="False" or value=="false": value=False
    #elif value=="True" or value=="true": value=True
    
    value = autoconvert(value)
    
    # add this dictionary if not already included
    if not name in dictionaries:
        dictionary = {parameter : value}
        setattr(ParameterDict, name, dictionary)
    else:
        # set existing parameter to specified value
        dictionary = getattr(ParameterDict, name)
        dictionary[parameter] = value

# Print out updated ParameterDict
outFile = sys.argv[1] #.replace(".","_").replace("%","").replace("-","_")
with open(outFile, 'w') as f:
    dictionaries = [item for item in dir(ParameterDict) if not item.startswith("__")]
    for name in dictionaries:
        dictionary = getattr(ParameterDict, name)
        f.write('%s = {\n' % name)
        for key, value in dictionary.items():
            f.write('    %s:%s,\n'
                   % ('{:<44}'.format(repr(key)),
                   '{:>44}'.format(repr(value))
                     )
                   )
        f.write('}\n\n')