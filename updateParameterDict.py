import ParameterDict
import sys

dictionaries = [item for item in dir(ParameterDict) if not item.startswith("__")]

# Loop over command-line arguments
for parameterChange in sys.argv[1:]:
    # parse this argument
    [name, parameter, value] = ((parameterChange.replace(":"," ")).replace("="," ")).split()
    
    # add this dictionary if not already included
    if not name in dictionaries:
         setattr(ParameterDict, name, name)
    
    # set parameter to specified value
    dictionary = getattr(ParameterDict, name)
    dictionary[parameter] = value

outFile = "updatedParameterDict.py"
with open(outFile, 'w') as f:
    for name in dictionaries:
        dictionary = getattr(ParameterDict, name)
        f.write('%s = {\n' % name)
        for key, value in dictionary.items():
            #f.write('\t\'%s\'\t\t:\t\t%s,\n' % (key, repr(value)))
            f.write('    %s:%s,\n'
                    % ('{:<44}'.format(repr(key)),
                       '{:>44}'.format(repr(value))
                      )
                   )
        f.write('}\n\n')