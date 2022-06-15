import sys
import re

comment_patt = re.compile(r'^\s*#')
class_patt = re.compile(r'^[ \t]*class ')
def_patt = re.compile(r'^[ \t]*def ')
fromimport_patt = re.compile(r'^([ \t]*)from\s*(metapathways[.]\S+)\s*import\s*(.*)$')

newnames = { 
    'metapathways.general_utils': 'gutils',
    'metapathways.metapathways_steps' : 'mpsteps',
    'metapathways.configuration' : 'configmod',
    'metapathways.pathwaytoolsutils' : 'ptoolsmod',
    'metapathways.diagnoze' : 'diagnoze',
    'metapathways.jobscreator': 'jobcreatormod',
    'metapathways.fastareader': 'fastreadrmod',
    'metapathways.sysutil':  'sysutils',
    'metapathways.LCAComputation' : 'lcacompute',
    'metapathways.parse': 'parsemod',
    'metapathways.metapathways_utils': 'mputils',
    'metapathways.execution': 'execmod',
    'metapathways.sequences': 'sequencemod',
    'metapathways.sampledata': 'sampledata',
    'metapathways.errorcodes': 'errormod',
    'metapathways.context': 'contextmod',
    'metapathways.MetaPathways_gather_run_stats': 'gatherstatmod',
    'metapathways.parameters': 'paramsmod'
}

function_module = {}

#with open(filename, 'r') as fin:
for line in sys.stdin:
    _line = line.rstrip()
    if comment_patt.search(line) or class_patt.search(line) or def_patt.search(line):      
        pass
        print(_line)
    else: # not a comment
        if fromimport_patt.search(_line): 
           res = fromimport_patt.search(_line)
           prefix = res.group(1) 
           module = res.group(2) 
           submodule = module.split('.')[1]
           functions = [x.strip() for x in res.group(3).strip().split(',') if x]
           for func in functions:
               function_module[func] = module
           print( prefix + 'from metapathways import '+ submodule + ' as ' + newnames[module])
        else: 
           for func in function_module.keys():
               res = re.search('[^.0-9a-zA-Z](' + func + ')[^.0-9a-zA-Z]?', _line)
                  
               if res:
                   replace =  res.group(1) 
                   _line = re.sub(replace, newnames[function_module[func]] + "." +  replace, _line)
           print(_line)

   
