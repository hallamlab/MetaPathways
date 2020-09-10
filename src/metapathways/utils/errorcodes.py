import sys
errorcodes={}
errorcodes[0]="Run completed successfully"
errorcodes[1]="Error during the PREPROCESSING step"
errorcodes[2]="Error during the ORF PREDICTION step"
errorcodes[3]="Error during the CREATE and FILTER AMINOS step"
errorcodes[4]="Error during the FUNC_SEARCH  step"
errorcodes[5]="Error during the PARSE_FUNC_SEARCH step"
errorcodes[6]="Error during the SCAN_rRNA step"
errorcodes[7]="Error during the SCAN_tRNA step"
errorcodes[8]="Error during the ANNOTATION step"
errorcodes[9]="Error during the BUILD_PGDB step"
errorcodes[10]="Error during the COMPUTE_RPKM step"
errorcodes[200]="Multiple errors found run again with \"-v\" or check the \"errors_warnings_log.txt\" for details"
errorcodes[11]="Error while creating biom format files"
errorcodes[12]="Error while creating Pathologic input"
errorcodes[11]="Error while creating biom format files"
errorcodes[15]="Error while computing refscores"
errorcodes[16]="Error during the CREATE_ANNOT_REPORTS step"
errorcodes[17]="Raw secuences for databases missing"
errorcodes[18]="Missing environment variable"

errorcodes[200]="Multiple errors found run again with \"-v\" or check the \"errors_warnings_log.txt\" for details"

errors={}
error_list=[]


def insert_error(i):
   global errors
   global errorcodes

   if i in errorcodes:
      errors[i] = errorcodes[i]
   else:
      errors[i] = "Unknown error"

   error_list.append(i)

def error_message(i) :
    global errorcodes
    if i in errorcodes:
       return errorcodes[i]
    return "Unkonwn error"

def get_error_list():
   global errors
   return errors

def get_recent_error():
   global errors
   if error_list:
      return error_list[-1]
   return 0

def return_code():
   global errors
   if errors:
      sys.exit(errors.keys()[0])
   sys.exit(0)
   
def exit_code(code):
   sys.exit(code)
   
