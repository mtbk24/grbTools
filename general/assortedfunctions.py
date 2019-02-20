from __future__ import division
import os, sys



def counted(f):
    '''
    counted(f)
    
    A counting decorator. Counts how many times a function is called.
    Place this at the front of a funciton and it will count the number of times it is called.
    
    @counted
    def pctDiff(val1, val2):
    
    After the loop is done, you can use pctDiff.calls to tell you how many times
    the pcdDiff function was called. This can be used when you want to plot data labels only once.
    
    '''
    def wrapped(*args, **kwargs):
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped



#######################################################################
old_stdout = sys.stdout  # SAVE ORIGINAL PRINT TO SCREEN SETTINGS.
# DISABLE PRINTING
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
# RESTORE PRINTING
def enablePrint():
    sys.stdout = old_stdout


#######################################################################
def nearest(array, value):
    ''' 
    nearest(array, value)
    
    Returns the index and the value in an array of floats that matches closest to an input value of interest.
    
    
    '''
    return min(enumerate(array), key=lambda x: abs(x[1]-value))


#######################################################################
def most_common(L):
    import itertools, operator
    # get an iterable of (item, iterable) pairs
    SL = sorted((x, i) for i, x in enumerate(L))
    # print 'SL:', SL
    groups = itertools.groupby(SL, key=operator.itemgetter(0))
    # auxiliary function to get "quality" for an item
    def _auxfun(g):
        item, iterable = g
        count = 0
        min_index = len(L)
        for _, where in iterable:
            count += 1
            min_index = min(min_index, where)
        # print 'item %r, count %r, minind %r' % (item, count, min_index)
        return count, -min_index
    # pick the highest-count/earliest item
    return max(groups, key=_auxfun)[0]


#######################################################################
def _fileExists(filename):
    '''
    Check to see if file exists.
    '''
    if os.path.exists(filename):
        raise Exception, "File Already Exists and I won't overwrite it.  Choose another file name."
    else:
        print("File does not exist.  Safe to create it. ")



#######################################################################
def _fileExistsIterate(filename):
    '''
    This will check to see if the file exists, and if so, it takes the most recent version number and adds one onto it to proivde a new filename.
    For example, if you pass /Users/KimiZ/GRBs/fakeit_results-01-grbm.fit and say all versions up to and including -04- already exist, this will produce a new file /Users/KimiZ/GRBs/fakeit_results-05-grbm.fit
    
    _fileExistsIterate(filename)
    
    filename: str, full filename including directory if you are outside the directory.  If not, pass only the filename. 
    
    The base of the filename must have a version number separated by -'s.  For example:  fakeit_results-01-grbm.fit
        Returns a new filename with the a version number 1 greater than the highest that already exists.  This is to ensure no files are written over.
    '''
    direc = os.path.dirname(filename)
    filebase = os.path.basename(filename)
    parts = filebase.split("-")[0], filebase.split("-")[-1]
    if parts[0] == parts[1]:
        raise IOError("need to have '-'s' in filename, separating the version number from the rest of the filename. ex: file-01-grbm.txt ")
    else:
        pass

    for i in range(1, 100):
        if not os.path.exists(direc+"/"+"%s-%02i-%s"%(parts[0], i, parts[1])):
            filename_new = os.path.join(direc, "%s-%02i-%s"%(parts[0], i, parts[1]))
            return filename_new
        else:
            pass


#######################################################################
def _findSimilarFiles(filename):
    '''
    _findSimilarFiles(filename)
    
    filename: str, full filename including directory if you are outside the directory.  If not, pass only the filename. 
    
    The base of the filename must have a version number separated by -'s.  For example:  fakeit_results-01-grbm.fit
        All files with *fakeit_results* and *grbm.fit* will be returned. Returns a list of all versions of the files that have the same start and end, separated by the -01- which represents the version number.
    '''
    direc = os.path.dirname(filename)
    filebase = os.path.basename(filename)
    parts = filebase.split("-")[0], filebase.split("-")[-1]
    if parts[0] == parts[1]:
        raise IOError("need to have '-' in filename to represent version. ex: file-01-grbm.txt ")
    else:
        pass
    possible_files = glob.glob(os.path.join(direc,"*%s*%s*" %(filebase.split("-")[0], filebase.split("-")[-1])))
    return possible_files


#######################################################################
# def _checkExistsCreate(directory, answer=None):
#     '''
#     _checkExistsCreate(directory)
    
#     directory:  str, pass directory name that you want to create if it does not already exist. 
#     Checks existence first so the origianl is not overwritten.
    
#     answer: str, "y" "ye" or "yes" to make directory without being prompted.  Use "no" "n" or any other response to contine without making the directory. 
#     Use this feature, answer = "yes" or answer = "no", if you DO NOT want the program to prompt you for a response. 
#     Default is None, so if answer is left as None, the program will prompt you for a response.
    
#     '''
#     if not os.path.exists(directory):
#         if answer is not None:
#             response = answer
#         else:
#             response = raw_input("Directory %s does not exist, would you like to create it now?  y or n:  "%directory)
        
#         if response in list(["y", "yes", "ye", "Y", "YE", "YES", "Ye", "Yes"]):
#             os.mkdir(directory)
#             print "\n"*3
#             print "***  Directory %s has been created."%directory
#             print "\n"*3
#         else:
#             print "\n"*3
#             print "***  Will not create directory."
#     else:
#         pass


#######################################################################
def _checkExistsCreate(directory, autoCreate=False):
    '''
    _checkExistsCreate(directory, autoCreate=False)
    
    directory:  str, pass directory name that you want to create.
    autoCreate: True or False. True indicates the directory will automatically
                be created without asking the user. 
                Default is False.
    
    '''
    response1 = "\n Directory %s created. \n"%directory
    response2 = "\n Directory %s was NOT created. \n"%directory

    if not os.path.exists(directory):
        if autoCreate:
            # if autoCreate = True, automatically create directory without asking.
            os.mkdir(directory)
            print(response1)
        else:
            resp = raw_input("\n Directory %s does not exist, would you like to create it now? \n y or n: \n"%directory)
            if resp in list(["y", "yes", "ye", "Y", "YE", "YES", "Ye", "Yes"]):
                os.mkdir(directory)
                print(response1)
            else:
                print(response2)
                pass
    else:
        print("\n Directory %s already exists. Passing... \n"%directory)
        pass



#######################################################################
def _checkExistsCreate2(directory, answer=None):
    '''
    _checkExistsCreate(directory)
    
    directory:  str, pass directory name that you want to create if it does not already exist. 
    Checks existence first so the origianl is not overwritten.
    
    answer: answer is always yes
    
    '''
    if not os.path.exists(directory):
        if answer is not None:
            response = answer
        else:
            response = raw_input("Directory %s does not exist, would you like to create it now?  y or n:  "%directory)
        
        if response in list(["y", "yes", "ye", "Y", "YE", "YES", "Ye", "Yes"]):
            os.mkdir(directory)
            print "\n"*3
            print "***  Directory %s has been created."%directory
            print "\n"*3
        else:
            print "\n"*3
            print "***  Will not create directory."
    else:
        pass


#######################################################################
def _getAbsPath(filename):
  return os.path.abspath(os.path.expanduser(filename))
pass


#######################################################################
def _getIterable(obj):
  '''
  If obj is not iterable, this method encapsulate it in a list and return such a list
  '''
  try:
      dumb                    = len(obj)
      return obj
  except:
      inList                  = [obj]
      return inList  
  pass
pass


#######################################################################


#def printdocs(f):
#    print("%s: %s"%(f.__name__, f.__doc__))
#
#@printdocs
#def get_detdir(det):
#    ''' Enter 'L' or 'G' for detector analysis, and 'GBM' or 'GBMwLAT' returned.'''
#    if 'L' in det:
#        det_dir = "GBMwLAT"
#    elif 'G' in det:
#        det_dir = "GBM"
#    else:
#        raise Exception, "Don't recognize the detector"
#    return det_dir

