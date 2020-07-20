#!/usr/bin/env python
# Installs GSAS-II from network using subversion and creates platform-specific shortcuts.
# works for Mac & Linux & Windows
import os, stat, sys, platform, subprocess, datetime

version = "$Id: bootstrap.py 3515 2018-07-30 02:14:14Z toby $"
g2home = 'https://subversion.xray.aps.anl.gov/pyGSAS/'
path2GSAS2 = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
now = str(datetime.datetime.now())
print('Running bootstrap from {} at {}\n\tId: {}'.format(path2GSAS2,now,version))
fp = open(os.path.join(path2GSAS2,'bootstrap.log'),'a')
fp.write('Running bootstrap from {} at {}\n\tId: {}\n'.format(path2GSAS2,now,version))
fp.close()
################################################################################
################################################################################
def GetConfigValue(*args): return True
# routines copied from GSASIIpath.py
proxycmds = []
'Used to hold proxy information for subversion, set if needed in whichsvn'
svnLocCache = None
'Cached location of svn to avoid multiple searches for it'

def MakeByte2str(arg):
    '''Convert output from subprocess pipes (bytes) to str (unicode) in Python 3.
    In Python 2: Leaves output alone (already str). 
    Leaves stuff of other types alone (including unicode in Py2)
    Works recursively for string-like stuff in nested loops and tuples.

    typical use::

        out = MakeByte2str(out)

    or::

        out,err = MakeByte2str(s.communicate())
    
    '''
    if isinstance(arg,str): return arg
    if isinstance(arg,bytes): return arg.decode()
    if isinstance(arg,list):
        return [MakeByte2str(i) for i in arg]
    if isinstance(arg,tuple):
        return tuple([MakeByte2str(i) for i in arg])
    return arg

def getsvnProxy():
    '''Loads a proxy for subversion from the file created by bootstrap.py
    '''
    proxyinfo = os.path.join(path2GSAS2,"proxyinfo.txt")
    if os.path.exists(proxyinfo):
        global proxycmds
        global host,port  # only in bootstrap.py
        proxycmds = []
        fp = open(proxyinfo,'r')
        host = fp.readline().strip()
        port = fp.readline().strip()
        fp.close()
        setsvnProxy(host,port)
        if not host.strip(): return '',''
        return host,port
    return '',''

def setsvnProxy(host,port):
    '''Sets the svn commands needed to use a proxy
    '''
    global proxycmds
    proxycmds = []
    host = host.strip()
    port = port.strip()
    if not host.strip(): return
    proxycmds.append('--config-option')
    proxycmds.append('servers:global:http-proxy-host='+host)
    if port.strip():
        proxycmds.append('--config-option')
        proxycmds.append('servers:global:http-proxy-port='+port)
        
def whichsvn():
    '''Returns a path to the subversion exe file, if any is found.
    Searches the current path after adding likely places where GSAS-II
    might install svn. 

    :returns: None if svn is not found or an absolute path to the subversion
      executable file.
    '''
    # use a previosuly cached svn location
    global svnLocCache
    if svnLocCache: return svnLocCache
    # prepare to find svn
    is_exe = lambda fpath: os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    svnprog = 'svn'
    if sys.platform.startswith('win'): svnprog += '.exe'
    host,port = getsvnProxy()
    if GetConfigValue('debug') and host:
        print('DBG_Using proxy host {} port {}'.format(host,port))
    # add likely places to find subversion when installed with GSAS-II
    pathlist = os.environ["PATH"].split(os.pathsep)
    pathlist.insert(0,os.path.split(sys.executable)[0])
    pathlist.insert(1,path2GSAS2)
    for rpt in ('..','bin'),('..','Library','bin'),('svn','bin'),('svn',),('.'):
        pt = os.path.normpath(os.path.join(path2GSAS2,*rpt))
        if os.path.exists(pt):
            pathlist.insert(0,pt)    
    # search path for svn or svn.exe
    for path in pathlist:
        exe_file = os.path.join(path, svnprog)
        if is_exe(exe_file):
            try:
                p = subprocess.Popen([exe_file,'help'],stdout=subprocess.PIPE)
                res = p.stdout.read()
                p.communicate()
                svnLocCache = os.path.abspath(exe_file)
                return svnLocCache
            except:
                pass        
    svnLocCache = None

def svnVersion(svn=None):
    '''Get the version number of the current subversion executable

    :returns: a string with a version number such as "1.6.6" or None if
      subversion is not found.

    '''
    if not svn: svn = whichsvn()
    if not svn: return

    cmd = [svn,'--version','--quiet']
    s = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = MakeByte2str(s.communicate())
    if err:
        print ('subversion error!\nout=%s'%out)
        print ('err=%s'%err)
        s = '\nsvn command:  '
        for i in cmd: s += i + ' '
        print(s)
        return None
    return out.strip()

def svnVersionNumber(svn=None):
    '''Get the version number of the current subversion executable

    :returns: a fractional version number such as 1.6 or None if
      subversion is not found.

    '''
    ver = svnVersion(svn)
    if not ver: return 
    M,m = ver.split('.')[:2]
    return int(M)+int(m)/10.

################################################################################
################################################################################
print(70*'*')
#testing for incorrect locale code'
try:
    import locale
    locale.getdefaultlocale()
except ValueError:
    print('Your location is not set properly. This causes problems for matplotlib')
    print('  (see https://github.com/matplotlib/matplotlib/issues/5420.)')
    print('Will try to bypass problem by setting LC_ALL to en_US.UTF-8 (US English)')
    os.environ['LC_ALL'] = 'en_US.UTF-8'
    locale.getdefaultlocale()
print('Preloading matplotlib to build fonts...')
try:
    import matplotlib
except:
    pass
print('Checking python packages...')
missing = []
for pkg in ['numpy','scipy','matplotlib','wx',]:
    try:
        exec('import '+pkg)
    except:
        missing.append(pkg)

if missing:
    msg = """Sorry, this version of Python cannot be used
for GSAS-II. It is missing the following package(s):
\t"""
    for pkg in missing: msg += " "+pkg
    print(msg)
    print("\nPlease install these package(s) and try running this again.")
    print("Showing first error: ")
    for pkg in ['numpy','scipy','matplotlib','wx',]:
        exec('import '+pkg)
    sys.exit()
try:
    import OpenGL
    install_with_easyinstall = None 
except:
    try:  	 	 
	    from setuptools.command import easy_install  	 	 
    except ImportError:
        print('You are missing the OpenGL Python package. This can be ')
        print('installed by this script if the setuptools package is installed')
        print('Please install either OpenGL (pyopengl) or setuptools')
        print("package and try running this again.")
        sys.exit()
    print("Missing the OpenGL Python package. Will attempt to install this later.")
    def install_with_easyinstall(package):  	 	 
        try:   	 	 
            print("trying system-wide ")  	 	 
            easy_install.main(['-f',os.path.split(__file__)[0],package])  	 	 
            return  	 	 
        except:  	 	 
            pass  	 	 
        try:   	 	 
            print("trying user level ") 	 	 
            easy_install.main(['-f',os.path.split(__file__)[0],'--user',package])  	 	 
            return  	 	 
        except:  	 	 
            print("\nInstall of "+package+" failed. Error traceback follows:")
            import traceback  	 	 
            print(traceback.format_exc())
            sys.exit() 

host = None
port = '80'
print('\nChecking for subversion...')
svn = whichsvn() # resets host & port if proxyinfo.txt is found
if not svn:
    print("Sorry, subversion (svn) could not be found on your system.")
    print("Please install this or place in path and rerun this.")
    #raise Exception('Subversion (svn) not found')
    sys.exit()
print(' found svn image: '+svn)

if install_with_easyinstall:  	 	 
    print('\nInstalling PyOpenGL. Lots of warnings will follow... ')
    install_with_easyinstall('PyOpenGl')  	 	 
    print('done.')
    
print('Ready to bootstrap GSAS-II from repository\n\t'+g2home+'\nto '+path2GSAS2)
proxycmds = []
proxyinfo = os.path.join(path2GSAS2,"proxyinfo.txt")
if os.path.exists(proxyinfo):
    fp = open(proxyinfo,'r')
    host = fp.readline().strip()
    port = fp.readline().strip()
    fp.close()
    os.remove(proxyinfo)
print(70*"=")
if sys.version_info[0] == 2:
    getinput = raw_input
else:
    getinput = input

# get proxy from environment variable
key = None
for i in os.environ.keys():
    if 'https_proxy' == i.lower():
        key = i
        break
else:
    for i in os.environ.keys():
        if 'http_proxy' == i.lower():
            key = i
            break
if key:
    val = os.environ[key].strip()
    if val[-1] == '/':
        val = val[:-1]
    if len(val.split(':')) > 2:
        host = ':'.join(val.split(':')[:-1])
        port = val.split(':')[-1]
    else:
        host = ':'.join(val.split(':')[:-1])
        port = val.split(':')[-1]

# get proxy from user, if terminal available
try:
    if host:
        ans = getinput("Enter the proxy address (type none to remove) ["+host+"]: ").strip()
        if ans.lower() == "none": host = ""
    else:
        ans = getinput("Enter your proxy address [none needed]: ").strip()
        if ans: host = ans
    if host:
        ans = getinput("Enter the proxy port ["+port+"]: ").strip()
        if ans == "": ans=port
        port = ans
except EOFError:
    host = ""
    port = ""
if host:
    proxycmds.append('--config-option')
    proxycmds.append('servers:global:http-proxy-host='+host.strip())
    if port:
        proxycmds.append('--config-option')
        proxycmds.append('servers:global:http-proxy-port='+port.strip())
    fp = open(proxyinfo,'w')
    fp.write(host.strip()+'\n')
    fp.write(port.strip()+'\n')
    fp.close()
    fp = open(os.path.join(path2GSAS2,'bootstrap.log'),'a')
    fp.write('Proxy info written: {} port\n'.format(host,port))
    print('Proxy info written: {} port'.format(host,port))
    fp.close()

# patch: switch GSAS-II location if linked to XOR server (removed May/June 2017)
cmd = [svn, 'info']
p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
res,err = p.communicate()
if '.xor.' in str(res):
    print('Switching previous install with .xor. download location to\n\thttps://subversion.xray.aps.anl.gov/pyGSAS')
    cmd = [svn, 'switch','--relocate','https://subversion.xor.aps.anl.gov/pyGSAS',
           'https://subversion.xray.aps.anl.gov/pyGSAS']
    if proxycmds: cmd += proxycmds
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res,err = p.communicate()
    if err:
        print('Please report this error to toby@anl.gov:')
        print(err)
        print(res)
# patch: switch GSAS-II location if switched to 2frame version (removed August 2017)
if '2frame' in str(res):
    print('Switching previous 2frame install to trunk\n\thttps://subversion.xray.aps.anl.gov/pyGSAS')
    cmd = [svn, 'switch',g2home + '/trunk',path2GSAS2,
           '--non-interactive','--trust-server-cert',
           '--accept','theirs-conflict','--force','--ignore-ancestry']
    if proxycmds: cmd += proxycmds
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    res,err = p.communicate()
    if err:
        print('Please report this error to toby@anl.gov:')
        print(err)
        print(res)

print('\n'+75*'*')
print('Now preparing to install GSAS-II')
tryagain = True
err = False
firstPass = 0
while(tryagain):
    tryagain = False
    if err:
        print('Retrying after a cleanup...')
        cmd = [svn, 'cleanup', path2GSAS2]
        s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
        out,err = s.communicate()    
        if err:
            print('subversion returned an error:')
            print(out)
            print(err)
    cmd = [svn, 'co', g2home+ 'trunk/', path2GSAS2, '--non-interactive', '--trust-server-cert']
    if proxycmds: cmd += proxycmds
    msg = 'svn load command: '
    for item in cmd: msg += " "+item
    print(msg)
    s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
    print('\nsubversion output:')
    out,err = s.communicate()
    if err:
        print('subversion returned an error:')
        print(out)
        print(err)
        if firstPass == 0: tryagain = True
    firstPass += 1
if err:
    print('Retrying with a command for older svn version...')
    cmd = [svn, 'co', g2home+ 'trunk/', path2GSAS2]
    if proxycmds: cmd += proxycmds
    msg = ""
    for item in cmd: msg += " " + item
    print(msg)
    s = subprocess.Popen(cmd,stderr=subprocess.PIPE)
    out,err = s.communicate()
    if err:
        print('subversion returned an error:')
        print(err)
        print('  *** GSAS-II failed to be installed. A likely reason is a network access')
        print('  *** problem, most commonly because you need to use a network proxy. Please')
        print('  *** check with a network administrator or use http://www.whatismyproxy.com/\n')
        sys.exit()
print('\n'+75*'*')

try:
    import GSASIIpath
    print('import of GSASIIpath completed')
except Exception as err:
    print('\n'+75*'=')
    print('Failed with import of GSASIIpath. This is unexpected.')
    print('GSAS-II will not run without correcting this. Contact toby@anl.gov')
    print(err)
    print(75*'=')
    sys.exit()

for a in sys.argv[1:]:
    if 'allbin' in a.lower() or 'server' in a.lower():
        print('Loading all binaries with command...')
        if not GSASIIpath.svnSwitchDir('AllBinaries','',g2home+ 'Binaries/',None,True):
            print('Binary load failed')
            sys.exit()
        break
else:
    GSASIIpath.DownloadG2Binaries(g2home)
        
#===========================================================================
# test if the compiled files load correctly
#===========================================================================

script = """  # commands that test each module can at least be loaded & run something in pyspg
try:
    import GSASIIpath
    GSASIIpath.SetBinaryPath()
    import pyspg
    import histogram2d
    import polymask
    import pypowder
    import pytexture
    pyspg.sgforpy('P -1')
    print('==OK==')
except:
    pass
"""
p = subprocess.Popen([sys.executable,'-c',script],stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                     cwd=path2GSAS2)
res,err = p.communicate()
if '==OK==' not in str(res) or p.returncode != 0:
    print('\n'+75*'=')
    print('Failed when testing the GSAS-II compiled files. GSAS-II will not run')
    print('without correcting this.\n\nError message:')
    #print(res)
    print(err)
    #print('\nAttempting to open a web page on compiling GSAS-II...')
    print('Please see web page\nhttps://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/CompileGSASII')
    #import webbrowser
    #webbrowser.open_new('https://subversion.xray.aps.anl.gov/trac/pyGSAS/wiki/CompileGSASII')
    print(75*'=')
#    if '86' in platform.machine() and (sys.platform.startswith('linux')
#                                        or sys.platform == "darwin"
#                                        or sys.platform.startswith('win')):
#        print('Platform '+sys.platform+' with processor type '+platform.machine()+' is supported')
#    else:
#        print('Platform '+sys.platform+' with processor type '+platform.machine()+' not is supported')
    sys.exit()
else:
    print('Successfully tested compiled routines')
#===========================================================================
# import all .py files so that .pyc files get created
print('Byte-compiling all .py files...')
import compileall
compileall.compile_dir(path2GSAS2,quiet=True)
print('done')
#===========================================================================
# platform-dependent stuff
#===========================================================================
if sys.version_info[0] > 2:
    def execfile(file):
        with open(file) as source_file:
            exec(source_file.read())
# on Windows, make a batch file with Python and GSAS-II location hard-coded
if sys.platform.startswith('win') and os.path.exists(
    os.path.join(path2GSAS2,"makeBat.py")):
    execfile(os.path.join(path2GSAS2,"makeBat.py"))
#===========================================================================
# on a Mac, make an applescript 
elif sys.platform.startswith('darwin') and os.path.exists(
    os.path.join(path2GSAS2,"makeMacApp.py")):
    sys.argv = [os.path.join(path2GSAS2,"makeMacApp.py")]
    print(u'running '+sys.argv[0])
    execfile(sys.argv[0])
#===========================================================================
# On linux, make desktop icon
elif sys.platform.startswith('linux') and os.path.exists(
    os.path.join(path2GSAS2,"makeLinux.py")):
    sys.argv = [os.path.join(path2GSAS2,"makeLinux.py")]
    print(u'running '+sys.argv[0])
    execfile(sys.argv[0])

