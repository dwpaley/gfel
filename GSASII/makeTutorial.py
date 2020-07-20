# -*- coding: utf-8 -*-
#makeTutorial.py
########### SVN repository information ###################
# $Date: 2019-09-02 08:29:27 -0700 (Mon, 02 Sep 2019) $
# $Author: toby $
# $Revision: 4131 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/makeTutorial.py $
# $Id: makeTutorial.py 4131 2019-09-02 15:29:27Z toby $
########### SVN repository information ###################
'''
*makeTutorial: Make Tutorial Web Page*
---------------------------------------------

Creates an HTML page (``GSASII/help/Tutorials.html``) listing all the tutorials in
:data:`GSASIIctrlGUI.tutorialIndex`. Run this after adding a new tutorial to that
index.
'''

from __future__ import print_function
import os
#import requests
import GSASIIpath
#import GSASIIctrl as G2G
onlineVideos = []
'''a list of videos that are in box, since I don't know how to check if they
are present anymore
'''
onlineVideos.append('https://anl.box.com/v/CalibrationofanareadetectorinG')
onlineVideos.append('https://anl.box.com/v/CalibrationTutorial')
onlineVideos.append('https://anl.box.com/v/CalibrationofaTOFpowderdiffrac')
onlineVideos.append('https://anl.box.com/v/Combinedrefinement')
onlineVideos.append('https://anl.box.com/v/TOFcombinedXNRietveldrefinemen')
onlineVideos.append('https://anl.box.com/v/NeutronCWPowderData')
onlineVideos.append('https://anl.box.com/v/FindProfParamCW')
onlineVideos.append('https://anl.box.com/v/DeterminingWavelength')
onlineVideos.append('https://anl.box.com/v/FitPeaks----')
onlineVideos.append('https://anl.box.com/v/LaboratoryX-')
onlineVideos.append('https://anl.box.com/v/FittingSmallAngleScatteringDat')
onlineVideos.append('https://anl.box.com/v/FitBkgTut---')
onlineVideos.append('https://anl.box.com/v/SmallAngleImageProcessing')
onlineVideos.append('https://anl.box.com/v/Integrationofareadetectordatai')
onlineVideos.append('https://anl.box.com/v/MerohedraltwinrefinementinGSAS')
onlineVideos.append('https://anl.box.com/v/ParametricFitting')
onlineVideos.append('https://anl.box.com/v/SequentialRefinementofSmallAng')
onlineVideos.append('https://anl.box.com/v/SequentialTutorial')
onlineVideos.append('https://anl.box.com/v/SimpleMagnetic')
onlineVideos.append('https://anl.box.com/v/SimTutorial-')
onlineVideos.append('https://anl.box.com/v/SmallAngleSizeDistribution')
onlineVideos.append('https://anl.box.com/v/StackingFaults-I')
onlineVideos.append('https://anl.box.com/v/StartingGSAS')
onlineVideos.append('https://anl.box.com/v/Strainfittingof2DdatainGSAS-II')
onlineVideos.append('https://anl.box.com/v/Textureanalysisof2DdatainGSAS-')
onlineVideos.append('https://anl.box.com/v/TOFSequentialSinglePeakFit')
#onlineVideos.append('

if __name__ == '__main__':
    GSASIIpath.SetBinaryPath()
    import GSASIIctrlGUI as G2G
    G2BaseURL = G2G.G2BaseURL
    tutorialIndex = G2G.tutorialIndex
    tutURL = G2BaseURL +'/Tutorials'
    outname = os.path.join(GSASIIpath.path2GSAS2,'help','Tutorials.html')

    dirList = [l[0] for l in tutorialIndex if len(l) >= 3]

    # loop through directories in Tutorials repository
    dirs = [d[:-1] for d in GSASIIpath.svnList(tutURL,False).split('\n') if d and d[-1] == '/']    
    for d in dirs:
        if d not in dirList: print(u"Tutorial directory not in GSASIIctrlGUI.tutorialIndex: "+d)

    #import sys
    #out = sys.stdout
    out = open(outname,'w')
    print('<!-- Do not edit this file. It is created by makeTutorial.py from info in GSASIIctrlGUI.py --!>',file=out)
    print('<h2>List of GSAS-II tutorials</H2><UL>',file=out)
    print('''
    <p> A list of available tutorials appears below. Each tutorial is a
    web page that can be opened using the link below, but most tutorials also need
    to have example data files downloaded. This can also be done with links included below, 
    but it can be easier to access tutorials using 
    <b>Help/Tutorials</b> menu item.
    When this menu entry is used from inside GSAS-II (unless "browse tutorial on web" is selected),
    the data files are downloaded to a local directory and GSAS-II will start from that directory
    for most file open commands. Most tutorials have also been recorded as videos of the computer screen
    along with naration. Links are provided below where videos are available. 
    </p>''',file=out)

    videolist = '<UL>'
    for l in tutorialIndex:
        if len(l) == 1:
            print("</UL><h4>{}</H4><UL>".format(l[0]),file=out)
        else:
            pageURL = tutURL+'/'+l[0]+'/'+l[1]
            dataURL = tutURL+'/'+l[0]+'/data'
            suffix = ''
            if l[2][0] == ' ':
                suffix = ' <A href="#prereq">*</A>'
            if suffix: 
                print('<UL><LI><A href="{}">{}</A>{}'.format(pageURL,l[2].strip(),suffix),file=out)
            else:
                print('<LI><A href="{}">{}</A>'.format(pageURL,l[2].strip()),file=out)
            
            # check for video tutorial
            videoName = '{:-<12s}'.format(
                os.path.splitext(l[1])[0].replace(' ','')[:30])            
            vname = 'https://anl.box.com/v/{}'.format(videoName)
            #if requests.get(vname).status_code == 200:
            if vname in onlineVideos:
                video = '<A href="{}">video</A>'.format(vname)
                #print(' [link: <A href="{}">video</A>]'.format(vname),file=out)
                #print('Found video',vname)
                videolist += '<LI><A href="{}">{}</A></LI>\n'.format(vname,l[2].strip())
            else:
                video =''
                print('No video for {:45s}{}'.format(videoName,l[2]))
            # check for data
            if GSASIIpath.svnList(dataURL,False):
                exampledata = '<A href="{}">Exercise files</A>'.format(dataURL)
                #print(' [link: <A href="{}">Exercise files</A>].'.format(dataURL),file=out)
            else:
                exampledata = ''
                #print(' [No exercise files].',file=out)
            if video and exampledata:
                print(' [links: {}, {}].'.format(video, exampledata),file=out)
            elif exampledata:
                print(' [link: {}].'.format(exampledata),file=out)
            elif video:
                print(' [link: {}, no example data].'.format(video),file=out)
            else:
                print(' [no example data or video].',file=out)
                
            if len(l) > 3:
                print("<blockquote><I>"+l[3]+"</I></blockquote>",file=out)
            if suffix: print('</UL>',file=out)
    #        if l[2][0] == ' ':
    #            print(' (Note that this tutorial requires previous as prerequisite)',file=out)

    videolist += '</UL>\n'
    print('</UL>\n<A name=prereq>* Indented tutorials require the previous unindented tutorial as a prerequisite',file=out)
    print('<h3>Tutorials with video-recorded examples</H3>', file=out)
    print(videolist, file=out)
    print("<P>The video tutorials are also <A href=https://pan.baidu.com/s/1C1jq1amfuVmcY2n91cQcsg> mirrored in China</A></P>",
                  file=out)
    out.close()
